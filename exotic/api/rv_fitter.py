# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #
from astropy import constants as const
from astropy import units as u
import copy
from itertools import cycle
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ultranest import ReactiveNestedSampler
from scipy.optimize import least_squares

try:
    from elca import lc_fitter
except ImportError:
    from .elca import lc_fitter


Mjup = const.M_jup.to(u.kg).value
Msun = const.M_sun.to(u.kg).value
Rsun = const.R_sun.to(u.m).value
Grav = const.G.to(u.m**3/u.kg/u.day**2).value

def planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0, mu=1):
    # please see original: https://github.com/ucl-exoplanets/pylightcurve/blob/master/pylightcurve/models/exoplanet_lc.py
    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * np.pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [mu*bb * np.sin(inclination), mu*sma_over_rs * np.sin(vv), - bb * mu*np.cos(inclination)]

    if periastron < np.pi / 2:
        aa = 1.0 * np.pi / 2 - periastron
    else:
        aa = 5.0 * np.pi / 2 - periastron
    bb = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * np.pi
    mid_time = float(mid_time) - (period / 2.0 / np.pi) * (bb - eccentricity * np.sin(bb))
    m = (time_array - mid_time - np.int_((time_array - mid_time) / period) * period) * 2.0 * np.pi / period
    u0 = m
    stop = False
    u1 = 0
    for ii in range(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (1 - eccentricity * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-6)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise RuntimeError('Failed to find a solution in 10000 loops')

    vv = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan((u1) / 2))
    rr = mu*sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv))

    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)

    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))

    return [x, y, z]

def rv_model(time, params, dt=0.0001):
    """
    Compute the radial velocity model for a planet orbiting a star

    Args:
        time (array): time array [BJD]
        params (dict): dictionary of parameters
        dt (float): time step for computing velocity

    Returns:
        array: radial velocity model [m/s]    
    """
    xp,yp,zp = planet_orbit(params['per'], params['ars'], params['ecc'], 
                            params['inc'], params['omega'], params['tmid'], 
                            time, mu=1-params['mu'], ww=0)
    # TODO optimize by only computing in X direction
    xp2,yp2,zp2 = planet_orbit(params['per'], params['ars'], params['ecc'],
                            params['inc'], params['omega'], params['tmid'],
                            time+dt, mu=1-params['mu'], ww=0)
    return (xp2-xp)*params['mu']*params['rstar']*80520833.33333333 # Rstar/day -> m/s

def acceleration_model(time, params, dt=0.0001):
    # estimate radial acceleration from velocity
    v1 = rv_model(time, params, dt=dt)
    v2 = rv_model(time+dt, params, dt=dt)
    return (v2-v1)/dt # m/s^2

def timetrans_to_timeperi(tc, per, ecc, omega):
    """
    Convert Time of Transit to Time of Periastron Passage - from RADVEL

    Args:
        tc (float): time of transit    
        per (float): period [days]
        ecc (float): eccecntricity
        omega (float): longitude of periastron (radians)
    
    Returns:
        float: time of periastron passage

    """
    try:
        if ecc >= 1:
            return tc
    except ValueError:
        pass
    
    f = np.pi/2 - omega
    ee = 2 * np.arctan(np.tan(f/2) * np.sqrt((1-ecc)/(1+ecc)))  # eccentric anomaly
    tp = tc - per/(2*np.pi) * (ee - ecc*np.sin(ee))      # time of periastron
    
    return tp


# simultaneously fit multiple data sets with global and local parameters
class rv_fitter(lc_fitter):
    # formality for lc_fitter
    # needed to get posterior plot
    ns_type = 'ultranest' 

    def __init__(self, input_rv_data, global_bounds, local_rv_bounds, verbose=False):
        # keys for input_data: time, flux, ferr, airmass, priors all numpy arrays
        self.data = copy.deepcopy(input_rv_data)
        self.global_bounds = global_bounds
        self.local_bounds = local_rv_bounds
        self.verbose = verbose

        self.fit_nested()

    def fit_nested(self):
        """ Fit the data using nested sampling """

        # create bound arrays for generating samples
        nobs = len(self.data)
        gfreekeys = list(self.global_bounds.keys())

        lfreekeys = []
        boundarray = [self.global_bounds[k] for k in gfreekeys]
        alltime = []
        for i in range(nobs):
            lfreekeys.append(list(self.local_bounds[i].keys()))
            boundarray.extend([self.local_bounds[i][k] for k in lfreekeys[-1]])
            alltime.extend(self.data[i]['time'])
        boundarray = np.array(boundarray)
        alltime = np.array(alltime)
        dalltime = alltime - alltime.min()
        # make global time array and masks for each data set
        tmask = np.zeros(len(alltime), dtype=bool)
        idx = 0
        for i in range(nobs):
            self.data[i]['time_mask'] = tmask.copy()
            self.data[i]['time_mask'][idx:idx+len(self.data[i]['time'])] = True
            idx += len(self.data[i]['time'])

        # transform unit cube to prior volume
        bounddiff = np.diff(boundarray,1).reshape(-1)
        def prior_transform(upars):
            return (boundarray[:,0] + bounddiff*upars)

        def loglike(pars):
            chi2 = 0

            # global keys
            for j, key in enumerate(gfreekeys):
                self.data[0]['priors'][key] = pars[j]

            # compute mass ratio
            mtotal = Msun*self.data[0]['priors']['mstar'] + \
                     Mjup*self.data[0]['priors']['mplanet'] # kg
            mu = self.data[0]['priors']['mplanet']*Mjup/mtotal

            # semi-major axis using Kepler's 3rd law
            semimajor = (Grav*mtotal*self.data[0]['priors']['per']**2/4/np.pi**2)**(1/3) # m

            # estimate a/Rs
            ars = semimajor/(self.data[0]['priors']['rstar']*Rsun)

            # apply precession correction
            #omega_offset = dalltime * self.data[0]['priors']['rv_linear'] + \
            #        self.data[0]['priors']['rv_quad']*dalltime**2

            # priors should be the same except for first idx, we'll fix later
            self.data[0]['priors']['ars'] = ars
            self.data[0]['priors']['mu'] = mu
            #self.data[0]['priors']['omega'] = omega + omega_offset
            orbit = rv_model(alltime, self.data[0]['priors'])
            
            # apply trend to RV data
            orbit += dalltime * self.data[0]['priors']['rv_linear'] + \
                    self.data[0]['priors']['rv_quad']*dalltime**2

            # for each light curve
            for i in range(nobs):
                # set for consistency purposes, not really needed though
                self.data[i]['priors']['mu'] = mu
                self.data[i]['priors']['ars'] = ars

                # local keys
                ti = sum([len(self.local_bounds[k]) for k in range(i)])
                for j, key in enumerate(lfreekeys[i]):
                    self.data[i]['priors'][key] = pars[j+ti+len(gfreekeys)]

                # extract relevant part of the rv model
                model = orbit[self.data[i]['time_mask']]

                # handle offset
                detrend = self.data[i]['vel'] - model
                model += np.mean(detrend)
                chi2 += np.mean(((self.data[i]['vel']-model)/(self.data[i]['velerr']))**2)/nobs
                # penalize for jitter - changes likelihood surface too much to parameterize
                #chi2 += np.log(self.data[i]['priors']['jitter']*np.ones_like(self.data[i]['velerr'])).sum()

            # maximization metric for nested sampling
            return -0.5*chi2

        freekeys = []+gfreekeys
        for n in range(nobs):
            for k in lfreekeys[n]:
                freekeys.append(f"local_{n}_{k}")

        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=6e5)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=6e5, show_status=self.verbose, viz_callback=self.verbose)

        self.parameters = {}
        self.quantiles = {}
        self.errors = {}

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # compute some ratios
        mtotal = Msun*self.parameters.get('mstar', self.data[0]['priors']['mstar']) + \
                 Mjup*self.parameters.get('mplanet', self.data[0]['priors']['mplanet']) # kg
        mu = self.parameters.get('mplanet', self.data[0]['priors']['mplanet'])*Mjup/mtotal
        self.parameters['mu'] = mu
        self.errors['mu'] = 0.01*mu

        # semi-major axis using Kepler's 3rd law
        semimajor = (Grav*mtotal*self.parameters['per']**2/4/np.pi**2)**(1/3) # m

        # estimate a/Rs
        ars = semimajor/(self.parameters.get('rstar', self.data[0]['priors']['rstar'])*Rsun)

        for n in range(nobs):
            self.data[n]['errors'] = {}
            self.data[n]['priors']['mu'] = mu
            self.data[n]['priors']['ars'] = ars

            for k in lfreekeys[n]:
                pkey = f"local_{n}_{k}"
                # replace with final parameters
                self.data[n]['priors'][k] = self.parameters[pkey]
                self.data[n]['errors'][k] = self.errors[pkey]

                if k == 'rprs' and 'rprs' not in freekeys:
                    self.parameters[k] = self.data[n]['priors'][k]
                    self.errors[k] = self.data[n]['errors'][k]

            dtime = self.data[n]['time'] - self.data[n]['time'].min()
            self.data[n]['model'] = rv_model(self.data[n]['time'], self.data[n]['priors']) + \
                                    self.data[n]['priors']['rv_linear']*dtime + \
                                    self.data[n]['priors']['rv_quad']*dtime**2

            detrend = self.data[n]['vel'] - self.data[n]['model']
            self.data[n]['priors']['offset'] = np.mean(detrend) # TODO use monte carlo
            self.data[n]['detrend'] = self.data[n]['vel'] - self.data[n]['priors']['offset']
            self.data[n]['residuals'] = self.data[n]['detrend'] - self.data[n]['model']

            # scale factor to get avg error to equal std of residuals
            self.data[n]['error_scale'] = np.std(self.data[n]['residuals'])/np.mean(self.data[n]['velerr'])

        # global up-scaled model
        self.alltime = np.linspace(min(alltime), max(alltime), 100000)
        dalltime = self.alltime - self.alltime.min()
        rv = rv_model(self.alltime, self.data[0]['priors'])

        # apply trend to RV data
        self.allmodel = rv + self.data[n]['priors']['rv_linear']*dalltime + \
                             self.data[n]['priors']['rv_quad']*dalltime**2
        self.allphase = (self.alltime-self.data[n]['priors']['tmid'])/self.parameters['per']

        # compute acceleration of the star
        self.acceleration = acceleration_model(self.alltime, self.data[0]['priors'])

        self.K = (rv.max()-rv.min())/2

    def planet_star_distance(self, time):
        planet= planet_orbit(self.data[0]['priors']['per'],  self.data[0]['priors']['ars'],  self.data[0]['priors']['ecc'], 
                            self.data[0]['priors']['inc'],  self.data[0]['priors']['omega'],  self.data[0]['priors']['tmid'], 
                            time, mu=1-self.data[0]['priors']['mu'], ww=0)
    
        star = planet_orbit(self.data[0]['priors']['per'],  self.data[0]['priors']['ars'],  self.data[0]['priors']['ecc'], 
                            self.data[0]['priors']['inc'],  self.data[0]['priors']['omega'],  self.data[0]['priors']['tmid'], 
                            time, mu=self.data[0]['priors']['mu'], ww=0)

        # flip the star so it's on the other side of center of mass      
        star[0] *= -1
        star[1] *= -1
        star[2] *= -1

        return np.sqrt((planet[0]-star[0])**2 + (planet[1]-star[1])**2 + (planet[2]-star[2])**2)

    def plot_bestfit(self, title=""):
        """Plot the best-fit model and residuals for each dataset."""

        fig, ax = plt.subplots(3,figsize=(10,11))
        
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','magenta','blue','cyan','lime','gold','red',])

        for n in range(len(self.data)):
            phase = (self.data[n]['time'] - self.data[n]['priors']['tmid'])/self.parameters['per']
            phase += 0.5
            phase %= 1
            phase -= 0.5
            ncolor = next(colors)
            nmarker = next(markers)

            ax[1].errorbar(phase, self.data[n]['detrend'], yerr=self.data[n]['velerr'], 
                        marker=nmarker, color=ncolor, ls='')
            ax[2].errorbar(phase, self.data[n]['residuals'], yerr=self.data[n]['velerr'], 
                        marker=nmarker, color=ncolor, ls='', 
                        label=rf"{self.data[n]['name']} $\sigma$ = {np.std(self.data[n]['residuals']):.2f} m/s")
            ax[0].errorbar(self.data[n]['time']-int(self.alltime.min()), self.data[n]['detrend'], yerr=self.data[n]['velerr'], 
                        marker=nmarker, color=ncolor, ls='', 
                        label=rf"{self.data[n]['name']} $\sigma$ = {np.std(self.data[n]['residuals']):.2f} m/s")

        ax[1].axhline(0, color='k', ls='--', alpha=0.5)
        ax[1].set_xlim([-0.5,0.5])
        ax[2].set_xlim([-0.5,0.5])
        ax[0].set_xlim([self.alltime.min()-int(self.alltime.min()), self.alltime.max()-int(self.alltime.min())])
        nphase = (self.allphase+0.5)%1-0.5
        si = np.argsort(nphase)

        label = rf"$K$ = {self.K:.2f} m/s" + "\n" \
                rf"$P$ = {self.parameters['per']:.4f} $\pm$ {self.errors['per']:.2e}" + "\n" \
                rf"$ecc$ = {self.parameters.get('ecc', self.data[0]['priors']['ecc']):.4f} $\pm$ {self.errors.get('ecc',0):.4f}" + "\n" \
                rf"$\omega$ = {self.parameters.get('omega',self.data[0]['priors']['omega']):.2f} $\pm$ {self.errors.get('omega',0):.2f}" + "\n" \
                rf"$M_p$ = {self.parameters['mplanet']:.4f} $\pm$ {self.errors['mplanet']:.4f}"+r"$M_{Jup}$"

        ax[1].plot(nphase[si], self.allmodel[si], 'k-', label=label)
        ax[0].plot(self.alltime-int(self.alltime.min()), self.allmodel, 'k-', label='', alpha=0.75)
        ax[0].set_xlabel(f'BJD-{int(self.alltime.min())}')
        ax[1].set_ylabel("RV (m/s)")
        ax[0].set_ylabel("RV (m/s)")
        ax[2].set_ylabel("Residuals (m/s)")
        ax[0].set_xlabel(f"Days Since {int(self.alltime.min())}")
        ax[1].set_xlabel("Phase")
        ax[2].set_xlabel("Phase")
        ax[0].set_title(title)
        ax[2].grid(ls='--')
        ax[2].legend(loc='upper right')
        ax[1].legend(loc='best')
        return fig, ax

    
    def plot_bestfit_acceleration(self, title="", phase_limits=[-0.5,0.5]):
        """Plot the best-fit model with acceleration"""

        fig, ax = plt.subplots(3,figsize=(10,11))
        
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','magenta','blue','cyan','lime','gold','red',])

        for n in range(len(self.data)):
            phase = (self.data[n]['time'] - self.data[n]['priors']['tmid'])/self.parameters['per']
            phase += 0.5
            phase %= 1
            phase -= 0.5
            ncolor = next(colors)
            nmarker = next(markers)

            #ax[1].errorbar(phase, self.data[n]['detrend'], yerr=self.data[n]['velerr'], 
            #            marker=nmarker, color=ncolor, ls='')

            ax[0].errorbar(self.data[n]['time']-int(self.alltime.min()), self.data[n]['detrend'], yerr=self.data[n]['velerr'], 
                        marker=nmarker, color=ncolor, ls='', 
                        label=rf"{self.data[n]['name']} $\sigma$ = {np.std(self.data[n]['residuals']):.2f} m/s")

        ax[1].axhline(0, color='k', ls='--', alpha=0.5)
        ax[1].set_xlim([-0.5,0.5])
        ax[2].set_xlim([-0.5,0.5])
        ax[0].set_xlim([self.alltime.min()-int(self.alltime.min()), self.alltime.max()-int(self.alltime.min())])
        nphase = (self.allphase+0.5)%1-0.5
        si = np.argsort(nphase)

        label = rf"$K$ = {self.K:.2f} m/s" + "\n" \
                rf"$P$ = {self.parameters['per']:.4f} $\pm$ {self.errors['per']:.2e}" + "\n" \
                rf"$ecc$ = {self.parameters['ecc']:.4f} $\pm$ {self.errors['ecc']:.4f}" + "\n" \
                rf"$\omega$ = {self.parameters['omega']:.2f} $\pm$ {self.errors['omega']:.2f}" + "\n" \
                rf"T$_{{mid}}$ = {self.data[n]['priors']['tmid']:.4f}"
                #rf"$M_p$ = {self.parameters['mplanet']:.4f} $\pm$ {self.errors['mplanet']:.4f}"+r"$M_{Jup}$" \

        ax[0].plot(self.alltime-int(self.alltime.min()), self.allmodel, 'k-', label='', alpha=0.75)
        #ax[1].plot(nphase[si], self.allmodel[si]*(-1/self.parameters['mu'])/1000, 'k-', label=label)
        #ax[1].plot((nphase[si]-max_phase)*self.parameters['per'], self.allmodel[si]*(-1/self.parameters['mu'])/1000, 'k-', label=label)


        max_idx = self.allmodel[si].argmax()
        # compute smaller phase limits +/- 0.01 around the max
        phase_limits = np.array([(nphase[si][max_idx]-0.01), (nphase[si][max_idx]+0.01)])

        max_time = self.alltime[si][max_idx]
        max_phase = nphase[si][max_idx]

        ax[1].plot((nphase[si]-max_phase)*self.parameters['per'], self.allmodel[si], 'k-', label=label)

        # plot the acceleration
        ax[2].plot((nphase[si]-max_phase)*self.parameters['per'], (1./24)*self.acceleration[si], 'k-', label=label)

        # plot vertical lines at the max
        ax[2].axvline((nphase[si][max_idx]-max_phase)*self.parameters['per'], color='r', ls='--', alpha=0.5, label=f"Periapsis ")

        # plot vertical lines at the eclipse
        emin = (nphase[si][max_idx]-max_phase)*self.parameters['per']-3.104/24-0.07
        emax = (nphase[si][max_idx]-max_phase)*self.parameters['per']-3.104/24+0.07

        # fill between vertical lines
        ax[1].fill_betweenx(ax[1].get_ylim(), emin, emax, color='c', alpha=0.5, label=f"Eclipse")
        ax[2].fill_betweenx(ax[2].get_ylim(), emin, emax, color='c', alpha=0.5, label=f"Eclipse")

        # plot eclipse
        ax[2].axvline((nphase[si][max_idx]-max_phase)*self.parameters['per']-3.104/24, color='g', ls='--', alpha=0.75, label=f"Mid-Eclipse")
        ax[1].axvline((nphase[si][max_idx]-max_phase)*self.parameters['per']-3.104/24, color='g', ls='--', alpha=0.75, label=f"Mid-Eclipse")
        
        # plot periapsis
        ax[1].axvline((nphase[si][max_idx]-max_phase)*self.parameters['per'], color='r', ls='--', alpha=0.5, label=f"Periapsis ")

        # plot the residuals
        ax[0].set_xlabel(f'BJD-{int(self.alltime.min())}')
        ax[0].set_ylabel("Radial Velocity of Star (m/s)")
        ax[0].set_xlabel(f"Days Since {int(self.alltime.min())}")
        ax[0].legend(loc='upper right')
        ax[0].set_title(title)
        #ax[1].set_xlabel("Phase")
        ax[1].set_xlabel(f"Time from Periapsis  [day]")

        ax[1].set_ylabel("Velocity of Star (m/s)")
        ax[1].legend(loc='best')

        #ax[1].set_xlim(phase_limits)
        ax[1].set_xlim((phase_limits-max_phase)*self.parameters['per'])

        ax[2].set_xlim((phase_limits-max_phase)*self.parameters['per'])
        ax[2].set_ylabel("Acceleration of Star (m/s/hr)")
        ax[2].set_xlabel(f"Time from Periapsis  [day]")
        ax[2].grid(ls='--')
        ax[1].grid(ls='--')

        return fig, ax


    def plot_orbit(self):
        """Plot the best-fit orbit in 3D."""
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')
        ax2 = fig.add_subplot(331,projection='3d')
        newtime = np.linspace(self.data[0]['priors']['tmid']+1, self.data[0]['priors']['tmid']+self.data[0]['priors']['per']+1, 10000)
        xs,ys,zs = planet_orbit( self.data[0]['priors']['per'],  self.data[0]['priors']['ars'],  self.data[0]['priors']['ecc'], 
                                self.data[0]['priors']['inc'],  self.data[0]['priors']['omega'],  self.data[0]['priors']['tmid'], 
                                newtime, mu= self.data[0]['priors']['mu'], ww=0)
        xp,yp,zp = planet_orbit( self.data[0]['priors']['per'],  self.data[0]['priors']['ars'],  self.data[0]['priors']['ecc'], 
                                self.data[0]['priors']['inc'],  self.data[0]['priors']['omega'],  self.data[0]['priors']['tmid'], 
                                newtime, mu=1- self.data[0]['priors']['mu'], ww=0)
        xs *=-1
        ys *=-1
        zs *=-1
        distance = np.sqrt((xs-xp)**2+(ys-yp)**2+(zs-zp)**2)
        #im = ax.scatter(xs, ys, zs, c=newtime, s=10, zorder=1, cmap='viridis')
        # highlight 1 day sections of orbit? 
        edata = eclipse(newtime, self.data[0]['priors'])
        tdata = transit(newtime, self.data[0]['priors'])
        mide = np.argmin(edata)
        midt = np.argmin(tdata)
        # time of periastron - integral is wrong...
        omega = np.linspace(0, -np.pi/2- self.data[0]['priors']['omega']*np.pi/180,10000)
        dt = (self.data[0]['priors']['per']/(2*np.pi*np.sqrt(1-self.data[0]['priors']['ecc']**2)))
        fn = (1-self.data[0]['priors']['ecc']**2)/(1+self.data[0]['priors']['ecc']*np.cos(omega))**2
        integral = np.trapz(fn,omega)*dt
        tperi = newtime[mide]-integral
        midp = np.argmin(np.abs(newtime-tperi))
        tperi2 = newtime[np.argmin(distance)]
        midp2 = np.argmin(np.abs(newtime-tperi2)) # numerical solution
        tperi3 = timetrans_to_timeperi(newtime[midt], self.data[0]['priors']['per'], self.data[0]['priors']['ecc'], self.data[0]['priors']['omega']*np.pi/180) + self.data[0]['priors']['per']

        print("mid eclipse: ", newtime[mide])
        print("mid transit: ", newtime[midt])
        print("mid periastron: ", newtime[midp2])
        print("numerical dt = ", newtime[mide]-newtime[midp2])
        print("integral dt = ", integral)
        print(f"Integral: {tperi}")
        print(f"Numerical: {tperi2}")
        print(f"RADVEL: {tperi3}")
        midp3 = np.argmin(np.abs(newtime-tperi3))

        cmap = plt.cm.jet  # define the colormap
        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)][::-1]
        # force the first&last color entry to be grey
        cmaplist[0] = (.5, .5, .5, 1.0)
        cmaplist[-1] = (.5, .5, .5, 1.0)

        # create the new map
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'Custom cmap', cmaplist, cmap.N)

        # define the bins and normalize
        # bins of 4 hours
        tmid = newtime[midt]
        emid = newtime[mide]
        DT = tmid - emid
        nbins = int(np.round(DT)+4)
        bounds = np.round(np.linspace(emid-1.5, tmid+2, nbins),1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        im = ax.scatter(xp, yp, zp, c=newtime, cmap=cmap, norm=norm,  s=10, zorder=1)
        ax.scatter(xp[mide], yp[mide], zp[mide]+1, c='white', edgecolor='black', marker='o', s=75, zorder=4, label='Eclipse')
        ax.scatter(xp[midp2], yp[midp2], zp[midp2]+1, c='white', edgecolor='black', marker='s', s=75, zorder=4, label='Periastron')
        ax.scatter(xp[midt], yp[midt], zp[midt]+1, c='white', edgecolor='black', marker='P', s=75, zorder=4, label='Transit')
        ax.scatter(1,0,0, c='white', edgecolor='black', marker='*', s=100, zorder=4, label='Star')
        ax.scatter(100,0,0, c='white', edgecolor='black', marker='v', s=150, zorder=4, label='View from Earth')
        # hide z-axi1 labels and ticks
        ax.set_zticklabels([])
        ax.set_zlim([-100,100])
        ax.set_xlim([-100,100])
        ax.set_ylim([-100,100])
        ax.set_xlabel('X [R_Star]')
        ax.set_ylabel('Y [R_Star]')
        #ax.set_zlabel('Z [R_Star]')
        ax.set_title("Orbit of HD 80606 b")
        ax.legend(loc='best')
        ax2.scatter(xs, ys, zs, c=newtime, cmap=cmap, norm=norm, s=10, zorder=3)
        ax2.set_zlim([-1,1])
        ax2.set_xlim([-1,1])
        ax2.set_ylim([-1,1])
        ax2.set_zticklabels([])

        ax.set_xlim([-20,80])
        ax.set_ylim([-50,50])
        ax.view_init(90, 0)
        ax2.view_init(90, 0)
        cbar = plt.colorbar(im)
        cbar.set_label('Time [BJD]')

        return fig
    
    