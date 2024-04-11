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
from copy import deepcopy
from itertools import cycle
import matplotlib.pyplot as plt
import numpy as np
from pylightcurve.models.exoplanet_lc import eclipse_mid_time, transit_flux_drop
from scipy import stats
try:
    from ultranest import ReactiveNestedSampler
except ImportError:
    import dynesty
    import dynesty.plotting
    from dynesty.utils import resample_equal
    from scipy.stats import gaussian_kde

try:
    from elca import glc_fitter, lc_fitter
except ImportError:
    from .elca import glc_fitter, lc_fitter

AU = const.au.to(u.m).value
Mjup = const.M_jup.to(u.kg).value
Msun = const.M_sun.to(u.kg).value
Rsun = const.R_sun.to(u.m).value
Grav = const.G.to(u.m**3/u.kg/u.day**2).value


def planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0, mu=1, W=0):
    # see original @ https://github.com/ucl-exoplanets/pylightcurve/blob/master/pylightcurve/models/exoplanet_lc.py
    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0
    W = W * np.pi / 180.0

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
    rr = mu*sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv+W))

    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)

    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))

    return [x, y, z]

def pytransit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
            mid_time, time_array, method='claret', precision=3):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                             method=method, precision=precision)

def transit(times, values):
    model = pytransit([values['u0'], values['u1'], values['u2'], values['u3']], 
                    values['rprs'], values['per'], values['ars'], 
                    values['ecc'], values['inc'], values['omega'],
                    values['tmid'], times, method='claret', precision=3)
    return model

from pylightcurve.models.exoplanet_lc import transit as pytransit
from pylightcurve.models.exoplanet_lc import eclipse_mid_time

def eclipse(times, values):
    tme = eclipse_mid_time(values['per'], values['ars'], values['ecc'], values['inc'], values['omega'], values['tmid'])
    model = pytransit([0,0,0,0], 
                    values['rprs']*values['fpfs']**0.5, values['per'], values['ars'], 
                    values['ecc'], values['inc'], values['omega']+180,
                    tme, times, method='quad', precision=3)
    return model + values['fpfs']*values['rprs']**2 # second part is eclipse depth

def get_phase(times, per, tmid):
    return (times - tmid + 0.25 * per) / per % 1 - 0.25

def phasecurve(times, values):
    tme = eclipse_mid_time(values['per'], values['ars'], values['ecc'], values['inc'], values['omega'], values['tmid'])
    edepth = values['fpfs']*values['rprs']**2 # eclipse depth
    emodel = pytransit([0,0,0,0], 
                    edepth**0.5, values['per'], values['ars'], 
                    values['ecc'], values['inc'], values['omega']+180,
                    tme, times, method='quad', precision=3)
    tmodel = pytransit([values['u1'], values['u2'], 0, 0], 
                    values['rprs'], values['per'], values['ars'], 
                    values['ecc'], values['inc'], values['omega'],
                    values['tmid'], times, method='quad', precision=3)
    c0 = edepth - values['c1'] - values['c3']
    brightness = 1 + c0 + values['c1']*np.cos(2*np.pi*(times-tme)/values['per']) + values['c2']*np.sin(2*np.pi*(times-tme)/values['per']) + values['c3']*np.cos(4*np.pi*(times-tme)/values['per']) + values['c4']*np.sin(4*np.pi*(times-tme)/values['per'])
    emask = np.floor(emodel)
    return (brightness*(emask) + (edepth+emodel)*(1-emask))*tmodel

def rv_model(time, params, dt=0.0001):
    xp,yp,zp = planet_orbit(params['per'], params['ars'], params['ecc'], 
                            params['inc'], params['omega'], params['tmid'], 
                            time, mu=1, ww=0)
    # TODO optimize by only computing in X direction
    xp2,yp2,zp2 = planet_orbit(params['per'], params['ars'], params['ecc'],
                            params['inc'], params['omega'], params['tmid'],
                            time+dt, mu=1, ww=0)
    # (const.R_sun / 0.0001 / u.day).to(u.m/u.s).value
    # scale in order to get sensible units
    return (xp2-xp)*params['mu']*params['rstar']*80520833.33333333 # Rstar/day -> m/s

    # something like a velocity vector
    #v = np.array([np.diff(xp), np.diff(yp), np.diff(zp)]).T 
    # transits are looking in -X direction
    # need to flip to get direction of stellar motion correct after scaling
    #view = np.array([1,0,0]) # looking in +X direction
    #dot = np.dot(v*params['mu'], view) # project with line of sight
    #rv = dot*(params['rstar']*const.R_sun / dt / u.day).to(u.m/u.s)
    #K = (rv.max() - rv.min())/2


# simultaneously fit multiple data sets with global and local parameters
class joint_fitter(glc_fitter):
    ns_type = 'ultranest'

    def __init__(self, lc_data, local_lc_bounds, rv_data, local_rv_bounds, global_bounds, ephemeris, individual_fit=True, verbose=False):
        # keys for input_data: time, flux, ferr, airmass, priors all numpy arrays
        self.lc_data = lc_data
        self.local_lc_bounds = local_lc_bounds
        self.rv_data = rv_data
        self.local_rv_bounds = local_rv_bounds
        self.individual_fit = individual_fit
        self.global_bounds = global_bounds
        self.ephemeris = ephemeris
        self.verbose = verbose

        self.fit_nested()

    def fit_nested(self):

        # create bound arrays for generating samples
        nobs_lc = len(self.lc_data)
        nobs_rv = len(self.rv_data)

        # fit individual light curves to constrain priors
        if self.individual_fit:
            for i in range(nobs_lc):

                # spit out individual t-mids

                print(f"Fitting individual light curve {i+1}/{nobs_lc}")
                mybounds = dict(**self.local_lc_bounds[i])#, **self.global_bounds)
                if 'inc' in self.global_bounds:
                    mybounds['inc'] = self.global_bounds['inc']
                if 'tmid' in self.global_bounds:
                    mybounds['tmid'] = self.global_bounds['tmid']

                myfit = lc_fitter(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'],
                    self.lc_data[i]['ferr'],
                    self.lc_data[i]['airmass'],
                    self.lc_data[i]['priors'],
                    mybounds
                )

                # update global bounds with +/- 5 sigma
                for j, key in enumerate(self.local_lc_bounds[i].keys()):
                    if 'rprs' in key:
                        self.local_lc_bounds[i][key] = [
                            max(0, myfit.parameters[key] - 10*myfit.errors[key]),
                            myfit.parameters[key] + 10*myfit.errors[key]
                        ]
                    else:
                        self.local_lc_bounds[i][key] = [
                            myfit.parameters[key] - 5*myfit.errors[key],
                            myfit.parameters[key] + 5*myfit.errors[key]
                        ]
                # if using gaussian prior on inclination
                if 'inc' in self.global_bounds:
                    self.global_bounds['inc'] = [
                        myfit.parameters['inc']- myfit.errors['inc']*10,
                        min(myfit.parameters['inc']+ myfit.errors['inc']*10, 90)
                    ]                
                del(myfit)

        # keys for global bounds
        gfreekeys = list(self.global_bounds.keys())

        # keys for local bounds
        lfreekeys = []
        boundarray = [self.global_bounds[k] for k in gfreekeys]
        alltime = [] # array for calling rv orbit eq once

        # collect RV data
        for i in range(nobs_rv):
            lfreekeys.append(list(self.local_rv_bounds[i].keys()))
            boundarray.extend([self.local_rv_bounds[i][k] for k in lfreekeys[-1]])
            alltime.extend(self.rv_data[i]['time'])
        alltime = np.array(alltime)
        dalltime = alltime - alltime.min()

        # collect LC data
        for i in range(nobs_lc):
            lfreekeys.append(list(self.local_lc_bounds[i].keys()))
            boundarray.extend([self.local_lc_bounds[i][k] for k in lfreekeys[-1]])
        boundarray = np.array(boundarray)

        print(boundarray)

        # make global time array and time masks for each rv data set
        tmask = np.zeros(len(alltime), dtype=bool)
        idx = 0
        for i in range(nobs_rv):
            self.rv_data[i]['time_mask'] = tmask.copy()
            self.rv_data[i]['time_mask'][idx:idx+len(self.rv_data[i]['time'])] = True
            idx += len(self.rv_data[i]['time'])

        # transform unit cube to prior volume
        bounddiff = np.diff(boundarray,1).reshape(-1)

        # create a gaussian priors, if any
        gaussian_priors = {}
        for k in ['rstar','mstar']:
            if k in gfreekeys:
                gaussian_priors[k] = stats.norm(
                    self.global_bounds[k][0],
                    self.global_bounds[k][1])
                gaussian_priors[k+"_idx"] = gfreekeys.index(k)
            #def transform_1d(quantile):
            #    return gaussdistribution.ppf(quantile)
            def prior_transform(upars):
                newpars = (boundarray[:,0] + bounddiff*upars) # uniform
                for k in gaussian_priors.keys():
                    if 'idx' in k:
                        continue
                    newpars[gaussian_priors[k+"_idx"]] = gaussian_priors[k].ppf(upars[gaussian_priors[k+"_idx"]])
                    #newpars[staridx] = transform_1d(upars[staridx])
                return newpars

        # estimate orbit numbers for ephemeris
        if self.ephemeris is not None:
            self.ephemeris['tmid_orbit'] = np.round((self.ephemeris['tmid'][:,0]-self.ephemeris['prior']['tmid'][0])/self.ephemeris['prior']['per'][0])
            self.ephemeris['historic_orbit'] = np.round((self.ephemeris['historic'][:,0]-self.ephemeris['prior']['tmid'][0])/self.ephemeris['prior']['per'][0])
            self.ephemeris['emid_orbit'] = np.round((self.ephemeris['emid'][:,0]-self.ephemeris['prior']['emid'][0])/self.ephemeris['prior']['per'][0])

        def loglike(pars):
            rv_chi2 = 0  # radial velocity
            lc_chi2 = 0  # transit light curves
            eph_chi2 = 0 # ephemeris

            # set parameters
            for j, key in enumerate(gfreekeys):
                self.rv_data[0]['priors'][key] = pars[j]

            # compute mass ratio
            mtotal = Msun*self.rv_data[0]['priors']['mstar'] + \
                     Mjup*self.rv_data[0]['priors']['mplanet'] # kg
            mu = self.rv_data[0]['priors']['mplanet']*Mjup/(Msun*self.rv_data[0]['priors']['mstar'])

            # semi-major axis using Kepler's 3rd law
            semimajor = (Grav*mtotal*self.rv_data[0]['priors']['per']**2/4/np.pi**2)**(1/3) # m

            # estimate a/Rs
            ars = semimajor/(self.rv_data[0]['priors']['rstar']*Rsun)

            # compute ecc and omega - old
            #omega_rad = np.arctan(self.rv_data[0]['priors']['esinw']/self.rv_data[0]['priors']['ecosw'])
            #omega = omega_rad*180/np.pi
            #ecc = self.rv_data[0]['priors']['ecosw']/np.cos(omega_rad)

            # priors should be the same, except for that first idx, we'll fix later
            self.rv_data[0]['priors']['ars'] = ars
            self.rv_data[0]['priors']['mu'] = mu

            # global rv model
            orbit = rv_model(alltime, self.rv_data[0]['priors'])

            # apply trend to RV data
            orbit += dalltime * self.rv_data[0]['priors']['rv_linear'] + \
                    self.rv_data[0]['priors']['rv_quad']*dalltime**2

            # for each RV dataset compute chi2
            for i in range(nobs_rv):

                # set global parameters
                self.rv_data[i]['priors']['mu'] = mu
                self.rv_data[i]['priors']['ars'] = ars

                for j, key in enumerate(gfreekeys):
                    self.rv_data[i]['priors'][key] = pars[j]

                # set local parameters - TODO remove or test
                ti = sum([len(self.local_rv_bounds[k]) for k in range(i)])
                for j, key in enumerate(lfreekeys[i]):
                    self.rv_data[i]['priors'][key] = pars[j+ti+len(gfreekeys)]

                # extract relevant part of the rv model
                model = orbit[self.rv_data[i]['time_mask']]

                # handle offset
                detrend = self.rv_data[i]['vel'] - model
                model += np.mean(detrend)

                # TODO add error scaling to chi2
                rv_chi2 += np.sum(((self.rv_data[i]['vel']-model)/(self.rv_data[i]['velerr']))**2)#/nobs_rv

            # for each LC dataset compute chi2
            for i in range(nobs_lc):

                # set global parameters
                self.lc_data[i]['priors']['mu'] = mu
                self.lc_data[i]['priors']['ars'] = ars

                for j, key in enumerate(gfreekeys):
                    self.lc_data[i]['priors'][key] = pars[j]

                # set local parameters
                ti = sum([len(self.local_lc_bounds[k]) for k in range(i)])
                for j, key in enumerate(lfreekeys[i+nobs_rv]):
                    self.lc_data[i]['priors'][key] = pars[j+ti+len(gfreekeys)]

                # compute lc model
                model = transit(self.lc_data[i]['time'], self.lc_data[i]['priors'])

                # handle offset
                detrend = self.lc_data[i]['flux']/model
                model *= np.mean(detrend)
                lc_chi2 += 2*np.sum(((self.lc_data[i]['flux']-model)/(self.lc_data[i]['ferr']))**2)

            # add ephemeris into chi2
            self.ephemeris['tmid_orbit'] = np.round((self.ephemeris['tmid'][:,0]-self.lc_data[0]['priors']['tmid'])/self.rv_data[0]['priors']['per'])

            # predict mid-transit time
            tmid_pred = self.ephemeris['tmid_orbit']*self.rv_data[0]['priors']['per'] + self.lc_data[0]['priors']['tmid']
            eph_chi2 += np.sum(((tmid_pred-self.ephemeris['tmid'][:,0])/(self.ephemeris['tmid'][:,1]+self.ephemeris['noise']))**2)

            # predict mid-eclipse time [slow but accurate]
            emid = eclipse_mid_time(
                self.rv_data[0]['priors']['per'], 
                self.rv_data[0]['priors']['ars'], 
                self.rv_data[0]['priors']['ecc'], 
                self.lc_data[0]['priors']['inc'], 
                self.rv_data[0]['priors']['omega'], 
                self.lc_data[0]['priors']['tmid'])

            self.ephemeris['emid_orbit'] = np.round((self.ephemeris['emid'][:,0]-emid)/self.rv_data[0]['priors']['per'])

            emid_pred = self.ephemeris['emid_orbit']*self.rv_data[0]['priors']['per'] + emid

            #eph_chi2 += np.sum(((emid_pred-self.ephemeris['emid'][:,0])/(self.ephemeris['emid'][:,1]+self.ephemeris['noise']))**2)
            eph_chi2 += np.sum(((emid_pred-self.ephemeris['emid'][:,0])/(self.ephemeris['emid'][:,1]))**2)

            # maximization metric for nested sampling
            return -0.5*(rv_chi2 + lc_chi2 + eph_chi2)

        # make labels for posterior plot
        freekeys = []+gfreekeys
        for n in range(nobs_rv+nobs_lc):
            for k in lfreekeys[n]:
                freekeys.append(f"local_{n}_{k}")

        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=2e5)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=2e5, show_status=self.verbose, viz_callback=self.verbose)

        try:
            self.parameters = deepcopy(self.lc_data[0]['priors'])
        except:
            self.parameters = deepcopy(self.rv_data[0]['priors'])

        self.parameters_median = {}
        self.quantiles = {}
        self.errors = {}

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            if key == 'rstar': self.parameters[key] = self.results['posterior']['median'][i]
            self.parameters_median[key] = self.results['posterior']['median'][i]

            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # compute some ratios
        mtotal = Msun*self.parameters.get('mstar', self.rv_data[0]['priors']['mstar']) + \
                 Mjup*self.parameters.get('mplanet', self.rv_data[0]['priors']['mplanet']) # kg
        mu = self.parameters.get('mplanet', self.rv_data[0]['priors']['mplanet'])*Mjup / \
                (self.parameters.get('mstar', self.rv_data[0]['priors']['mstar'])*Msun)

        # semi-major axis using Kepler's 3rd law
        semimajor = (Grav*mtotal*self.parameters['per']**2/4/np.pi**2)**(1/3) # m

        # estimate a/Rs
        ars = semimajor/(self.parameters.get('rstar', self.rv_data[0]['priors']['rstar'])*Rsun)

        for n in range(nobs_rv): # TODO
            self.rv_data[n]['errors'] = {}
            self.rv_data[n]['priors']['mu'] = mu
            self.rv_data[n]['priors']['ars'] = ars

            # set global parameters
            for k in gfreekeys:
                self.rv_data[n]['priors'][k] = self.parameters[k]
                self.rv_data[n]['errors'][k] = self.errors[k]

            # set local parameters
            for k in lfreekeys[n]:
                pkey = f"local_{n}_{k}"
                # replace with final parameters
                self.rv_data[n]['priors'][k] = self.parameters[pkey]
                self.rv_data[n]['errors'][k] = self.errors[pkey]

            dtime = self.rv_data[n]['time'] - alltime.min()
            self.rv_data[n]['model'] = rv_model(self.rv_data[n]['time'], self.rv_data[n]['priors']) + \
                dtime * self.parameters.get('rv_linear',0) + \
                dtime**2 * self.parameters.get('rv_quad',0)

            detrend = self.rv_data[n]['vel'] - self.rv_data[n]['model']
            self.rv_data[n]['priors']['offset'] = np.mean(detrend) # TODO use monte carlo
            self.rv_data[n]['detrend'] = self.rv_data[n]['vel'] - self.rv_data[n]['priors']['offset']
            self.rv_data[n]['residuals'] = self.rv_data[n]['detrend'] - self.rv_data[n]['model']

        # global up-scaled rv model
        self.rv_time = np.linspace(min(alltime), max(alltime), 100000)
        self.dalltime = self.rv_time - self.rv_time.min()
        rv = rv_model(self.rv_time, self.rv_data[0]['priors'])
        self.rv_model = rv + self.parameters.get('rv_linear',0)*self.dalltime + \
                             self.parameters.get('rv_quad',0)*self.dalltime**2
        self.rv_phase = (self.rv_time-self.parameters['tmid'])/self.parameters['per']

        time = []
        data = []
        detrended = []
        models = []
        for n in range(nobs_lc):
            self.lc_data[n]['errors'] = {}
            self.lc_data[n]['priors']['mu'] = mu
            self.lc_data[n]['priors']['ars'] = ars

            for k in gfreekeys:
                self.lc_data[n]['priors'][k] = self.parameters[k]
                self.lc_data[n]['errors'][k] = self.errors[k]

            for k in lfreekeys[n+nobs_rv]:
                pkey = f"local_{n+nobs_rv}_{k}"
                self.lc_data[n]['priors'][k] = self.parameters[pkey]
                self.lc_data[n]['errors'][k] = self.errors[pkey]
                if k == 'rprs' and 'rprs' not in freekeys:
                    self.parameters[k] = self.lc_data[n]['priors'][k]
                    self.errors[k] = self.lc_data[n]['errors'][k]

            model = transit(self.lc_data[n]['time'], self.lc_data[n]['priors'])
            airmass = np.exp(self.lc_data[n]['airmass']*self.lc_data[n]['priors']['a2'])
            detrend = self.lc_data[n]['flux']/(model*airmass)
            self.lc_data[n]['priors']['a1'] = np.mean(detrend)
            self.lc_data[n]['residuals'] = self.lc_data[n]['flux'] - model*airmass*self.lc_data[n]['priors']['a1']
            self.lc_data[n]['detrend'] = self.lc_data[n]['flux']/(airmass*self.lc_data[n]['priors']['a1'])
            self.lc_data[n]['model'] = model
            self.lc_data[n]['phase'] = get_phase(self.lc_data[n]['time'], self.parameters['tmid'], self.parameters['per'])
            time.extend(self.lc_data[n]['time'])
            data.extend(self.lc_data[n]['flux'])
            detrended.extend(detrend)
            models.extend(model)

        self.time = np.array(time)
        self.data = np.array(data)
        self.models = np.array(models)
        self.phase = (self.time-self.parameters['tmid'])/self.parameters['per']
        self.detrended = np.array(detrended)
        self.residuals = self.detrended - self.models

        self.mc_extras()

    def mc_extras(self):
        prior = self.rv_data[0]['priors'].copy()
        orbits = [] # rv model
        transits = []
        semi = []
        arss = []
        mratio = []
        eccs = []
        omegas = []
        emids = []

        offsets = {}
        for n in range(len(self.rv_data)):
            offsets[self.rv_data[n]['name']] = []

        # monte carlo over posteriors to derive offsets and a/Rs
        for n in range(1000):
            
            # TODO fix error on omega to reduce OC uncertainty for e-mid

            # randomize each parameter
            for i, key in enumerate(self.errors.keys()):
                #if key == 'omega':
                #    continue
                try:
                    prior[key] = np.random.normal(self.parameters_median[key], self.errors[key])
                except:
                    pass

            # compute some ratios
            mtotal = Msun*prior['mstar'] + Mjup*prior['mplanet'] # kg
            mu = prior['mplanet']*Mjup/(Msun*prior['mstar'])

            # semi-major axis using Kepler's 3rd law
            semimajor = (Grav*mtotal*prior['per']**2/4/np.pi**2)**(1/3) # m

            # estimate a/Rs
            ars = semimajor/(prior['rstar']*Rsun)
            prior['ars'] = ars
            prior['mu'] = mu

            emid = eclipse_mid_time(
                    prior['per'], 
                    prior['ars'], 
                    prior['ecc'], 
                    prior['inc'], 
                    prior['omega'], 
                    prior['tmid']-prior['per'])
            
            emids.append(emid)
            orbits.append(rv_model(np.linspace(self.rv_time.min(), self.rv_time.min()+prior['per'], 10000), prior))
            semi.append(semimajor)
            arss.append(ars)
            mratio.append(mu)

            # derive rv offsets
            for n in range(len(self.rv_data)):

                dtime = self.rv_data[n]['time'] - self.rv_data[n]['time'].min()
                model = rv_model(self.rv_data[n]['time'], prior) + \
                    dtime * prior.get('rv_linear',0) + \
                    dtime**2 * prior.get('rv_quad',0)

                detrend = self.rv_data[n]['vel'] - model
                offsets[self.rv_data[n]['name']].append(np.mean(detrend))
            
            # estimate transit duration
            for n in range(len(self.lc_data)):
                model = transit(self.lc_data[n]['time'], prior)
                tmask = model < 1
                dt = np.diff(self.lc_data[n]['time']).mean()
                duration = tmask.sum()*dt
                transits.append(duration)
                # airmass = np.exp(self.lc_data[n]['airmass']*self.lc_data[n]['priors']['a2'])
                # detrend = self.lc_data[n]['flux']/(model*airmass)
                # self.lc_data[n]['priors']['a1'] = np.mean(detrend)
                # self.lc_data[n]['residuals'] = self.lc_data[n]['flux'] - model*airmass*self.lc_data[n]['priors']['a1']

        for n in range(len(self.rv_data)):
            key = self.rv_data[n]['name']+"_offset"
            self.parameters[key] = np.mean(offsets[self.rv_data[n]['name']])
            self.errors[key] = np.std(offsets[self.rv_data[n]['name']])

        # a/Rs
        self.parameters['ars'] = np.median(arss)
        self.errors['ars'] = np.std(arss)

        # semimajor axis in AU
        self.parameters['a'] = np.median(semi)/AU
        self.errors['a'] = np.std(semi)/AU

        # mass ratio Mp/Ms
        self.parameters['mu'] = np.median(mratio)
        self.errors['mu'] = np.std(mratio)

        # mid -eclipse 
        self.parameters['emid'] = np.median(emids)
        self.errors['emid'] = np.std(emids)

        # rv semi-amplitude m/s
        orbit_K = (np.max(orbits,1) - np.min(orbits,1))/2
        self.parameters['K'] = np.median(orbit_K)
        self.errors['K'] = np.std(orbit_K)
        
        # duration in days
        self.parameters['T14'] = np.median(transits)
        self.errors['T14'] = np.std(transits)

        self.parameters['rprs2'] = self.parameters['rprs']**2
        self.errors['rprs2'] = 2*self.parameters['rprs']*self.errors['rprs']

    def plot_rv_bestfit(self):
        fig, ax = plt.subplots(3,figsize=(10,11))
        
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','magenta','blue','cyan','lime','gold','red',])

        for n in range(len(self.rv_data)):
            phase = (self.rv_data[n]['time'] - self.rv_data[n]['priors']['tmid'])/self.parameters['per']
            phase += 0.5
            phase %= 1
            phase -= 0.5
            ncolor = next(colors)
            nmarker = next(markers)

            ax[1].errorbar(phase, self.rv_data[n]['detrend'], yerr=self.rv_data[n]['velerr'], 
                        marker=nmarker, color=ncolor,alpha=0.75, ls='')
            ax[2].errorbar(phase, self.rv_data[n]['residuals'], yerr=self.rv_data[n]['velerr'], 
                        marker=nmarker, color=ncolor,alpha=0.75, ls='', label=rf"{self.rv_data[n]['name']} $\sigma$ = {np.std(self.rv_data[n]['residuals']):.2f} m/s")
            ax[0].errorbar(self.rv_data[n]['time']-int(self.rv_time.min()), self.rv_data[n]['detrend'], yerr=self.rv_data[n]['velerr'], 
                        marker=nmarker, color=ncolor,alpha=0.75, ls='', label=rf"{self.rv_data[n]['name']} $\sigma$ = {np.std(self.rv_data[n]['residuals']):.2f} m/s")
        ax[1].axhline(0, color='k', ls='--', alpha=0.5)
        ax[1].set_xlim([-0.5,0.5])
        ax[2].set_xlim([-0.5,0.5])
        ax[0].set_xlim([self.rv_time.min()-int(self.rv_time.min()), self.rv_time.max()-int(self.rv_time.min())])
        nphase = (self.rv_phase+0.5)%1-0.5
        si = np.argsort(nphase)
        label = rf"$K$ = {self.parameters['K']:.2f} $\pm$ {self.errors['K']:.2f} m/s" + "\n" \
                rf"$P$ = {self.parameters['per']:.5f} $\pm$ {self.errors['per']:.1e} day" + "\n" \
                rf"$e$ = {self.parameters.get('ecc',self.rv_data[0]['priors']['ecc']):.5f} $\pm$ {self.errors.get('ecc',0):.1e}" + "\n" \
                rf"$\omega$ = {self.parameters.get('omega',self.rv_data[0]['priors']['ecc']):.3f} $\pm$ {self.errors.get('omega',0):.3f} deg"   
        ax[1].plot(nphase[si], self.rv_model[si], 'k-', label=label, alpha=0.75, zorder =2)
        ax[0].plot(self.rv_time-int(self.rv_time.min()), self.rv_model, 'k-', label='', alpha=0.75)
        ax[0].set_xlabel(f'BJD-{int(self.rv_time.min())}')
        ax[1].set_ylabel("RV (m/s)")
        ax[0].set_ylabel("RV (m/s)")
        ax[2].set_ylabel("Residuals (m/s)")
        ax[0].set_xlabel(f"Days Since {int(self.rv_time.min())} BJD")
        ax[1].set_xlabel("Phase")
        ax[2].set_xlabel("Phase")
        ax[2].grid(ls='--')
        ax[0].legend(loc='upper right')
        ax[1].legend()
        ax[2].legend()
        plt.tight_layout()

        return fig, ax

    def plot_oc_transits(self):

        ############### O-C plot
        fig,ax = plt.subplots(1, figsize=(10,7))

        self.ephemeris['tmid_orbit'] = np.round((self.ephemeris['tmid'][:,0]-self.lc_data[0]['priors']['tmid'])/self.rv_data[0]['priors']['per'])
        self.ephemeris['historic_orbit'] = np.round((self.ephemeris['historic'][:,0]-self.lc_data[0]['priors']['tmid'])/self.rv_data[0]['priors']['per'])

        # predict mid-transit time
        tmid_pred = self.ephemeris['tmid_orbit']*self.rv_data[0]['priors']['per'] + self.lc_data[0]['priors']['tmid']
        tmid_residual = self.ephemeris['tmid'][:,0] - tmid_pred
        historic_pred = self.ephemeris['historic_orbit']*self.rv_data[0]['priors']['per'] + self.lc_data[0]['priors']['tmid']
        historic_residual = self.ephemeris['historic'][:,0] - historic_pred

        print("Transit Residual [min]:", tmid_residual*24*60)
        ax.errorbar(self.ephemeris['historic_orbit'], historic_residual*24*60, yerr=self.ephemeris['historic'][:,1]*24*60, ls='none', marker='^',label='Not Used',color='gray')
        
        ax.errorbar(self.ephemeris['tmid_orbit'], tmid_residual*24*60, yerr=self.ephemeris['tmid'][:,1]*24*60, ls='none', marker='s',label='Mid-Transit Measurement',color='black')
        ylower = (tmid_residual.mean()-3*np.std(tmid_residual))*24*60
        yupper = (tmid_residual.mean()+3*np.std(tmid_residual))*24*60

        # TODO plot each point individually
        #for i in range(len(self.lc_data)):
        ax.errorbar(0, 0, yerr=self.errors['tmid']*24*60, ls='none', marker='s', color='black')

        # upsample data
        epochs = (np.linspace(self.ephemeris['tmid_orbit'].min()-7, self.ephemeris['tmid_orbit'].max()+7, 1000)) 

        depoch = epochs.max() - epochs.min()
        ax.set_xlim([epochs.min()-depoch*0.01, epochs.max()+depoch*0.01])

        # best fit solution
        model = epochs*self.rv_data[0]['priors']['per'] + self.lc_data[0]['priors']['tmid']

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['per'], self.errors['per'], size=10000)
        mc_b = np.random.normal(self.parameters['tmid'], self.errors['tmid'], size=10000)
        mc_model = np.expand_dims(epochs,-1) * mc_m + mc_b

        # create a fill between area for uncertainty of new ephemeris
        diff = mc_model.T - model
        ax.fill_between(epochs, np.percentile(diff,16,axis=0)*24*60, np.percentile(diff,84,axis=0)*24*60, alpha=0.2, color='k', label=r'Uncertainty ($\pm$ 1$\sigma$)')

        # duplicate axis and plot days since mid-transit
        ax2 = ax.twiny()
        ax2.set_xlabel(f"Time [BJD - {self.parameters['tmid']:.1f}]",fontsize=14)
        ax2.set_xlim(ax.get_xlim())
        xticks = ax.get_xticks()
        dt = np.round(xticks*self.parameters['per'],1)
        ax2.set_xticklabels(dt)
        show_2sigma = False

        if self.ephemeris['prior'] is not None:
            # create fill between area for uncertainty of old/prior ephemeris
            epochs_p = ((epochs*self.rv_data[0]['priors']['per'] + self.lc_data[0]['priors']['tmid']) - self.ephemeris['prior']['tmid'][0])/self.ephemeris['prior']['per'][0]
            #epochs_p = (np.linspace(self.ephemeris['tmid'][:,0].min()-7, self.ephemeris['tmid'][:,0].max()+7, 1000) 
            prior_p = epochs_p*self.ephemeris['prior']['per'][0] + self.ephemeris['prior']['tmid'][0]
            mc_m_p = np.random.normal(self.ephemeris['prior']['per'][0], self.ephemeris['prior']['per'][1], size=10000)
            mc_b_p = np.random.normal(self.ephemeris['prior']['tmid'][0], self.ephemeris['prior']['tmid'][1], size=10000)
            mc_model_p = np.expand_dims(epochs_p,-1) * mc_m_p + mc_b_p
            diff_p = mc_model_p.T - model

            # plot an invisible line so the 2nd axes are happy
            #ax2.plot(epochs, (model-prior_p)*24*60, ls='--', color='r', alpha=1)

            if show_2sigma:
                ax.fill_between(epochs, np.percentile(diff_p,2,axis=0)*24*60, np.percentile(diff_p,98,axis=0)*24*60, alpha=0.1, color='r', label=r'Prior ($\pm$ 2$\sigma$)')
            else:
                # show ~1 sigma
                ax.fill_between(epochs, np.percentile(diff_p,36,axis=0)*24*60, np.percentile(diff_p,64,axis=0)*24*60, alpha=0.1, color='r', label=r'Prior ($\pm$ 1$\sigma$)')

            #if ylim == 'prior':
            #    ax.set_ylim([ min(np.percentile(diff_p,1,axis=0)*24*60),
            #                max(np.percentile(diff_p,99,axis=0)*24*60)])
            #elif ylim == 'average':
            #    ax.set_ylim([ 0.5*(min(np.percentile(diff,1,axis=0)*24*60) + min(np.percentile(diff_p,1,axis=0)*24*60)),
            #                0.5*(max(np.percentile(diff,99,axis=0)*24*60) + max(np.percentile(diff_p,99,axis=0)*24*60))])

        ax.axhline(0,color='black',alpha=0.5,ls='--',
                    label="Period: {:.5f}+-{:.5f} days\nT_mid: {:.4f}+-{:.4f} BJD".format(self.parameters['per'], self.errors['per'], np.round(self.parameters['tmid'],4), np.round(self.errors['tmid'],4)))

        # TODO sig figs
        #lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
        #    str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
        #    str(round_to_2(self.errors.get('tmid',0)))
        #)

        ax.legend(loc='best')
        ax.set_xlabel("Epoch [number]",fontsize=14)
        ax.set_ylabel("Residuals [min]",fontsize=14)
        ax.grid(True, ls='--')
        plt.tight_layout()
        return fig

    def plot_oc_eclipses(self):

        # O-C plot eclipse times
        fig,ax = plt.subplots(1, figsize=(10,7))

        # predict mid-eclipse time
        emid = self.parameters['emid']

        self.ephemeris['emid_orbit'] = np.round((self.ephemeris['emid'][:,0]-emid)/self.rv_data[0]['priors']['per'])

        # predict mid-transit time
        emid_pred = self.ephemeris['emid_orbit']*self.rv_data[0]['priors']['per'] + emid
        emid_residual = self.ephemeris['emid'][:,0] - emid_pred
        
        # eclipse_mid_time(
        #     self.rv_data[0]['priors']['per'], 
        #     self.rv_data[0]['priors']['ars'], 
        #     self.rv_data[0]['priors']['ecc'], 
        #     self.lc_data[0]['priors']['inc'], 
        #     self.rv_data[0]['priors']['omega'], 
        #     self.lc_data[0]['priors']['tmid'] - self.lc_data[0]['priors']['per'])

        self.ephemeris['emid_orbit'] = np.round((self.ephemeris['emid'][:,0]-emid)/self.rv_data[0]['priors']['per'])

        emid_pred = self.ephemeris['emid_orbit']*self.rv_data[0]['priors']['per'] + emid
        emid_residual = self.ephemeris['emid'][:,0] - emid_pred

        print("Eclipse Residual [min]:", emid_residual*24*60)
        
        ax.errorbar(self.ephemeris['emid_orbit'], emid_residual*24*60, yerr=self.ephemeris['emid'][:,1]*24*60, ls='none', marker='s',label='Mid-Eclipse Measurement',color='black')
        ylower = (emid_residual.mean()-3*np.std(emid_residual))*24*60
        yupper = (emid_residual.mean()+3*np.std(emid_residual))*24*60

        # upsample data
        epochs = (np.linspace(self.ephemeris['emid_orbit'].min()-7, max(7,self.ephemeris['emid_orbit'].max()+7), 1000))

        depoch = epochs.max() - epochs.min()
        ax.set_xlim([epochs.min()-depoch*0.01, epochs.max()+depoch*0.01])

        # best fit solution
        model = epochs*self.parameters['per'] + emid

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['per'], self.errors['per'], size=10000)
        mc_b = np.random.normal(self.parameters['emid'], self.errors['emid'], size=10000)
        mc_model = np.expand_dims(epochs,-1) * mc_m + mc_b

        # create a fill between area for uncertainty of new ephemeris
        diff = mc_model.T - model
        ax.fill_between(epochs, np.percentile(diff,16,axis=0)*24*60, np.percentile(diff,85,axis=0)*24*60, alpha=0.2, color='k', label=r'Uncertainty ($\pm$ 1$\sigma$)')

        # duplicate axis and plot days since mid-transit
        ax2 = ax.twiny()
        ax2.set_xlabel(f"Time [BJD - {emid:.1f}]",fontsize=14)
        ax2.set_xlim(ax.get_xlim())
        xticks = ax.get_xticks()
        dt = np.round(xticks*self.parameters['per'],1)
        ax2.set_xticklabels(dt)
        show_2sigma = False

        if self.ephemeris['prior'] is not None:
            # create fill between area for uncertainty of old/prior ephemeris
            epochs_p = ((epochs*self.rv_data[0]['priors']['per'] + emid) - self.ephemeris['prior']['emid'][0])/self.ephemeris['prior']['per'][0]
            #epochs_p = (np.linspace(self.ephemeris['tmid'][:,0].min()-7, self.ephemeris['tmid'][:,0].max()+7, 1000) 
            prior_p = epochs_p*self.ephemeris['prior']['per'][0] + self.ephemeris['prior']['tmid'][0]
            mc_m_p = np.random.normal(self.ephemeris['prior']['per'][0], self.ephemeris['prior']['per'][1], size=10000)
            mc_b_p = np.random.normal(self.ephemeris['prior']['emid'][0], self.ephemeris['prior']['emid'][1], size=10000)
            mc_model_p = np.expand_dims(epochs_p,-1) * mc_m_p + mc_b_p
            diff_p = mc_model_p.T - model

            # plot an invisible line so the 2nd axes are happy
            ax2.plot(epochs, (model-prior_p)*24*60, ls='--', color='r', alpha=0)

            if show_2sigma:
                ax.fill_between(epochs, np.percentile(diff_p,2,axis=0)*24*60, np.percentile(diff_p,98,axis=0)*24*60, alpha=0.1, color='r', label=r'Prior ($\pm$ 2$\sigma$)')
            else:
                # show ~1 sigma
                ax.fill_between(epochs, np.percentile(diff_p,16,axis=0)*24*60, np.percentile(diff_p,85,axis=0)*24*60, alpha=0.1, color='r', label=r'Prior ($\pm$ 1$\sigma$)')

            #if ylim == 'prior':
            #    ax.set_ylim([ min(np.percentile(diff_p,1,axis=0)*24*60),
            #                max(np.percentile(diff_p,99,axis=0)*24*60)])
            ax.set_ylim([ 0.5*(min(np.percentile(diff,1,axis=0)*24*60) + min(np.percentile(diff_p,1,axis=0)*24*60)),
                            0.5*(max(np.percentile(diff,99,axis=0)*24*60) + max(np.percentile(diff_p,99,axis=0)*24*60))])

        ax.axhline(0,color='black',alpha=0.5,ls='--',
                    label="Period: {:.5f}+-{:.5f} days\nE_mid: {:.4f}+-{:.4f} BJD".format(self.parameters['per'], self.errors['per'], np.round(emid,4), np.round(self.errors['emid'],4)))

        # TODO sig figs
        #lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
        #    str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
        #    str(round_to_2(self.errors.get('tmid',0)))
        #)

        ax.legend(loc='best')
        ax.set_xlabel("Epoch [number]",fontsize=14)
        ax.set_ylabel("Residuals [min]",fontsize=14)
        ax.grid(True, ls='--')
        plt.tight_layout()
        return fig