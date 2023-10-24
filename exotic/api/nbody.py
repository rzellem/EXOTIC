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
# ########################################################################### #
# Code for running N-body simulations and computing transit timing variations.
#  Adapted from: https://github.com/pearsonkyle/Nbody-ai
# ########################################################################### #
# ########################################################################### #
import os
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
import rebound
from exotic.api.plotting import corner
from ultranest import ReactiveNestedSampler
from astropy.io import fits
from astropy import units as u
from scipy.interpolate import interp1d
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

mearth = u.M_earth.to(u.kg)
msun = u.M_sun.to(u.kg)
mjup = u.M_jup.to(u.kg)

# outlier robust "average"
maxavg = lambda x: (np.percentile(np.abs(x),75)+np.max(np.abs(x))  )*0.5

def empty_data(N): # orbit parameters in Nbody simulation
    return {
        'x':np.zeros(N),
        'y':np.zeros(N),
        'z':np.zeros(N),
        'P':np.zeros(N),
        'a':np.zeros(N),
        'e':np.zeros(N),
        'inc':np.zeros(N),
        'Omega':np.zeros(N),
        'omega':np.zeros(N),
        'M':np.zeros(N),
    }

def generate(objects, Ndays=None, Noutputs=None):
    # create rebound simulation 
    # for object parameters see: 
    # https://rebound.readthedocs.io/en/latest/_modules/rebound/particle.html
    sim = rebound.Simulation()
    sim.units = ('day', 'AU', 'Msun')
    for i in range(len(objects)):
        sim.add( **objects[i] ) 
    sim.move_to_com() 

    if Ndays and Noutputs:
        return integrate(sim, objects, Ndays, Noutputs)
    else:
        return sim

def integrate(sim, objects, Ndays, Noutputs):
    # ps is now an array of pointers and will change as the simulation runs
    ps = sim.particles
    
    times = np.linspace(0., Ndays, Noutputs) # units of days 

    pdata = [empty_data(Noutputs) for i in range(len(objects)-1)]
    star = {'x':np.zeros(Noutputs),'y':np.zeros(Noutputs),'z':np.zeros(Noutputs) }

    # run simulation time steps 
    for i,time in enumerate(times):
        sim.integrate(time)

        # record data 
        for k in star.keys():
            star[k][i] = getattr(ps[0], k)

        for j in range(1,len(objects)):
            for k in pdata[j-1].keys():
                pdata[j-1][k][i] = getattr(ps[j],k)

    sim_data = ({
        'pdata':pdata,
        'star':star,
        'times':times,
        'objects':objects,
        'dt': Noutputs/(Ndays*24*60*60) # conversion factor to get seconds for RV semi-amp
    })

    return sim_data 

def analyze(m, ttvfast=False):

    if ttvfast:
        tt = transit_times( m['pdata'][0]['x'], m['star']['x'], m['times'] )
        ttv,per,t0 = TTV(np.arange(len(tt)),tt )
        return np.arange(len(ttv)),ttv,tt

    RV = np.diff(m['star']['x'])*1.496e11*m['dt'] # measurements -> seconds

    freq,power,peak_periods = lomb_scargle(m['times'][1:], RV, min_freq=1./365, max_freq=1.)

    data = {
        'times':m['times'],
        'RV':{'freq':freq, 'power':power, 'signal':RV, 'max':maxavg(RV), 'peak_periods':peak_periods},
        'mstar': m['objects'][0]['m'],
        'planets':[],
        'objects':m['objects']
    }

    # parse planet data 
    for j in range(1,len(m['objects'])):
        key = label='Planet {}'.format(j)
        pdata = {}

        # compute transit times 
        tt = transit_times( m['pdata'][j-1]['x'], m['star']['x'], m['times'] )
        if len(tt)>=3:
            ttv,per,t0 = TTV(np.arange(len(tt)),tt )
        else:
            per = m['objects'][j]['P']
            t0 = 0 
            ttv = [0]

        # periodogram of o-c                
        freq,power,peak_periods = lomb_scargle(np.arange(len(ttv)),ttv)

        # save data 
        for k in ['e','inc','a', 'omega']:
            pdata[k] = np.mean( m['pdata'][j-1][k] )
        pdata['mass'] = m['objects'][j]['m']
        pdata['P'] = per
        pdata['t0'] = t0
        pdata['tt'] = tt
        pdata['ttv'] = ttv
        pdata['max'] = maxavg(ttv)
        pdata['freq'] = freq
        pdata['power'] = power
        pdata['peak_periods'] = peak_periods
        
        # down sample data for plotting
        pdata['x'] = m['pdata'][j-1]['x'][::4]
        pdata['y'] = m['pdata'][j-1]['z'][::4]
        data['planets'].append(pdata)

    return data

def lomb_scargle(t,y,dy=None, min_freq=None, max_freq=None, npeaks=0, peaktol=0.05):

    # perform quick periodicity analysis on residuals + find peaks in power plot
    if dy is not None:
        ls = LombScargle(t, y, dy)
    else:
        ls = LombScargle(t, y)

    if min_freq is None:
        min_freq = 1./(1.5*(max(t)-min(t)))
    
    if max_freq is None:
        max_freq = 1./2.5

    # compute power spectrum
    freq,power = ls.autopower(maximum_frequency=max_freq, minimum_frequency=min_freq, 
                              nyquist_factor=2, samples_per_peak=5, method='cython')

    # find peaks in power spectrum
    peaks,amps = find_peaks(power,height=peaktol)

    # sort by biggest to smallest amplitude
    peaks = peaks[np.argsort(amps['peak_heights'])[::-1]]
    peak_periods = 1./freq[peaks]
    return freq,power,peak_periods

def report(data, savefile=None):

    # set up simulation summary report 
    f = plt.figure( figsize=(12,8) ) 
    plt.subplots_adjust()
    ax = [  plt.subplot2grid( (2,3), (0,0) ), # x,y plot
            plt.subplot2grid( (2,3), (1,0) ), # table data 
            plt.subplot2grid( (2,3), (0,1) ), # # RV semi-amplitude plot  #include theoretical limit for 2body system
            plt.subplot2grid( (2,3), (1,1) ), # O-C plot # TODO add second axis for time (day)
            plt.subplot2grid( (2,3), (0,2) ), # lomb scargle for RV semi-amplitude
            plt.subplot2grid( (2,3), (1,2) )  # lomb scargle for o-c
        ]
    plt.subplots_adjust(top=0.96, bottom=0.07, left=0.1, right=0.98, hspace=0.3, wspace=0.3)

    ax[2].plot(data['times'][1:],data['RV']['signal'],'k-' )
    ax[2].set_xlim([0, 2.5*data['planets'][-1]['P'] ])
    ax[2].set_ylabel('RV semi-amplitude (m/s)')
    ax[2].set_xlabel('time (day)')

    ax[4].plot(1./data['RV']['freq'],data['RV']['power'] )
    # plot reconstructed signal from sin series 

    ax[4].set_xlabel('Period (day)')
    ax[4].set_ylabel('|RV semi-amplitude| Power')
    ax[4].set_xlim([1, 1.5*data['planets'][-1]['P'] ])
    #u = (m['objects'][1]['m']/m['objects'][0]['m']) * np.sqrt( G*(m['objects'][0]['m']+m['objects'][1]['m'])/m['objects'][1]['a']  )*1.496e11/(24*60*60)
    #print('expected RV:',u)

    # create table stuff 
    keys = ['mass','a','P','inc','e','max']
    units = [msun/mearth,1,1,1,1,24*60]
    rounds = [2,3,2,2,3,1,1]

    # make 2d list for table 
    tdata=[]
    for i in range(len(units)):
        tdata.append( [0 for i in range(len(data['planets']))] )
    
    tinfo = { 
        'mass':'Mass (Mearth)',
        'a':'Semi-major Axis (au)',
        'P':'Period (day)', 
        'inc':'Inclination (deg)',
        'e':'Eccentricity',
        'max':'TTV Max Avg (min)',
    }
    row_labels = [ tinfo[k] for k in keys ]
    col_labels = [ "P {}".format(i+1) for i in range(len(data['planets']))]

    # for each planet in the system 
    for j in range(1,len(data['objects'])):
        
        # plot orbits
        ax[0].plot( data['planets'][j-1]['x'], data['planets'][j-1]['y'],label='Planet {}'.format(j),lw=0.5,alpha=0.5 )

        if len(data['planets'][j-1]['ttv']) >= 2:
            # plot O-C
            ax[3].plot(data['planets'][j-1]['ttv']*24*60,label='Planet {}'.format(j) )  # convert days to minutes

            # lomb-scargle for O-C
            ax[5].semilogx( 1./data['planets'][j-1]['freq'], data['planets'][j-1]['power'], label='Planet {}'.format(j) )
            
            # add peak periods to legend
            for i,per in enumerate(data['planets'][j-1]['peak_periods']):
                # vertical line with label
                if per > 0:
                    ax[5].axvline(per, ls='--', lw=0.5, label=f"{np.round(per,2)} day")

        # populate table data 
        for i,k in enumerate(keys):
            tdata[i][j-1] = np.round( data['planets'][j-1][k] * units[i], rounds[i])

            if k == 'inc':
                tdata[i][j-1] = np.round( np.rad2deg(data['planets'][j-1][k]) * units[i], rounds[i])

    table = ax[1].table(cellText=tdata, 
                        colLabels=col_labels,
                        rowLabels=row_labels,
                        colWidths=[0.1]*len(col_labels),
                        loc='center')

    ax[1].set_title('Stellar mass: {:.2f} Msun'.format(data['mstar']) )
    table.scale(1.5,1.5)
    
    ax[1].axis('tight')
    ax[1].axis('off')

    ax[0].set_ylabel('(au)')
    ax[0].set_xlabel('(au)')
    ax[0].legend(loc='best')

    ax[1].set_ylabel('position (au)')
    ax[1].set_xlabel('time (day)')

    ax[3].set_ylabel('O-C (minutes)')
    ax[3].set_xlabel('transit epoch')
    ax[3].legend(loc='best')
    ax[3].grid(True)
        
    ax[5].set_xlabel('Period (epoch)')
    ax[5].set_ylabel('|O-C| Power' )
    #ax[5].set_xlim([1,30])
    ax[5].legend(loc='best')

    if savefile:
        plt.savefig(savefile)
        plt.close()
    else:
        plt.show()

def find_zero(t1,x1, t2,x2):
    """
    Find the zero of a line given two points.

    Parameters
    ----------
    t1 : float
        Time of first point.
    
    x1 : float
        Position of first point.

    t2 : float
        Time of second point.
    
    x2 : float
        Position of second point.

    Returns
    -------
    T0 : float
        Time of zero.
    """
    m = (x2-x1)/(t2-t1)
    T0 = -x1/m + t1
    return T0

def transit_times(xp,xs,times):
    """
    Find transit times from position data.

    Parameters
    ----------
    xp : array
        Planet x positions.

    xs : array
        Star x positions.

    times : array
        Simulation times.

    Returns
    -------
    tt : array
        Transit times.
    """
    # check for sign change in position difference
    dx = xp-xs 
    tt = []
    for i in range(1,len(dx)):
        if dx[i-1] >= 0 and dx[i] <= 0:
            tt.append( find_zero(times[i-1],dx[i-1], times[i],dx[i]) )
    return np.array(tt)

def TTV(epochs, tt):
    """
    Fit a line to the data and subtract it from the data.

    Parameters
    ----------
    epochs : array
        Orbits epochs of the data.

    tt : array
        Transit times.

    Returns
    -------
    ttv : array
        Transit times with the linear trend subtracted.

    m : float
        Slope of the linear trend.

    b : float
        Intercept of the linear trend.
    """
    N = len(epochs)
    A = np.vstack([np.ones(N), epochs]).T
    b, m = np.linalg.lstsq(A, tt, rcond=None)[0]
    ttv = (tt-m*np.array(epochs)-b)
    return [ttv,m,b]


# directory path of file
fpath = os.path.dirname(os.path.abspath(__file__))
def estimate_prior( amplitude, periodicity, amp_lim=2, per_lim=4, file_name = os.path.join(fpath,"ttv_grid.fits"), show=True):

    with fits.open(file_name) as ttv_grid_list:

        # The grid is stored in multiple extensions
        Npts = ttv_grid_list[0].header['NAXIS1']

        # in Earth
        mass_grid = np.linspace(ttv_grid_list[0].header['CRVAL1'], ttv_grid_list[0].header['CRVAL1'] + Npts*ttv_grid_list[0].header['CDELT1'], Npts)
        mass_grid_2 = np.linspace(ttv_grid_list[0].header['CRVAL1'], ttv_grid_list[0].header['CRVAL1'] + Npts*ttv_grid_list[0].header['CDELT1'], int(0.5*Npts))
        period_grid = np.linspace(ttv_grid_list[0].header['CRVAL2'], ttv_grid_list[0].header['CRVAL2'] + Npts*ttv_grid_list[0].header['CDELT2'], Npts)
        #period_grid_2 = np.linspace(ttv_grid_list[0].header['CRVAL2'], ttv_grid_list[0].header['CRVAL2'] + Npts*ttv_grid_list[0].header['CDELT2'], int(0.5*Npts))

        # extract data
        pg, mg = np.meshgrid(period_grid, mass_grid)
        amplitude_grid = ttv_grid_list[0].data
        periodicity1_grid = ttv_grid_list[1].data
        periodicity2_grid = ttv_grid_list[2].data
        # combine periods by min
        #periodicity_grid = np.minimum(periodicity1_grid, periodicity2_grid)
        rv_grid = ttv_grid_list[3].data

        # amplitude grid
        fig, ax = plt.subplots(3, 2, figsize=(9, 11))
        fig.suptitle("N-body Estimates from Transit Timing Variations", fontsize=16)
        plt.subplots_adjust(left=0.0485, right=0.985, bottom=0.15, top=0.9, wspace=0.3)
        im = ax[0,0].imshow(amplitude_grid, origin='lower', extent=[mass_grid[0], mass_grid[-1],period_grid[0], period_grid[-1]], vmin=0,vmax=30, aspect='auto', cmap='jet', interpolation='none')
        ax[0,0].set_ylabel('Period Ratio')
        ax[0,0].set_xlabel('Mass [Earth]')
        cbar = fig.colorbar(im, ax=ax[0,0])
        cbar.set_label('TTV Amplitude [min]')

        # periodoicity 1 grid
        im = ax[0,1].imshow(periodicity1_grid, origin='lower', extent=[mass_grid[0], mass_grid[-1],period_grid[0], period_grid[-1]], vmin=0,vmax=0.5*np.nanmax(ttv_grid_list[1].data), aspect='auto', cmap='jet', interpolation='none')
        ax[0,1].set_ylabel('Period Ratio')
        ax[0,1].set_xlabel('Mass [Earth]')
        cbar = fig.colorbar(im, ax=ax[0,1])
        cbar.set_label('TTV Periodicity 1st Order [epoch]')

        # simulate some results based on the periodogram
        mask = (amplitude_grid > amplitude-amp_lim) & (amplitude_grid < amplitude+amp_lim) & \
                (((periodicity1_grid > periodicity-per_lim) & (periodicity1_grid < periodicity+per_lim)) | \
                 ((periodicity2_grid > periodicity-per_lim) & (periodicity2_grid < periodicity+per_lim)))
        # plot mask
        im = ax[1,0].imshow(mask, origin='lower', extent=[mass_grid[0], mass_grid[-1],period_grid[0], period_grid[-1]], vmin=0,vmax=1, aspect='auto', cmap='binary_r', interpolation='none')
        ax[1,0].set_ylabel('Period Ratio')
        ax[1,0].set_xlabel('Mass [Earth]')
        cbar = fig.colorbar(im, ax=ax[1,0])
        cbar.set_label('N-body Prior from O-C Data')

        # compare to original
        # TODO find why the last bin has double the number of points
        masses = mg.T[mask].flatten()
        ax[1,1].hist(masses, bins=mass_grid_2, alpha=0.5)
        ax[1,1].set_xlabel('Mass [Earth]')
        ax[1,1].set_ylabel('Prior Probability')
        ax[1,1].axes.yaxis.set_ticklabels([])
        ax[1,1].set_xlim([mg.min(), mg.max()-10])

        # compare to original
        periods = pg.T[mask].flatten()
        ax[2,0].hist(periods, bins=period_grid, density=True, alpha=0.5)
        ax[2,0].set_xlabel('Period Ratio')
        ax[2,0].set_ylabel('Prior Probability')
        ax[2,0].axes.yaxis.set_ticklabels([])
        ax[2,0].set_xlim([pg.min(), pg.max()])

        # plot histogram for rv data
        rvs = rv_grid.T[mask].flatten()
        ax[2,1].hist(rvs, bins=50, density=True, alpha=0.5)
        ax[2,1].set_xlabel('RV Semi-Amplitude [m/s]')
        ax[2,1].set_ylabel('Prior Probability')
        ax[2,1].axes.yaxis.set_ticklabels([])

    # m_earth, days, m/s
    return masses, periods, rvs, fig, ax

def interp_distribution(values, nbins=50):
    value_grid = np.linspace(np.min(values), np.max(values), nbins)
    heights, edges = np.histogram(values, bins=value_grid, density=True)
    edge_center = (edges[:-1] + edges[1:]) / 2

    # Step 2: Normalize the histogram
    bin_widths = np.diff(edges)
    total_area = np.sum(bin_widths * heights)
    normalized_heights = heights / total_area

    # Step 3: Interpolate to create a continuous PDF
    return interp1d(edge_center, normalized_heights, kind='linear', bounds_error=False, fill_value=0)




class nbody_fitter():

    def __init__(self, data, prior=None, bounds=None, verbose=True):
        self.data = data
        self.bounds = bounds
        self.prior = prior
        self.verbose = verbose
        self.fit_nested()

    def fit_nested(self):

        # set up some arrays for mapping sampler output
        freekeys = []
        boundarray = []
        for i,planet in enumerate(self.bounds):
          for bound in planet:
            if '_logl' in bound or '_fn' in bound:
                continue
            freekeys.append(f"{i}_{bound}")
            boundarray.append(planet[bound])

        # find min and max time for simulation
        min_time = np.min(self.data[1]['Tc'])
        max_time = np.max(self.data[1]['Tc'])
        # todo extend for multiplanet systems

        sim_time = (max_time-min_time)*1.05
        # TODO extend for multiplanet systems
        Tc_norm = self.data[1]['Tc'] - min_time  # normalize the data to the first observation
        self.orbit = np.rint(Tc_norm / self.prior[1]['P']).astype(int)  # number of orbits since first observation (rounded to nearest integer)

        # numpify
        boundarray = np.array(boundarray)
        bounddiff = np.diff(boundarray,1).reshape(-1)

        # create queue and save simulation to
        def loglike(pars):
            chi2 = 0

            # set parameters
            for i,par in enumerate(pars):
                idx,key = freekeys[i].split('_')
                idx = int(idx)
                if key == 'tmid':
                    continue
                # this dict goes to REBOUND and needs specific keys
                self.prior[idx][key] = par

            # prior likelihood function for P, m of outer planet
            likelihood_m = self.bounds[2]['m_logl'](self.prior[2]['m'])
            likelihood_P = self.bounds[2]['P_logl'](self.prior[2]['P'])
            if likelihood_m <= 0 or likelihood_P <= 0:
                return -1e6

            # take difference between data and simulation
            # run N-body simulation
            sim_data = generate(self.prior, sim_time, int(sim_time*24)) # uses linspace behind the scenes

            # json structure with analytics of interest from the simulation
            # ttv_data = analyze(sim_data) # slow

            sim_shift = 0
            # loop over planets, check for data
            for i,planet in enumerate(self.prior):
                if self.data[i]:
                    # compute transit times from N-body simulation
                    Tc_sim = transit_times( sim_data['pdata'][i-1]['x'], sim_data['star']['x'], sim_data['times'] )

                    # derive an offset in time from the first planet
                    if i-1==0:
                        sim_shift = Tc_sim.min()

                    # shift the first mid-transit in simulation to observation
                    Tc_sim -= sim_shift # first orbit has tmid at 0

                    # scale Tc_sim to data
                    try:
                        residual = self.data[i]['Tc'] - Tc_sim[self.orbit]
                    except:
                        #import pdb; pdb.set_trace
                        #ttv,m,b = TTV(np.arange(Tc_sim.shape[0]), Tc_sim)
                        #ttv1,m1,b1 = TTV(self.orbit, self.data[i]['Tc'])
                        # could be unstable orbit or not enough data
                        # switch to average error and clip by max epoch?
                        # print(self.prior)
                        chi2 += -1e6
                        continue

                    import pdb; pdb.set_trace
                    Tc_sim += residual.mean() # why?

                    # take difference between data and simulation
                    try:
                        chi2 += -0.5*np.sum(((self.data[i]['Tc'] - Tc_sim[self.orbit])/self.data[i]['Tc_err'])**2)
                    except:
                        chi2 += -1e6
                        # print(self.prior)
                        # usually unstable orbit

            return chi2


        def prior_transform(upars):
            return (boundarray[:,0] + bounddiff*upars)

        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e5)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e5, show_status=self.verbose, 
viz_callback=self.verbose)

        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        # TODO finish saving final results + posteriors

        #for i, key in enumerate(freekeys):
        #    self.parameters[key] = self.results['maximum_likelihood']['point'][i]
        #    self.errors[key] = self.results['posterior']['stdev'][i]
        #    self.quantiles[key] = [
        #        self.results['posterior']['errlo'][i],
        #        self.results['posterior']['errup'][i]]


if __name__ == '__main__':

    # create some sample data
    objects = [
        # units: Msun, Days, au
        {'m':0.95}, # stellar mass
        {'m':1.169*mjup/msun, 'P':2.797436, 'inc':3.14159/2, 'e':0, 'omega':0 }, 
        {'m':0.1*mjup/msun, 'P':2.797436*1.9, 'inc':3.14159/2, 'e':0.0,  'omega':0  }, 
    ] # HAT-P-37

    # create REBOUND simulation
    n_orbits = 2000

    # time the simulation
    t1 = time.time()
    # inputs: object dict, length of simulation in days, number of timesteps [1hr] (should be at least 1/20 orbital period)
    sim_data = generate(objects, objects[1]['P']*n_orbits, int(n_orbits*objects[1]['P']*24) )
    t2 = time.time()
    print(f"Simulation time: {t2-t1:.2f} seconds")

    # collect the analytics of interest from the simulation
    # lomb-scargle can be a lil slow
    ttv_data = analyze(sim_data)

    # plot the results
    report(ttv_data)

    # create a fake dataset
    tmids = 2459150 + ttv_data['planets'][0]['tt']
    # add random noise to observations
    tmids += np.random.normal(0,0.5,len(tmids))/(24*60)
    # randomly select 50 observations without repeat
    tmids = np.random.choice(tmids,50,replace=False)
    # add random error to observations between
    err = 1/24/60 + np.random.random()*0.25/24/60 + np.random.normal(0,0.1,len(tmids))/(24*60)
    # estimate orbital epochs
    orbit = np.rint((tmids-tmids.min())/ttv_data['planets'][0]['P']).astype(int)

    # estimate period from data
    ttv,P,T0 = TTV(orbit, tmids)

    # run though linear fitter to estimate prior
    from exotic.api.ephemeris import ephemeris_fitter
    
    # min and max values to search between for fitting
    bounds = {
        'P': [P - 0.1, P + 0.1],    # orbital period
        'T0': [T0 - 0.1, T0 + 0.1]  # mid-transit time
    }

    # used to plot red overlay in O-C figure
    prior = {
        'P': [P, 0.00001],   # value from linear lstq
        'T0': [T0, 0.001]  # value from linear lstq
    }

    lf = ephemeris_fitter(tmids, err, bounds, prior=prior)

    fig,ax = lf.plot_oc()
    plt.tight_layout()
    plt.show()

    # search for periodic signals in the data
    fig,ax = lf.plot_periodogram()
    plt.tight_layout()
    plt.savefig('periodogram.png')
    plt.show()
    plt.close()

    # estimate ttv amplitude
    amp = lf.amplitudes[0]*24*60
    per = lf.best_periods[0] # 1st order solution

    # estimate prior using periods from linear fit periodogram
    masses, per_ratio, rvs, fig, ax = estimate_prior(amp, per)
    masses *= mearth/msun
    periods = per_ratio*lf.parameters['P']
    prior_fn_mass = interp_distribution(masses)
    prior_fn_per = interp_distribution(periods)
    plt.tight_layout()
    plt.savefig('ttv_prior.png')
    plt.show()
    plt.close()

    # Parameters for N-body Retrieval
    nbody_prior = [
        # star
        {'m':0.95},

        # inner planet
        {'m':1.169*mjup/msun,
        'P':lf.parameters['P'],
        'inc':3.14159/2,
        'e':0,
        'omega':0},

        # outer planet
        {'m':masses.mean(),
        'P':periods.mean(),
        'inc':3.14159/2,
        'e':0,
        'omega':0,},
    ]

    # specify data to fit
    data = [
        {},                            # data for star (e.g. RV)
        {'Tc':tmids, 'Tc_err':err},    # data for inner planet (e.g. Mid-transit times)
        {}                             # data for outer planet (e.g. Mid-transit times)
    ]

    # TODO break P,m modes into individual runs

    # set up where to look for a solution
    nbody_bounds = [
        {}, # no bounds on stellar parameters

        {   # bounds for inner planet
            'P': [nbody_prior[1]['P']-0.025, nbody_prior[1]['P']+0.025],  # based on solution from linear fit\
            # 'T0': [None, None]  # mid-transit is automatically adjusted, no need to include in bounds
        },
        {  # bounds for outer planet
            'P':[periods.min(), periods.max()], # period [day]
            'P_logl': prior_fn_per,             # prior likelihood function
            'm':[masses.min(), masses.max()],   # mass [msun
            'm_logl': prior_fn_mass,            # prior likelihood function
            'omega': [-np.pi, np.pi]            # argument of periastron [rad]
        }
    ]

    # run the fitter
    nfit = nbody_fitter(data, nbody_prior, nbody_bounds)

    # print(nfit.parameters)
    # print(nfit.errors)