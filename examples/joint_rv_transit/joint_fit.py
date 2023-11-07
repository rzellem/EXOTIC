import pickle
import numpy as np
from pandas import read_csv
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from exotic.api.rv_fitter import rv_fitter
from exotic.api.joint_fitter import joint_fitter
from pylightcurve.models.exoplanet_lc import eclipse_mid_time

Mjup = const.M_jup.to(u.kg).value
Msun = const.M_sun.to(u.kg).value
Rsun = const.R_sun.to(u.m).value
Grav = const.G.to(u.m**3/u.kg/u.day**2).value


if __name__ == "__main__":
    df = read_csv('HD80606_RV.csv')

    # format keys for input
    prior = {
        'rprs':0.1,
        'per':111.4367,
        'inc':89.269, # https://ui.adsabs.harvard.edu/abs/2017AJ....153..136S/abstract
        'u0': 0.49387060813646527, 'u1': -0.07294561715247563,  # exotethys - limb darkening (nonlinear)
        'u2': 0.4578497817617948, 'u3': -0.21582471000247333,   # TESS bandpass
        'ecc': 0.93132,
        'omega':-58.699,
        'tmid':2458888.0757, # TESS
        'a1':1, # transit airmass - not used
        'a2':0,
        'fpfs':0.5, # F_p/F_s - for eclipse depth
        #'mu':((mplanet*const.M_jup) / (mplanet*const.M_jup + mstar*const.M_sun)).value,
        'rstar':1.05, # R_sun
        'mstar':1.05,  # M_Sun
        'mplanet':4.20, # M_Jupiter
        'rv_linear':0,
        'rv_quad':0
    }

    # historical values
    ephemeris = {
        'tmid': np.array([
            [2455210.6420, 0.004],  # Hebrard et al 2010 https://www.aanda.org/articles/aa/pdf/2010/08/aa14327-10.pdf
            [2454987.7842, 0.0049],  # Winn et al. 2010 https://arxiv.org/pdf/0907.5205.pdf 
            [2459556.6942, 0.0035], # 07/08 Dec, exoplanet watch
            [2455210.6502, 0.0064],  # Shporer et al. 2010 https://arxiv.org/pdf/1008.4129.pdf
        ]),
        'emid': np.array([
            [2454424.736, 0.00025] # spitzer - Laughlin2009, smaller uncertainty to give more weight
        ]),
        'historic': np.array([ # get plotted but not used
            [2454876.316, 0.023],  # Pont et al. 2009 https://arxiv.org/pdf/0906.5605.pdf - partial transit
            [2455099.196, 0.026],   # Shporer et al. 2010 https://arxiv.org/pdf/1008.4129.pdf - partial transit
            [2454876.338, 0.017], # OA, exoclock - partial transit
            [2457439.401, 0.012], # BP, exoclock - partial transit
            [2459222.401, 0.016], # RO, exoclock - partial transit
        ]),
        'prior': { 
            'per':[111.43670,0.0004],
            'tmid':[2455210.6428,0.001], # bonomo et al 2017
            'emid':[2455204.7938334458, 0.09805812418414957], # from propagating parameters below
            #'ecc':[0.93226,0.00066],
            #'omega':[301.03-360,0.2],
            #'a':[0.4565,0.0053],
            #'rstar':[1.037,0.032]
        }, 
        'noise': 0,  # days - additive noise for each measurement
    }

    # check emid key exists
    if 'emid' not in ephemeris.keys():
        # estimate the historical eclipse time
        emids = []
        default = prior.copy()
        print("computing mid-eclipse prior + uncertainty...")
        for i in range(4444):
            for k in ephemeris['prior'].keys():
                default[k] = np.random.normal(ephemeris['prior'][k][0], ephemeris['prior'][k][1])
            default['ars'] = default['a'] * AU / (default['rstar'] * Rsun)
            emid = eclipse_mid_time(default['per'], default['ars'], default['ecc'], default['inc'], default['omega'], default['tmid'])
            emids.append(emid-default['per'])
        ephemeris['prior']['emid'] = [np.mean(emids), np.std(emids)]
        print(ephemeris['prior']['emid'])

    # add a constant noise term to reduce influence
    ephemeris['noise'] = ephemeris['tmid'][:,1].mean()/2
    # my personal belief is that some of tmid values on NEA are underestimated, this helps to reduce their influence

    # compute some ratios
    prior['mu'] = prior['mplanet']*Mjup / (prior['mstar']*Msun)
    mtotal = Msun*prior['mstar'] + \
            Mjup*prior['mplanet'] # kg

    # semi-major axis using Kepler's 3rd law
    semimajor = (Grav*mtotal*prior['per']**2/4/np.pi**2)**(1/3) # m

    # estimate a/Rs
    prior['ars'] = semimajor/(prior['rstar']*Rsun)

    # precompute mid-eclipse time
    emid = eclipse_mid_time(prior['per'], prior['ars'], prior['ecc'], prior['inc'], prior['omega'], prior['tmid'])
    emid -= prior['per'] # subtract one period so it's close to the real tmid

    rv_data = []
    local_rv_bounds = []

    # add RV data and separate into individual datasets
    for tele in df.Telescope.unique():
        
        # special case
        if tele == 'HRS': # break into two distributions
            df_tel = df[df.Telescope == tele]
            errors = df_tel['ErrVel(m/s)'].values # plot histogram to see multiple modes
            pop1 = errors < 5

            rv_data.append({
                'time':df_tel['BJD'].values[pop1],
                'vel':df_tel['Vel(m/s)'].values[pop1],
                'velerr':df_tel['ErrVel(m/s)'].values[pop1],
                'priors':prior.copy(),
                'name':tele+"_1"
            })
            local_rv_bounds.append({}) # not used but required data structure

            rv_data.append({
                'time':df_tel['BJD'].values[~pop1],
                'vel':df_tel['Vel(m/s)'].values[~pop1],
                'velerr':df_tel['ErrVel(m/s)'].values[~pop1],
                'priors':prior.copy(),
                'name':tele+"_2"
            })
            local_rv_bounds.append({}) # not used but required data structure
        else:
                
            df_tel = df[df.Telescope == tele]
            rv_data.append({
                'time':df_tel['BJD'].values,
                'vel':df_tel['Vel(m/s)'].values,
                'velerr':df_tel['ErrVel(m/s)'].values,
                'priors':prior.copy(),
                'name':tele
            })

        #print(f"{tele} has {len(rv_data[-1]['time'])} points")
        local_rv_bounds.append({}) # not used but required data structure, sometimes used for fitting jitter


    # specify bounds for fitting - all uniform distributions
    global_bounds = {
        'per':[111,112],
        'tmid':[prior['tmid']-0.2, prior['tmid']+0.2],

        'inc':[87,90],
        'omega':[-61,-55], # ecosw, esinw will be used after rvfit 
        'ecc':[0.92,0.94],
        'mplanet':[3,4.5], # M_Jupiter, uniform prior
        #'rv_linear':[-0.01,0.01], # m/s/day

        #'mstar':[1.05, 0.001], # M_Sun, makes posteriors too flat
        #'rstar':[1.05, 0.005], # perfectly correlated to mass ratio -> impose gaussian prior
    }

    # set up transit data
    lc_data = []
    local_lc_bounds = []

    # TESS extracted light curve
    tess = np.loadtxt('2458888_08_lightcurve.csv', delimiter=',', skiprows=1)

    # tess data
    lc_data.append({
       'time':tess[:,0],
       'flux':tess[:,1],
       'ferr':tess[:,2],
       'airmass':np.zeros(tess.shape[0]),
       'priors':prior.copy(),
       'name':'TESS'
    })

    # perform a local LC fit to constrain inclination, Rp/Rs and Tmid
    local_lc_bounds.append({'rprs':[0.0,0.2]})

    # perform a global fit to constrain RV
    rv_only_bounds = {
        'per':[111,112],
        'omega':[-61,-55],
        'ecc':[0.92,0.94],
        #'rv_linear':[-0.01,0.01], # m/s/day
        'mplanet':[3,4.5], # M_Jupiter
    }

    # constrain the priors by fitting only the rv data first
    rvfit = rv_fitter(rv_data, rv_only_bounds, local_rv_bounds, verbose=True)

    for key in rvfit.parameters.keys():
        prior[key] = rvfit.parameters[key]

    # update global bounds with +/- 5 sigma from RV only fit
    for key in rvfit.errors:
        if key == 'mu':
            continue
        if key == 'ecc':
            global_bounds[key] = [max(0,rvfit.parameters[key]-5*rvfit.errors[key]), min(1,rvfit.parameters[key]+5*rvfit.errors[key])]
        else:
            global_bounds[key] = [rvfit.parameters[key]-5*rvfit.errors[key], rvfit.parameters[key]+5*rvfit.errors[key]]
    
    # inflate uncertainties for each data set such that the average error is similar to the stdev of residuals
    for i, data in enumerate(rv_data):
        diff = max(0, rvfit.data[i]['residuals'].std() - rvfit.data[i]['velerr'].mean())
        rv_data[i]['velerr'] = rvfit.data[i]['velerr'] + 0.25*diff # should really be 0.5
        print(f"Inflated RV uncertainties for {rv_data[i]['name']} by {diff:2f} m/s")
        rv_data[i]['diff'] = diff

    # update light curve priors based on rv fit
    for lc in lc_data:
        for key in rvfit.parameters:
            lc['priors'][key] = rvfit.parameters[key]        

    # update rv priors based on rv fit
    for rv in rv_data:
        for key in rvfit.parameters:
            rv['priors'][key] = rvfit.parameters[key]

    # run global fit
    myfit = joint_fitter(lc_data, local_lc_bounds, rv_data, local_rv_bounds, 
                        global_bounds, ephemeris,
                        verbose=True, individual_fit=True)

    # save to object to disk
    pickle.dump(myfit, open('grv_fit.pkl', 'wb'))
    
    # make some plots
    myfit.plot_bestfit(phase_limits='none')
    plt.savefig("joint_transit_fit.png")
    plt.close()

    myfit.plot_triangle()
    plt.savefig("joint_posterior.png")
    plt.close()

    myfit.plot_rv_bestfit()
    plt.savefig("joint_rv_fit.png")
    plt.close()

    myfit.plot_oc_transits()
    plt.savefig("joint_oc_transits.png")
    plt.close()

    myfit.plot_oc_eclipses()
    plt.savefig("joint_oc_eclipses.png")
    plt.close()

    for key in myfit.errors:
        print(f"{key} & {myfit.parameters[key]:.5f} $\pm$ {myfit.errors[key]:.5f}")

    for i, data in enumerate(myfit.rv_data):
        print(f"Inflated RV uncertainties for {data['name']} by {data['diff']:2f} m/s")
