import numpy as np
from pandas import read_csv
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from exotic.api.rv_fitter import rv_fitter

Mjup = const.M_jup.to(u.kg).value
Msun = const.M_sun.to(u.kg).value
Rsun = const.R_sun.to(u.m).value
Grav = const.G.to(u.m**3/u.kg/u.day**2).value

if __name__ == "__main__":
    df = read_csv('HD80606_RV.csv')
    # columns for: BJD, Vel(m/s), ErrVel(m/s), Telescope

    # format keys for input
    prior = {
        'rprs':0.1,
        'per':111.4367,
        'inc':89.269, # https://ui.adsabs.harvard.edu/abs/2017AJ....153..136S/abstract
        'u0': 0.49387060813646527, 'u1': -0.07294561715247563,  # limb darkening (nonlinear)
        'u2': 0.4578497817617948, 'u3': -0.21582471000247333,   # TESS bandpass generated from exotethys
        'ecc': 0.93132,
        'omega':-58.699,
        'tmid':2459556.6942,
        'a1':1, # transit airmass - not used
        'a2':0,
        'fpfs':0.5, # F_p/F_s - for eclipse depth
        #'mu':((mplanet*const.M_jup) / (mplanet*const.M_jup + mstar*const.M_sun)).value,
        'rstar':1.066, # R_sun
        'mstar':1.05,  # M_Sun
        'mplanet':4.20, # M_Jupiter
        'rv_linear':0,
        'rv_quad':0
    }
    # estimate some ratios and semi-major axis
    prior['mu'] = prior['mplanet']*Mjup / (prior['mplanet']*Mjup + prior['mstar']*Msun)
    mtotal = Msun*prior['mstar'] + Mjup*prior['mplanet'] # kg

    # semi-major axis using Kepler's 3rd law
    semimajor = (Grav*mtotal*prior['per']**2/4/np.pi**2)**(1/3) # m

    # estimate a/Rs
    prior['ars'] = semimajor/(prior['rstar']*Rsun)

    # alloc data
    rv_data = []
    local_rv_bounds = []

    # loop over telescopes
    for tele in df.Telescope.unique():
        # get data for this telescope
        df_tel = df[df.Telescope == tele]

        # add to data
        rv_data.append({
            'time':df_tel['BJD'].values,
            'vel':df_tel['Vel(m/s)'].values,
            'velerr':df_tel['ErrVel(m/s)'].values,
            'priors':prior.copy(),
            'name':tele
        })

        # local bounds are applied to each dataset separately
        print(f"{tele} has {len(rv_data[-1]['time'])} points")
        local_rv_bounds.append({
            #"jitter":[0,10], # don't fit for this, too degenerate
        })

    # bounds for optimization
    global_bounds = {
        'per':[111.3,111.5],
        'omega':[-65,-55],
        'ecc':[0.92,0.94],
        #'rv_linear':[-0.01,0.01], # m/s/day
        'mplanet':[3,4.5], # M_Jupiter
    }

    myfit = rv_fitter(rv_data, global_bounds, local_rv_bounds, verbose=True)
    
    # print Bayesian evidence
    print(myfit.results['logz'], myfit.results['logzerr'])

    # corner plot
    myfit.plot_triangle()
    plt.tight_layout()
    plt.savefig('RV_triangle.png')
    plt.show()

    # best fit
    myfit.plot_bestfit()
    plt.tight_layout()
    plt.savefig('RV_bestfit.png')
    plt.show()

    # show acceleration, useful for exposure time calculations
    # can occasionally show jitter in solution due to ODE solver
    # fig,axs = myfit.plot_bestfit_acceleration(phase_limits=[-0.06,-0.02])
    # plt.tight_layout()
    # plt.savefig('RV_timeseries_acceleration.png')
    # plt.show()