from exotic.api.ld import LimbDarkening, test_ld
from exotic.api.elca import transit, lc_fitter
from ldtk.filters import create_tess
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    prior = {
        'rprs': 0.02,                               # Rp/Rs
        'ars': 14.25,                               # a/Rs
        'per': 3.33,                                # Period [day]
        'inc': 88.5,                                # Inclination [deg]
        'u0': 0, 'u1': 0, 'u2': 0, 'u3': 0,         # limb darkening (nonlinear)
        'ecc': 0.5,                                 # Eccentricity
        'omega': 120,                               # Arg of periastron
        'tmid': 0.75,                               # Time of mid transit [day],
        'a1': 50,                                   # Airmass coefficients
        'a2': 0.,                                   # trend = a1 * np.exp(a2 * airmass)

        'T*':5000,
        'FE/H': 0,
        'LOGG': 3.89, 
    }

    """
    # example generating LD coefficients
    tessfilter = create_tess()

    # Test custom-entered ld coefficients
    filter_info = {
        'filter': "TESS",
        'wl_min': tessfilter.wl.min(),
        'wl_max': tessfilter.wl.max(),
        'u0': {"value": None, "uncertainty": None},
        'u1': {"value": None, "uncertainty": None},
        'u2': {"value": None, "uncertainty": None},
        'u3': {"value": None, "uncertainty": None}
    }

    ld_obj = LimbDarkening(prior)
    ld_obj.wl_min = filter_info['wl_min']
    ld_obj.wl_max = filter_info['wl_max']
    test_ld(ld_obj, filter_info)
    ld = [ld_obj.ld0[0], ld_obj.ld1[0], ld_obj.ld2[0], ld_obj.ld3[0]]
    prior['u0'],prior['u1'],prior['u2'],prior['u3'] = ld
    """

    from pylightcurve import exotethys

    # generate limb darkening coefficients for TESS
    get_prior = lambda key: float(prior[key])
    u0,u1,u2,u3 = exotethys(get_prior('LOGG'), get_prior('T*'), get_prior('FE/H'), 'TESS', method='claret', stellar_model='phoenix')
    prior['u0'],prior['u1'],prior['u2'],prior['u3'] = u0,u1,u2,u3

    time = np.linspace(0.7, 0.8, 1000)  # [day]

    # simulate extinction from airmass
    stime = time-time[0]
    alt = 90 * np.cos(4*stime-np.pi/6)
    #airmass = 1./np.cos(np.deg2rad(90-alt))
    airmass = np.zeros(time.shape[0])

    # GENERATE NOISY DATA
    data = transit(time, prior)*prior['a1']*np.exp(prior['a2']*airmass)
    data += np.random.normal(0, prior['a1']*250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time)) + np.random.normal(300e-6, 50e-6, len(time))

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, 0.1],
        'tmid': [prior['tmid']-0.01, prior['tmid']+0.01],
        'ars': [13, 15],
        #'a2': [0, 0.3] # uncomment if you want to fit for airmass
        # never list 'a1' in bounds, it is perfectly correlated to exp(a2*airmass)
        # and is solved for during the fit
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    fig, axs = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()

    fig = myfit.plot_triangle()
    plt.tight_layout()
    plt.show()