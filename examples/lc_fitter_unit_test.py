from exotic.exotic import LimbDarkening
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
        'ecc': 0.5,                                   # Eccentricity
        'omega': 120,                                 # Arg of periastron
        'tmid': 0.75,                               # Time of mid transit [day],
        'a1': 50,                                   # Airmass coefficients
        'a2': 0.,                                   # trend = a1 * np.exp(a2 * airmass)

        'teff':5000,
        'tefferr':50,
        'met': 0,
        'meterr': 0,
        'logg': 3.89, 
        'loggerr': 0.01
    }

    # example generating LD coefficients
    tessfilter = create_tess()

    ld_obj = LimbDarkening(
        teff=prior['teff'], teffpos=prior['tefferr'], teffneg=prior['tefferr'],
        met=prior['met'], metpos=prior['meterr'], metneg=prior['meterr'],
        logg=prior['logg'], loggpos=prior['loggerr'], loggneg=prior['loggerr'],
        wl_min=tessfilter.wl.min(), wl_max=tessfilter.wl.max(), filter_type="Clear")

    ld0, ld1, ld2, ld3, filt, wlmin, wlmax = ld_obj.nonlinear_ld()

    prior['u0'],prior['u1'],prior['u2'],prior['u3'] = [ld0[0], ld1[0], ld2[0], ld3[0]]

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