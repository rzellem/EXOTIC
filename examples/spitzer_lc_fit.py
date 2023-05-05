from exotic.api.elca import transit, lc_fitter
from pylightcurve import exotethys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    prior = {
        'ecc': 0.0062,
        'inc': 82.33,
        'omega': 90.0,
        'duration': 0.0504,
        'tmid': 2456401.7295,
        'a1': 1.138,
        'a2': -0.086,
        'ars': 4.93,
        'rprs': 0.155,
        'per': 0.8134749,

        'TEFF': 4400,
        'LOGG': 4.65,
        'FE/H': -0.05,
    }

    # generate limb darkening coefficients for Spitzer IRAC1 (3.6 micron)
    get_prior = lambda key: float(prior[key])
    u0,u1,u2,u3 = exotethys(get_prior('LOGG'), get_prior('TEFF'), get_prior('FE/H'), 'irac1', method='claret', stellar_model='atlas')
    prior['u0'],prior['u1'],prior['u2'],prior['u3'] = u0,u1,u2,u3

    # load data from file
    data = np.loadtxt('wasp43b_irac1_transit.txt')

    time = data[:,0]
    flux = data[:,1]
    fluxerr = data[:,2]

    # wx, wy, npp, ramp
    syspars = data[:,3:]

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, prior['rprs']*2],
        'tmid': [prior['tmid']-0.01, prior['tmid']+0.01],
        'inc': [prior['inc']-5, min(90,prior['inc']+5)],
    }

    # run the fit
    myfit = lc_fitter(time, flux, fluxerr, syspars, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    fig, axs = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()

    if myfit.mode == 'ns':
        fig = myfit.plot_triangle()
        plt.tight_layout()
        plt.show()