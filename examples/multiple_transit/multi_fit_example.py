import numpy as np
import matplotlib.pyplot as plt

from exotic.api.elca import transit, glc_fitter

if __name__ == "__main__":

    # simulate input data
    epochs = np.random.choice(np.arange(100), 8, replace=False)
    input_data = []
    local_bounds = []

    for i, epoch in enumerate(epochs):

        nobs = np.random.randint(50) + 100
        phase = np.linspace(-0.02-0.01*np.random.random(), 0.02+0.01*np.random.random(), nobs)
        
        prior = {
            'rprs':0.1,         # Rp/Rs
            'ars':14.25,        # a/Rs
            'per':3.5,          # Period [day]
            'inc':np.random.random()+87.5,         # Inclination [deg]
            'u0': 1.349, 'u1': -0.709,  # exotethys - limb darkening (nonlinear)
            'u2': 0.362, 'u3': -0.087,
            'ecc':0,            # Eccentricity
            'omega':0,          # Arg of periastron
            'tmid':1,        # time of mid transit [day],

            'a1':5000 + 2500*np.random.random(),   # airmass coeffcients
            'a2':-0.25 + 0.1*np.random.random()
        }

        time = prior['tmid'] + prior['per']*(phase+epoch)
        stime = time-time[0]
        alt = 90* np.cos(4*stime-np.pi/6)
        airmass = 1./np.cos( np.deg2rad(90-alt))
        model = transit(time, prior)*prior['a1']*np.exp(prior['a2']*airmass)
        flux = model*np.random.normal(1, np.mean(np.sqrt(model)/model)*0.25, model.shape)
        ferr = flux**0.5

        input_data.append({
            'time':time,
            'flux':flux,
            'ferr':ferr,
            'airmass':airmass,
            'priors':prior
        })

        # individual properties
        local_bounds.append({
            #'rprs':[0,0.2], # will overwrite global bounds if included
            # a2 is used for individual airmass detrending using: a1*exp(airmass*a2)
            'a2':[-0.75,0.25] 
            # a1 is solved for automatically using mean(data/model) and does not need
            # to be included as a free parameter. A monte carlo process is used after
            # fitting to derive uncertainties on it. It acts like a normalization factor.
        })

        #plt.plot(time,flux,marker='o')
        #plt.plot(time, model,ls='-')
        #plt.show()

    # shared properties between light curves
    global_bounds = {
        'rprs':[0,0.2],

        'per':[3.5-0.001,3.5+0.001],
        'tmid':[1-0.01,1+0.01],
        'inc':[87,90],
    }

    print('epochs:',epochs)
    myfit = glc_fitter(input_data, global_bounds, local_bounds, individual_fit=False, verbose=True)

    myfit.plot_bestfit()
    plt.tight_layout()
    plt.savefig('glc_fit.png')
    plt.close()
    #plt.show()

    myfit.plot_triangle()
    plt.tight_layout()
    plt.savefig('glc_triangle.png')
    plt.close()
    #plt.show()

    myfit.plot_bestfits()
    plt.tight_layout()
    plt.savefig('glc_mosaic.png')
    plt.show()
    plt.close()