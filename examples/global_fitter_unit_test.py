import numpy as np
import matplotlib.pyplot as plt

from exotic.api.elca import transit, glc_fitter

if __name__ == "__main__":

    # simulate input data
    epochs = np.random.choice(np.arange(100), 5, replace=False)
    input_data = []

    for i, epoch in enumerate(epochs):

        nobs = np.random.randint(50) + 100
        phase = np.linspace(-0.02-0.01*np.random.random(), 0.02+0.01*np.random.random(), nobs)
        
        prior = {
            'rprs':0.1, # Rp/Rs
            'ars':14.25,        # a/Rs
            'per':3.5,          # Period [day]
            'inc':87.5,         # Inclination [deg]
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

        #plt.plot(time,flux,marker='o')
        #plt.plot(time, model,ls='-')
        #plt.show()

    # shared properties between light curves
    global_bounds = {
        'per':[3.5-0.0001,3.5+0.0001],
        'tmid':[1-0.01,1+0.01],
        'ars':[14,14.5],
    }

    # individual properties
    local_bounds = {
        'rprs':[0,0.2],
        'a1':[0,1e4],
        'a2':[-0.5,0]
    }

    print('epochs:',epochs)
    myfit = glc_fitter(input_data, global_bounds, local_bounds, individual_fit=True, verbose=True)

    myfit.plot_triangle()
    plt.show()

    myfit.plot_bestfit()
    plt.show()