import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

from exotic.api.nbody import report, generate, nbody_fitter, analyze, estimate_prior, TTV, interp_distribution

mearth = u.M_earth.to(u.kg)
msun = u.M_sun.to(u.kg)
mjup = u.M_jup.to(u.kg)

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
