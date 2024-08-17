# EXOTIC Use Cases

This repository contains examples of how to use the EXOTIC software to perform a variety of tasks related to exoplanet transit science. The package is designed to be used with FITS images, photometric data, radial velocity data, and ephemeris data. The examples below are organized by the type of data used in the analysis.

## [Programmatic Access to Exoplanet Watch Results](Exoplanet_Watch_API.ipynb) 

Download parameters that are derived from photometric data and light curves for a given target. Over 400 targets and 6000 publicly available light curves from ground-based and space-based telescopes.

![](epw_results.png)

```python
import numpy as np
import matplotlib.pyplot as plt

from exotic.api.elca import transit
from exotic.api.ew import ExoplanetWatch

# This will load the results JSON from the link above
EW = ExoplanetWatch()
print(EW.target_list)

# names are case and space sensitive
result = EW.get('WASP-33 b')

# list the result properties
print(result.__dict__.keys())

# show the light curve
print(result.observations[0].lightcurve_url)

# extract the light curve data
time, flux, fluxerr, airmass, airmasscorr = result.observations[0].get_data()

# plot the data and best fit transit
plt.plot(time, flux/airmasscorr, 'ko')
plt.plot(time, transit(time, result.observations[0].parameters), 'r-')
plt.xlabel("Time [BJD]")
plt.ylabel("Rel. Flux")
```

Visit our [notebook](Exoplanet_Watch_API.ipynb) for more details on how to download the light curves from our results page programmatically.

## Acknowledgements

If you use any Exoplanet Watch data or in your publication, you are required to include the observers of those data as co-authors on your paper. To get in touch with your anonymous observer, contact the [AAVSO](https://www.aavso.org/) with their observer code.

If you make use of Exoplanet Watch in your work, please cite the papers [Zellem et al. 2020](https://ui.adsabs.harvard.edu/abs/2020PASP..132e4401Z/abstract) and [Pearson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022AJ....164..178P/abstract) and include the following standard acknowledgment in any published material that makes use of Exoplanet Watch data: **“This publication makes use of data products from Exoplanet Watch, a citizen science project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. This work is supported by NASA under award number NNX16AC65A to the Space Telescope Science Institute, in partnership with Caltech/IPAC, Center for Astrophysics|Harvard & Smithsonian, and NASA Jet Propulsion Laboratory."**

## [Fit a transit light curve](single_transit/transit_fit_example.py)

Start from a photometric timeseries and derive transit parameters like planetary radius, inclination and mid-transit time. Optimized in a Bayesian framework ([Ultranest](https://johannesbuchner.github.io/UltraNest/index.html)) with posteriors to assess degeneracies and uncertainties.

![](single_transit/bestfit.png)
![](single_transit/triangle.png)

```python
from exotic.api.elca import transit, lc_fitter
from pylightcurve import exotethys

# 
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

# generate limb darkening coefficients for TESS
get_prior = lambda key: float(prior[key])
u0,u1,u2,u3 = exotethys(get_prior('LOGG'), get_prior('T*'), get_prior('FE/H'), 'TESS', method='claret', stellar_model='phoenix')
prior['u0'],prior['u1'],prior['u2'],prior['u3'] = u0,u1,u2,u3

# create fake data if you don't have any
time = np.linspace(0.7, 0.8, 1000)  # [day]

# simulate extinction from airmass
stime = time-time[0]
alt = 90 * np.cos(4*stime-np.pi/6)
#airmass = 1./np.cos(np.deg2rad(90-alt)) # uncomment to simulate airmass
airmass = np.zeros(time.shape[0])

# GENERATE NOISY DATA
data = transit(time, prior)*prior['a1']*np.exp(prior['a2']*airmass)
data += np.random.normal(0, prior['a1']*250e-6, len(time))
dataerr = np.random.normal(300e-6, 50e-6, len(time)) + np.random.normal(300e-6, 50e-6, len(time))

# add optimization bounds for free parameters only
mybounds = {
    'rprs': [0, 0.1],
    'tmid': [prior['tmid']-0.01, prior['tmid']+0.01],
    'inc': [87,90],
    #'a2': [0, 0.3] # uncomment if you want to fit for airmass

    # a2 is used for individual airmass detrending using: a1*exp(airmass*a2)
    # a1 is solved for automatically using mean(data/model) and does not need
    # to be included as a free parameter. A monte carlo process is used after
    # fitting to derive uncertainties on it. It acts like a normalization factor.
    # never list 'a1' in bounds, it is perfectly correlated to exp(a2*airmass)
    # and is solved for during the fit
}

# call the fitting routine
myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

for k in myfit.bounds.keys():
    print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

# plot the best fit
fig, axs = myfit.plot_bestfit()
plt.tight_layout()
plt.savefig('bestfit.png')
plt.show()

# plot the posteriors
fig = myfit.plot_triangle()
plt.tight_layout()
plt.savefig('triangle.png')
plt.show()
```

## Fit multiple light curves simultaneously with shared and individual parameters

- [Simultaneous airmass detrending](multiple_transit/Multiple_Lightcurve_fit.ipynb) (more robust but takes much longer)

- [Airmass detrending prior to simultaneous fit](multiple_transit/Multiple_Lightcurve_Fit_Detrended.ipynb)

The notebooks above are also compatible with TESS data! Just don't include an `a2` parameter for airmass detrending in the local bounds.

![](multiple_transit/glc_fit.png)
![](multiple_transit/glc_mosaic.png)

![](tess/2458415_06_wasp-77ab_lightcurve.png)

```python
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
            'rprs':0.1, # Rp/Rs
            'ars':14.25, # a/Rs
            'per':3.5, # Period [day]
            'inc':np.random.random()+87.5, # Inclination [deg]
            'u0': 1.349, 'u1': -0.709, # exotethys - limb darkening (nonlinear)
            'u2': 0.362, 'u3': -0.087,
            'ecc':0, # Eccentricity
            'omega':0, # Arg of periastron
            'tmid':1, # time of mid transit [day],

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
```

## [Ephemeris fitting](ephemeris/fit_ephemeris.py)
- Observed - Calculated plot with colors coded to data source

![](ephemeris/oc.png)

- [Orbital Decay](ephemeris/fit_decay.py)

![](ephemeris/decay_oc.png)
![](ephemeris/decay_posterior.png)

- Periodogram for transit timing search with up to two orders

![](ephemeris/periodogram.png)

```python
import numpy as np
import statsmodels.api as sm
from exotic.api.ephemeris import ephemeris_fitter
import matplotlib.pyplot as plt

if __name__ == "__main__":
    Tc = np.array([  # measured mid-transit times
       2461656.170979  , 2460683.06352087, 2461312.22680483,
       2461840.72721957, 2461404.50457126, 2459614.88352437,
       2459967.2136158 , 2461250.70825625, 2460196.5097846 ,
       2459444.30884179, 2460297.17833986, 2460842.44956614,
       2461270.28536414, 2459447.10519731, 2459497.43952496,
       2460023.14013537, 2460445.37816892, 2461605.83844971,
       2459161.88372907, 2460878.80452072, 2460766.95256289,
       2460926.34200058, 2461712.09519907, 2460355.90036206,
       2459489.04986809, 2459379.99494715, 2459187.04825632,
       2460224.47242434, 2460979.4719767 , 2460406.23304265,
       2461700.91035602, 2460448.17416082, 2461611.43282057,
       2460132.19662872, 2460876.00842332, 2461446.45102159,
       2460395.04580418, 2460920.74932793, 2459463.87955256,
       2461756.83564536
    ])

    Tc_error = np.array([
       0.00086083, 0.00078861, 0.00086634, 0.00093534, 0.00075317,
       0.0008555 , 0.0007527 , 0.00078389, 0.00075229, 0.00076776,
       0.00042222, 0.00098135, 0.00075493, 0.00022053, 0.00038568,
       0.00070884, 0.00044426, 0.00043817, 0.00079945, 0.00084854,
       0.00053842, 0.00078082, 0.00079287, 0.00018737, 0.00018376,
       0.00062508, 0.00015116, 0.00034341, 0.00061096, 0.00083818,
       0.00076459, 0.00076027, 0.00081125, 0.00030851, 0.00015512,
       0.0007488 , 0.00027584, 0.00027871, 0.00080871, 0.00086118
    ])

    # labels for a legend
    labels = np.array([
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar'
    ])

    P = 2.7962868  # orbital period for your target

    Tc_norm = Tc - Tc.min()  # normalize the data to the first observation
    # print(Tc_norm)
    orbit = np.rint(Tc_norm / P)  # number of orbits since first observation (rounded to nearest integer)
    # print(orbit)

    # make a n x 2 matrix with 1's in the first column and values of orbit in the second
    A = np.vstack([np.ones(len(Tc)), orbit]).T

    # perform the weighted least squares regression
    res = sm.WLS(Tc, A, weights=1.0 / Tc_error ** 2).fit()
    # use sm.WLS for weighted LS, sm.OLS for ordinary LS, or sm.GLS for general LS

    params = res.params  # retrieve the slope and intercept of the fit from res
    std_dev = np.sqrt(np.diagonal(res.normalized_cov_params))

    slope = params[1]
    slope_std_dev = std_dev[1]
    intercept = params[0]
    intercept_std_dev = std_dev[0]

    # 3 sigma clip based on residuals
    calculated = orbit * slope + intercept
    residuals = (Tc - calculated) / Tc_error
    mask = np.abs(residuals) < 3
    Tc = Tc[mask]
    Tc_error = Tc_error[mask]
    labels = labels[mask]

    # print(res.summary())
    # print("Params =",params)
    # print("Error matrix =",res.normalized_cov_params)
    # print("Standard Deviations =",std_dev)

    print("Weighted Linear Least Squares Solution")
    print("T0 =", intercept, "+-", intercept_std_dev)
    print("P =", slope, "+-", slope_std_dev)

    # min and max values to search between for fitting
    bounds = {
        'P': [P - 0.1, P + 0.1],  # orbital period
        'T0': [intercept - 0.1, intercept + 0.1]  # mid-transit time
    }

    # used to plot red overlay in O-C figure
    prior = {
        'P': [slope, slope_std_dev],  # value from WLS (replace with literature value)
        'T0': [intercept, intercept_std_dev]  # value from WLS (replace with literature value)
    }

    lf = ephemeris_fitter(Tc, Tc_error, bounds, prior=prior, labels=labels)

    lf.plot_triangle()
    plt.subplots_adjust(top=0.9, hspace=0.2, wspace=0.2)
    plt.savefig("posterior.png")
    plt.close()

    fig, ax = lf.plot_oc()
    plt.tight_layout()
    plt.savefig("oc.png")
    plt.show()
    plt.close()

    fig, ax = lf.plot_periodogram()
    plt.tight_layout()
    plt.savefig("periodogram.png")
    plt.show()
    plt.close()

```

## [N-body interpretation of periodogram](nbody/README.md)

N-body simulations can be used to interpret the periodogram of transit timing variations (TTVs) and determine the masses of the planets in the system. The example below shows how to use the `nbody` module to generate a periodogram from an N-body simulation and compare it to the observed periodogram.

![](nbody/report.png)

```python
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
```

## [Radial Velocity](radial_velocity/rv_example.py)

![](radial_velocity/RV_bestfit.png)

```python
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

```

## [Joint Fit of Transit Photometry, Radial Velocity, and Ephemeris data (transit/eclipse times)](joint_rv_transit/joint_fit.py)

![](joint_rv_transit/joint_posterior.png)

## [Joint Fit of Transit Photometry, Radial Velocity, and Ephemeris data (transit/eclipse times)](joint_rv_transit/joint_fit.py)

![](joint_rv_transit/joint_posterior.png)

## [for_exotic_py_candidate_inits_maker.py](tess/candidates/for_exotic_py_candidate_inits_maker.py)

This script automates and enhances the process of generating initialization (inits) JSON files for candidate exoplanets specifically for use with exotic.py. It streamlines the preparation of inits files by interacting with the ExoFOP database and organizing the output data for easy integration with the EXOTIC pipeline.

Key Features:

Automated JSON File Download: Prompts the user to input a TIC ID, then automatically downloads the corresponding JSON file from ExoFOP.
Data Extraction and CSV Generation: Extracts relevant stellar and planetary parameters from the downloaded JSON files and compiles them into a comprehensive CSV file.
Inits JSON File Creation: Generates a JSON inits file tailored for use with exotic.py, with options to estimate missing parameters if necessary.
Directory Management: Organizes the storage of output files by dynamically creating directories for each TIC ID and candidate designation, ensuring consistent and clean file management.
File Location:
/EXOTIC/examples/tess/candidates/
The generated output files will be located in:
/EXOTIC/examples/tess/candidates/output_inits_files/for_exotic_py_candidate_inits_output_files/{tic_id}_file/{candidate_designation}_inits.json


## [for_toi_py_candidate_inits_maker.py](tess/candidates/for_toi_py_candidate_inits_maker.py)

This script enhances the generation of initialization (inits) JSON files specifically for use with toi.py and toi_indiv_lc.py. It automates the process by interacting with ExoFOP, organizing the output data, and streamlining the preparation of inits files.

Key Features:

Automated JSON File Download: Prompts the user to input a TIC ID and downloads the corresponding JSON file from ExoFOP.
Data Extraction and CSV Generation: Extracts and compiles relevant stellar and planetary parameters into a comprehensive CSV file.
Directory Management: Dynamically organizes output files for each TIC ID and candidate designation, ensuring consistent and clean file management.
File Location:

The generated output files will be located in:
/EXOTIC/examples/tess/candidates/output_inits_files/for_toi_py_candidate_inits_output_files/{tic_id}_file/{candidate_designation}_inits.json
TESS Light Curve Analysis


## [toi.py](tess/candidates/toi.py)
This script is used to process and analyze TESS light curve data for candidate exoplanets. It automates the process of downloading, fitting, and analyzing transit data. It is a modification from the script tess.py.

Key Features:

Automated Data Handling: Downloads and processes light curve data from TESS based on the provided initialization file (can be created with "for_toi_py_candidate_inits_maker.py").
Transit Fitting: Performs transit fitting and analysis using the provided parameters, generating visual and data outputs.
Output Management: Saves the results, including processed light curves and fitting results, in the specified output directory.
Usage:
To run the script:

Navigate to the directory containing the script:
cd /EXOTIC/examples/tess/candidates
Run the script with the following command:
python toi.py -i <input_inits_file_path> -o <output_directory_path>




## [toi_indiv_lc.py](tess/candidates/toi_indiv_lc.py)
This script is used to process and analyze TESS light curve data for candidate exoplanets. Unlike toi.py, which performs a global fit, toi_indiv_lc.py focuses on fitting individual light curves using the priors from an existing initialization file.

Key Features:

Automated Data Handling: Downloads and processes light curve data from TESS based on the provided initialization file (can be created with for_toi_py_candidate_inits_maker.py).
Transit Fitting: Uses provided priors to perform transit fitting and analysis on individual light curves, generating visual and data outputs without performing a global fit.
Output Management: Saves the results, including processed light curves and fitting results, in the specified output directory.
Usage:
To run the script:

Navigate to the directory containing the script:
cd /EXOTIC/examples/tess/candidates
Run the script with the following command:
python toi_indiv_lc.py -i <input_inits_file_path> -o <output_directory_path>
Example:

python toi_indiv_lc.py -i input_inits.json -o output_results/



## [TESS light curve generation](tess/tess_individ_lc.py)

The `tess_individ_lc.py` script allows for the extraction and fitting of individual exoplanet transits from TESS data without performing a global fit. This script is beneficial for users needing quicker access to individual light curve data for joint fit or individual light curve analysis with other Python tools like EXOTIC example notebooks.

- **Works with confirmed exoplanets.**
- **Created from `tess.py`. In `tess.py`, the global fit is performed first, and the derived parameters are used as priors in the individual fits. In `tess_individ_lc.py`, the NEA priors are used instead, potentially giving slightly different results.**
- **Outputs go to `tess/tess_individ_lc_output` folder.**