# EXOTIC (EXOplanet Transit Interpretation Code)

A python 3 package for reducing photometric data of transiting exoplanets into lightcurves, and retrieving transit epochs and planetary radii.

The EXOplanet Transit Interpretation Code relies upon the transit method for exoplanet detection. This method detects exoplanets by measuring the dimming of a star as an orbiting planet transits, which is when it passes between its host star and the Earth. If we record the host star’s emitted light, known as the flux, and observe how it changes as a function of time, we should observe a small dip in the brightness when a transit event occurs. A graph of host star flux vs. time is known as a lightcurve, and it holds the key to determining how large the planet is, and how long it will be until it transits again.

![Light Curve Graph displaying brightness versus time.](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Images/transitsimple.jpg)

The objective of this pipeline is to help you reduce your images of your transiting exoplanet into a lightcurve, and fit a model to your data to extract planetary information that is crucial to increasing the efficiency of larger observational platforms, and futhering our astronomical knowledge.

## Installation (instalação)

- [English](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English)

- [Português do Brasil](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Brazilian_Portuguese)


The easiest way to install exotic is with pip: 

`$ pip install exotic`

**Depending on your version of python you may need to use a different pip command (e.g. pip3).** If you're having trouble installing exotic from pip, please see our documentation for additional installation instructions including setting up dependencies for [Mac](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/EXOTIC-Installation-Instructions-for-Mac-Users.pdf), [Windows](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/EXOTIC-Installation-Instructions-for-Windows-Users.pdf) and [Linux](exotic_installation_linux.sh)

## Examples 
- [Google Collab](https://colab.research.google.com/drive/1W1vrgEp9_IjEN16WFxmopLDYvHeYfxpw?usp=sharing)
- [Jupyter Notebook](examples/Exotic_Notebook.ipynb)
- Command Line: `$ exotic`

or if you have an inits file already:

`$ exotic -i inits.json`

FITS files with a modern header including parameters for UT time, exposure time, WCS coordinations (optional) are required for EXOTIC. We provide a sample dataset consisting of 142 `fits` files taken by a 6” telescope of the exoplanet HAT-P-32 b (VMag = 11.44) observed on December 20, 2017. The telescope used to collect this dataset is part of the MicroObservatory Robotic Telescope Network operated by the Harvard-Smithsonian Center for Astrophysics.

[Sample Data](https://github.com/rzellem/EXOTIC_sampledata)

A lightcurve from the sample dataset is shown below:

![Lightcurve graph showing relative flux versus phase with error bars and interpolated curve.](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Images/HAT-P-32bExample.png)

For the full output of EXOTIC please see the [example output](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/example_output.txt)

```
*********************************************************
FINAL PLANETARY PARAMETERS

              Mid-Transit Time [BJD]: 2458107.714007 +- 0.000856 
  Radius Ratio (Planet/Star) [Rp/Rs]: 0.1503 +- 0.0009 
 Semi Major Axis/ Star Radius [a/Rs]: 5.146 +- 0.059 
               Airmass coefficient 1: 7397.280 +- 19.7116 
               Airmass coefficient 2: -0.1161 +- 0.0021 
The scatter in the residuals of the lightcurve fit is: 0.5414 %

*********************************************************
```



## Initializaton File

Get EXOTIC up and running faster with a json file. Please see the included file ([inits.json](inits.json)) meant for the [sample data](https://github.com/rzellem/EXOTIC_sampledata). The initialization file has the following fields: 

```json
{
    "user_info": {
            "Directory with FITS files": "sample-data/HatP32Dec202017",
            "Directory to Save Plots": "sample-data/",
            "Directory of Flats": null,
            "Directory of Darks": null,
            "Directory of Biases": null,

            "AAVSO Observer Code (N/A if none)": "RTZ",
            "Secondary Observer Codes (N/A if none)": "N/A",

            "Observation date": "December 17, 2017",
            "Obs. Latitude": "+31.68",
            "Obs. Longitude": "-110.88",
            "Obs. Elevation (meters)": 2616,
            "Camera Type (CCD or DSLR)": "CCD",
            "Pixel Binning": "1x1",
            "Filter Name (aavso.org/filters)": "V",
            "Observing Notes": "Weather, seeing was nice.",

            "Plate Solution? (y/n)": "n",

            "Target Star X & Y Pixel": [424, 286],
            "Comparison Star(s) X & Y Pixel": [[465, 183], [512, 263]]
    },
    "planetary_parameters": {
            "Target Star RA": "02:04:10",
            "Target Star Dec": "+46:41:23",
            "Planet Name": "HAT-P-32 b",
            "Host Star Name": "HAT-P-32",
            "Orbital Period (days)": 2.1500082,
            "Orbital Period Uncertainty": 1.3e-07,
            "Published Mid-Transit Time (BJD-UTC)": 2455867.402743,
            "Mid-Transit Time Uncertainty": 4.9e-05,
            "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.14856104152345367,
            "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.004688608636917226,
            "Ratio of Distance to Stellar Radius (a/Rs)": 5.344,
            "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.04,
            "Orbital Inclination (deg)": 88.98,
            "Orbital Inclination (deg) Uncertainity": 0.68,
            "Orbital Eccentricity (0 if null)": 0.159,
            "Star Effective Temperature (K)": 6001.0,
            "Star Effective Temperature (+) Uncertainty": 88.0,
            "Star Effective Temperature (-) Uncertainty": -88.0,
            "Star Metallicity ([FE/H])": -0.16,
            "Star Metallicity (+) Uncertainty": 0.08,
            "Star Metallicity (-) Uncertainty": -0.08,
            "Star Surface Gravity (log(g))": 4.22,
            "Star Surface Gravity (+) Uncertainty": 0.04,
            "Star Surface Gravity (-) Uncertainty": -0.04
    },
    "optional_info": {
            "Pixel Scale (Ex: 5.21 arcsecs/pixel)": null,
            "Filter Minimum Wavelength (nm)": null,
            "Filter Maximum Wavelength (nm)": null
    }
}
```




## Features/ Pipeline Architecture

- Automatic Plate Solution from http://nova.astrometry.net

- Resolve targets with [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/) + retrieve light curve priors

- Hot Pixel Masking

![](Documentation/Images/Hot_pixel_mask.png)

- Aperture Photometry with PSF centroiding (2D Gaussian + rotation)

![HAT-P-32 b Centroid Position Graph, X-Pixel versus Time in Julian Date.](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Images/centroids.png)

- Stellar masking in background estimate

![](Documentation/Images/Background_Estimate.png)

- Multiple comparison star + aperture size optimization

- Non-linear 4 parameter limb darkening with [LDTK](https://github.com/hpparvi/ldtk)

- Light curve parameter optimization with [Nested Sampling](https://dynesty.readthedocs.io/en/latest/index.html)

![Chart showing how Nested Sampling iterations reveal light curve optimization results.](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Images/posterior_sample.png)



## Contributing to EXOTIC

EXOTIC is an open source project that welcomes contributions. Please fork the repository and submit a pull request to the develop branch for your addition(s) to be reviewed. 

## Citation
If you use any of these algorithms in your work, please cite our 2020 paper: [Zellem, Pearson, Blaser, et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200309046Z/abstract) 

![https://exoplanets.nasa.gov/exoplanet-watch/about-exoplanet-watch/](https://github.com/rzellem/EXOTIC/raw/main/Documentation/Images/ExoplanetWatch.png)

Contribute to [Exoplanet Watch](https://exoplanets.nasa.gov/exoplanet-watch/about-exoplanet-watch/), a citizen science project that improves the properties of exoplanets and their orbits using observations processed with EXOTIC. Register with [AAVSO](https://www.aavso.org/exoplanet-section) and input your Observer Code to help track your contributions allowing for proper credit on future publications using those measurements. Ask about our Exoplanet Watch Slack Channel!

## Acknowledgements 
Exoplanet Watch is a project by NASA's Universe of Learning. NASA's Universe of Learning materials are based upon work supported by NASA under award number NNX16AC65A to the Space Telescope Science Institute, working in partnership with Caltech/IPAC, Center for Astrophysics | Harvard & Smithsonian, Jet Propulsion Laboratory, and Sonoma State University.
