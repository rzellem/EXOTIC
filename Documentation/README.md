## Installation and Running

EXOTIC can run on a Windows, Macintosh, or Linux/Unix computer. You can also use EXOTIC via the free Google Colab, which features cloud computing, many helpful plotting functions, and a simplified installation. However, if you are a user with many images or large images, we recommend running EXOTIC locally on your own computer.

 **Google Colab Cloud**
  - Features: does not require the user to install any software locally on their own computer.
  - Limitations: Requires user to upload their images to a free Gdrive account.
  - Recommendations: If you run out of space on your default Google/Gdrive account, you can sign up for a new, free account to use. Some users even make a new Google account for every new dataset to avoid running out of space.
  - [How to use EXOTIC on the Colab video](https://drive.google.com/file/d/10zlQRgT8iV3dSe0FVW7tiL-V86ewai_1/view)
  - [How to use EXOTIC on the Colab written instructions](http://docs.google.com/document/d/1GLnfX1DdGPpd1ArKNcoF2GGV6pwKR3aEYuwjSQlhiZQ/edit?usp=sharing)
  - [EXOTIC: Google Colab Cloud Version](https://colab.research.google.com/drive/1UcDfm3z1WnfdOpRwjCQYwDgK9Wh2cU6x?usp=sharing) (includes step-by-step instructions)
  

 **Locally On Your Own Computer**
  - Features: Images are read off of the user's harddrive- nothing is uploaded to Gdrive. This method can be helpful for those with large filesizes, many files, or a slow internet connection.
  - Limitations: Requires user to install Python3 and multiple subpackages.
  
  - Installation Instructions:
    1. ​[Download and install the latest release of Python.](https://www.python.org/downloads/)
        **NOTE FOR WINDOWS USERS:** make sure to check the box "Add Python to PATH" when installing.
        **NOTE FOR ALL USERS:** please download and install the latest release of Python, even if you have a previous installation already on your computer, to ensure that all Python packages are properly installed.
    2. [Download the latest release of EXOTIC.](https://github.com/rzellem/EXOTIC/releases)
    3. Unzip this file.
    4. Double-click on the appropriate installer for your operating system:
        - Windows: run_exotic_windows.bat
        - Macintosh: run_exotic_macintosh.command
        - Linux: run_exotic_linux.sh
    5. If you get a security warning about the software being from an unidentified, unsigned, or non-trusted developer, you can bypass it by:
        - Windows: click "More info" and then the "Run away" box at the bottom of the window.
        - Macintosh: Please follow [these instructions](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac).

- **We also recommend that you download our [sample transiting exoplanet dataset](https://github.com/rzellem/EXOTIC_sampledata)** to confirm that EXOTIC is running correctly on the Google Colab Cloud or your own computer.
- How EXOTIC Works
  - [Document](https://github.com/rzellem/EXOTIC/blob/main/Documentation/English/How-EXOTIC-Works.pdf)
  - [Video](https://drive.google.com/file/d/1x0kl8WtpEw9wS0JInbjVWvdzuTc9TTvS/view)

- Lastly, we offer these documents [in other languages](https://github.com/rzellem/EXOTIC/raw/main/Documentation/)

## Requirements
FITS files with a modern header including parameters for UT time, exposure time, WCS coordinations (optional) are required for EXOTIC.

## Sample Data and Outputs
We provide a [sample dataset](https://github.com/rzellem/EXOTIC_sampledata) consisting of 142 `fits` files taken by a 6” telescope of the exoplanet HAT-P-32 b (V-mag = 11.44) observed on December 20, 2017. The telescope used to collect this dataset is part of the [MicroObservatory Robotic Telescope Network](http://microobservatory.org) operated by the Harvard-Smithsonian Center for Astrophysics.

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
