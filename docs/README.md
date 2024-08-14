- [Installation and Running](#installation-and-running)
- [Requirements](#requirements)
- [Sample Data and Outputs](#sample-data-and-outputs)
- [How to Run EXOTIC with the Sample Data Locally on your Computer](#how-to-run-exotic-with-the-sample-data-locally-on-your-computer)
- [If you run into any issues with EXOTIC](#if-you-run-into-any-issues-with-exotic)
- [Initializaton File](#initializaton-file)


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

[Sample Data](https://github.com/rzellem/EXOTIC_sampledata/releases/)

A lightcurve from the sample dataset is shown below:

![Lightcurve graph showing relative flux versus phase with error bars and interpolated curve.](https://github.com/rzellem/EXOTIC/raw/main/docs/images/HAT-P-32bExample.png)

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

## How to Run EXOTIC with the Sample Data Locally on your Computer
1. Download the [Sample Data](https://github.com/rzellem/EXOTIC_sampledata/releases/) and unzip the file.

2. Launch EXOTIC by double-clicking on the appropriate installer for your operating system:
      - Windows: run_exotic_windows.bat
      - Macintosh: run_exotic_macintosh.command
      - Linux: run_exotic_linux.sh

3. A new loading screen will pop up where you will see many required packages being automatically installed and updated. If this is your first time running the code, this can take up to a few minutes, depending on your internet connection.

![EXOTIC Loading Screen](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_loading.png)

4. After EXOTIC finishes loading, you will see a welcome screen. Click Next to proceed.

  ![EXOTIC Welcome Screen](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_welcome.png)

5. Select how you want to use EXOTIC - either in "Real Time Reduction" (for quickly analyzing your data while simulatenously observing) or in "Complete Reduction" (for analyzing your data after an observing run). Since we already have our sample data on our computer harddrive, we will select "Complete Reduction" and click Next.

![EXOTIC Reduction Mode](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_redmode.png)

6. Select how you want to run EXOTIC: either by starting with raw image files (typically in FITS format or similar) or by starting with pre-reduced data (i.e., data already reduced with other software). Since we want to fully analyze our sample data, which are image files, we will select "Analyze your image files" and click Next.

![EXOTIC Ru Mode](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_runmode.png)

7. Select how you would like to enter information about your observing run (e.g., location, sky conditions, position of the target and comparison star(s) on the detector, etc.). For our sample data run, select "Manually".

![EXOTIC Input Observation Information](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_inputobsinfo.png)

8. Now enter the information about your observing run.
    - Folder with your FITS files - the folder that contains your image files from your observation run; typically, this is a file called "Downloads"
    - Folder to save EXOTIC output - select where you want EXOTIC to write your output, such as final lightcurves and the file you can upload directly to the AAVSO Exoplanet Database
    - Please enter your calibration file information - enter the location of your flats, darks, and biases
        - NOTE: these three calibrations should be in DIFFERENT/SEPARATE folders and should *not* be in the same folder as your observation data
        - If you do not have any of these calibrations, enter `null`
    - AAVSO Observer Code - if you do not have one, leave as N/A
    - Secondary Observer Codes - the AAVSO observer codes of anyone who helped out with your observations; if you do not have one, leave as N/A
    - Observation date - the date of your observation in DAY-MONTH-YEAR format
    - Obs. Latitude - the latitude of your observations, where North is denoted with a + and South is denoted with a -
    - Obs. Longitude - the longitude of your observations, where East is denoted with a + and West is denoted with a -
    - Obs. Elevation - the elevation of your observations in meters; note: if you do not know your elevation, leave this value as 0 and EXOTIC will look it up for you!
    - Camera Type - the type of camera you are using (typically, CCD, DSLR, or CMOS)
    - Pixel Binning - the binning of your pixels; if you do not bin your pixels, then please enter 1x1
    - Pixel Scale - how large each pixel is in units of arcseconds per pixel; if you do not know your pixel scale, please enter null and EXOTIC will calculate it for you
    - Filter - select the filter that you used for your observations; if you used a custom filter, please select "N/A" and EXOTIC will ask you for your filter specifics later; if you did not use a filter, please select either select "N/A" and then you will be prompted for wavelength information later or select the filter that most closely corresponds to your camera
    - Observing Notes - enter notes about your observing run, such as weather conditions
    - Plate solve my images - select if you want EXOTIC to calibrate the right ascenscion and declination of your pixels in your image via Astrometry.net; it is recommended that this option is selected
    - Align my images - select this option for EXOTIC to align all of your images to provide better tracking of your stars in your images; it is recommended that this option is selected
    - Target Star X & Y Pixel Position - the pixel location of your target exoplanet host star in [x-position, y-position] format
    - Comparison Star(s) X & Y Pixel Position - the pixel location of your comparision star(s) in [x-position, y-position] format; it is recommended that you input at least 2 comparision stars and EXOTIC will automatically select the "best" comparision by the one that produces the least amount of scatter in your data
    - *NOTE:* In the screenshot below, Rob has already entered all of the information for you for the sample data (with the exception that you'll need to point to the correct directory for your FITS files and your EXOTIC Output)
    
    ![EXOTIC Input Observation Information](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_inputobs.png)
    
9. Select how you want to enter your planetary system parameters. For first-time users, it is recommended that you enter these "Manually". Advanced users might want to simply adopt all planetary parameters from the NASA Exoplanet Archive. If you have already input your parameters in a pre-existing initialization file, then please select the last option.

![EXOTIC Input Planetary Parameters](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_inputplanetparams.png)

10. Enter the planetary parameters for your exoplanetary system. You can find this information on a variety of resources, but we recommend you look up your planet on the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu). Note to users who are reducing the sample HAT-P-32b dataset- Rob has already filled out the system parameters for HAT-P-32b.

![EXOTIC Input Planetary Parameters](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_planetparams.png)

11. EXOTIC will now create an initialization file for you. This file essentially saves all of the answers to the questions and prompts you just filled out- that way, you can just bypass all of these questions the next time you reduce these observations by simply loading in this initialization file. Please select the location on your harddrive where you want to save this file. It will be then named as `inits_MM_DD_YYYY__HH_MM_SS.json` where MM-DD-YYYY is the month, day, and year and HH_MM_SS is the current time.

![EXOTIC Save Initialization File](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_saveinits.png)
    
12. Congratulations! EXOTIC has successfully saved your initialization file! Click "Run EXOTIC" to launch EXOTIC and start analyzing the sample data!

![EXOTIC Initialization File Saved](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_initssaved.png)

13. EXOTIC is now running!

![EXOTIC is now Running! Saved](https://github.com/rzellem/EXOTIC/blob/develop/docs/images/exotic_running.png)

## If you run into any issues with EXOTIC
- Please send a message on our Slack Workspace on the #data-reduction channel to get help from the Exoplanet Watch Community. You can sign up for your own free Slack account by clicking on [this link](https://join.slack.com/t/uol-ets/shared_invite/zt-mvb4ljbo-LRBgpk3uMmUokbs4ge2JlA).

- Alternatively, you can email the Exoplanet Watch team directly at exoplanetwatch@jpl.nasa.gov.


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
