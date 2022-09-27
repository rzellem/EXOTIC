# EXOTIC Installation Instructions for Windows 10

## Running EXOTIC
There are two ways you can run EXOTIC: firstly through a web browser interface with a so-called Python Notebook running in Google Colab or in a Jupyter Notebook running locally on your own personal computer. 

Alternatively, you can install both Python and the EXOTIC program to your own computer and run them there. The program uses a combination of a graphical front-end to enter values and a command line window where you can view progress as the program runs.

Read the details of these two approaches below and select the best option for you:

### Google Colab

Google Colab is a free facility (available to anyone with a Google account) where you may write and run Python programs within your web browser without installing anything on your local machine. Running EXOTIC in the Google Colab is the better option for all users who are seeking a more user-friendly and interactive experience. It is especially useful for new EXOTIC users, students and those analysing data via remote observatories. 

To use Google Colab to run EXOTIC follow the links on the home page of the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC).

Once you are happy running running EXOTIC with the Colab and you want to run your own Python Notebook, please visit the instructions [‘How-to-Run-EXOTIC-with-the-Python-Notebook](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/How-to-Run-EXOTIC-with-the-Python-Notebook.pdf)’ in the [Documentation folder](https://github.com/rzellem/EXOTIC/tree/main/Documentation) on GitHub.

*Note: in order to be able to click on the links and select text in this PDF document, you must download it from GitHub. The GitHub preview simply shows you an image of the document, which does not allow for those functions.*

### GUI Version of EXOTIC for Windows

Early versions of EXOTIC had to be run wholly on the "command line" whereby responses to the prompts from the program had to be typed in line by line as the program ran each step. EXOTIC itself still runs in a command line window but you can now enter your choices for the various program options via a series of pop up screens before commencing the data-reduction run.

## Installing EXOTIC for Windows

1. #### Install Python for Windows 10 

If you do not already have a copy of the latest version of Python software on your computer, then you will need to install it before you can run EXOTIC. Even if you have an older copy of Python then do install another copy to get the latest release, which is needed to run EXOTIC.

Go to the [Python Releases for Windows](https://www.python.org/downloads/windows/) page on the Python website, where you will find the download files for Windows. Look at the list of "Stable Releases" and take the link to the most recent release of Python towards the top of that list. At the time of writing this was Python 3.9.5 - but it may be a later release for you. 

If you have a modern computer (i.e. one running the 64-bit version of Windows) then choose the "Windows Installer (64-bit )" link, otherwise if you are using an older machine you can use the "Windows Installer (32-bit )" version.

Right click to download and save the EXE file to your computer and double click it to start the installation process. 

- Choose the "Modify" option on the "Modify Setup" screen that pops up first. 
- On the next "Optional Features" screen, make sure that the boxes for "pip" and for "py launcher" are ticked then click "Next". 
- On the "Advanced Options" screen make sure that the boxes "Associate files with Python", "Create shortcuts for installed applications" and "Add Python to environment variables" are all ticked before clicking "Install". 
- Do not change the default install location or you may run into problems later on when using EXOTIC.

*Note: You can install Python direct from the Microsoft store, where you will find a recent copy of Python but we advise against this approach. Although this will install the program like any other windows software, it will not be installed in a similar way to the Python installed on other systems that do not run Windows. EXOTIC runs on many platforms, including Linux as well as Windows, and expects to find a similar set up on all of them to run without problems.*

2. #### Downloading EXOTIC and supporting files

2.1 Go to the home page of the [EXOTIC project](https://pypi.org/project/exotic/) at the Python Package Index (PyPI). This site hosts many major Python projects for a variety of computer platforms, including both Linux and Windows machines. Here you will find links to many components of the EXOTIC project, including copies of the documentation and the EXOTIC software itself. 

2.2 Click the "Download files" link in the left hand margin of the EXOTIC home page on PyPi. On the next page you will find files in various formats for the latest release of EXOTIC, each with a different file extension. The archive file you need will have a name beginning "exotic-" followed by the release number and finally the extension ".zip" (for example "exotic-1.3.3.zip)". Right click the file name link to download the compressed archive containing the full set of files. 

2.3 Go to the folder where you downloaded the file and right click the zip file that you have just saved and choose the "Extract All..." option. In the window that pops up, choose where to install the EXOTIC files, make sure the "Show extracted files when complete" box is ticked and click "Extract". After extraction, a new window will open up showing the folder where the files are located, named after the version of the EXOTIC program (e.g. EXOTIC1.3.3). 

3. #### Installing EXOTIC

In the file manager window that has just appeared, open the EXOTIC program folder to find the batch file named "run_exotic_windows.bat". You may find that there is just another folder named after the EXOTIC version inside the first folder, in which case open this as well. Having found it, double click the batch file to start the installation process.

You may see a warning window popping up from Microsoft Defender reading: 

###### *Windows Protected your PC* 

*Microsoft Defender SmartScreen prevented an unrecognised app from starting. Running this app might put your PC at risk.*

*Application: run_exotic_windows.bat* 
*Publisher:  Unknown publisher* 

If so you need not worry, this is just a precaution from Microsoft when it does not know the details of the origin of software, so on this occasion you can safely ignore the message and click "Run anyway"

A Windows command line window with a black background will then open up, showing you the progress of the installation for each of the files required to run EXOTIC. This process may take some time and you will see many lines of yellow and white text scroll by until the whole installation is finished. Do not worry when the process pauses on a line starting with "Building wheel for ...", since each of these steps will take some time. 

The final stage, beginning "Installing collected packages:..." followed by a long list of packages, will take several minutes to run and the screen will not change at all during this time - please be patient and just wait for the installation to complete on its own.

When the script completes, a new command line window in black appears that will show the progress of the data-reduction once EXOTIC is running. In front of this you will see the EXOTIC GUI window showing "Welcome to EXOTIC!" plus the version of EXOTIC you are running and other details. This GUI window is used to enter the values required to define how the data-reduction run will be performed.

You can move each window around separately to see what is happening in each of them.

## Running EXOTIC next time

To run the EXOTIC program in the future, simply click locate and run the "run_exotic_windows.bat" file again. This will firstly check whether any EXOTIC files need updating and update them automatically if there has been another release of the software. In this case you will see the same installation process as outlined above, though it will run very much faster, as most files will not have changed between releases. 

After the installation of any new or updated EXOTIC components has completed (or is skipped through if nothing has changed) the program will run as before...

## Install DS9 (Astronomical Image Viewing Software)

After installing EXOTIC itself, you may wish to install the SAOImageDS9 (Astronomical Image Viewing Software) to view the “.FITS” images you obtain during observations. For more information on SAOImageDS9, check out the [User Guide](http://ds9.si.edu/doc/user/index.html)


- Follow [SAOImageDS9 download link](https://sites.google.com/cfa.harvard.edu/saoimageds9/download)
- Click on "***Windows***" then select "***Windows 64 bit***".
- Run the Windows installer in the usual way once the file has been downloaded.
- Follow the instructions in the installer to complete the installation.

## Further Help and Support

If you have any questions or comments, please feel free to reach out to us on Slack or email at [exoplanetwatch@jpl.nasa.gov](exoplanetwatch@jpl.nasa.gov)



