# EXOTIC Installation Instructions for Ubuntu

## Running EXOTIC
There are two ways you can run EXOTIC: firstly through a web browser interface with a so-called Python Notebook running in Google Colab or in a Jupyter Notebook running locally on your own personal computer. 

Alternatively, you can install both Python and the EXOTIC program to your own computer and run them there. The program uses a combination of a graphical front-end to enter values and a command line window where you can view progress as the program runs.

Read the details of these two approaches below and select the best option for you:

### Google Colab

Google Colab is a free facility (available to anyone with a Google account) where you may write and run Python programs within your web browser without installing anything on your local machine. Running EXOTIC in the Google Colab is the better option for all users who are seeking a more user-friendly and interactive experience. It is especially useful for new EXOTIC users, students and those analysing data via remote observatories. 

To use Google Colab to run EXOTIC follow the links on the home page of the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC).

Once you are happy running running EXOTIC with the Colab and you want to run your own Python Notebook, please visit the instructions [‘How-to-Run-EXOTIC-with-the-Python-Notebook](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/How-to-Run-EXOTIC-with-the-Python-Notebook.pdf)’ in the [Documentation folder](https://github.com/rzellem/EXOTIC/tree/main/Documentation) on GitHub.

*Note: in order to be able to click on the links and select text in this PDF document, you must download it from GitHub. The GitHub preview simply shows you an image of the document, which does not allow for those functions.*

### GUI Version of EXOTIC for Ubuntu

Early versions of EXOTIC had to be run wholly on the "command line" whereby responses to the prompts from the program had to be typed in line by line as the program ran each step. EXOTIC itself still runs in a command line window but you can now enter your choices for the various program options via a series of pop up screens before commencing the data-reduction run.

## Installing EXOTIC for Ubuntu

1. #### Installing/updating Python and supporting files

   Execute the commands shown in the boxes below by typing them on the Ubuntu command line after the prompt symbol "$". After entering each command you will need to press the Enter (or Return) key and wait for it to complete, which is when the prompt reappears.

   Start by updating your Ubuntu system to the latest version and then proceed to install the additional software necessary to run the EXOTIC program.

   1.1 Update the list of software packages available in the repository of Ubuntu software. The "sudo" command enables you to run packages normally requiring root privileges, so you will need to enter the root password when prompted. This command will show you all the packages that can be upgraded and any available essential new packages and then download them for you.


   ```python
     $ sudo apt update
   ```

   1.2 Update each package  to its latest version. You will see a list of the packages being upgraded or installed, as well as a progress bar showing the progress of the overall upgrade to your Ubuntu system.


   ```python
      $ sudo apt upgrade -y
   ```

   1.3  Next install pip, which is a tool needed to install the rest of the EXOTIC software.

   ```python
   $ sudo apt install python3-pip
   ```

   1.4 When the installation is complete, verify the installation by checking the python and pip versions - enter:

   ```python
   $ pip3 --version
   ```

   The version number may vary, but it will look something like this:

   ```output
   $ pip 21.1.2 from /usr/lib/python3/dist-packages/pip (python 3.9)
   ```

   

2. #### Downloading EXOTIC and supporting files

   2.1 Go to the home page of the [EXOTIC project](https://pypi.org/project/exotic/) at the Python Package Index (PyPI). This site hosts many major Python projects for a variety of computer platforms, including both Linux and Windows machines. Here you will find links to all aspects of the EXOTIC project including the documentation and the software itself. Follow the "Homepage" link in the left hand margin that will take you to the [EXOTIC home page on GitHub](https://github.com/rzellem/EXOTIC), where the software for the project is developed. 

   2.2 Click the "Releases" link in the right hand margin of the EXOTIC home page on GitHub. On the next page you will find a series of releases of EXOTIC, each with a header outlining the key change for that version. The most recent version of EXOTIC (labelled "Latest release" in the left hand margin) is shown at the top. Go to the "Assets" section to get a copy of the latest release and right click the "Source code (tar.gz)" link to download the compressed archive containing the full set of files. This archive file will be named "EXOTIC" followed by the release number e.g. "EXOTIC-0.45.3.zip". Save the file to the usual folder on your computer where your browser stores downloaded  files. 

3. #### Installing EXOTIC

   3.1 Create the folder you want to locate the EXOTIC files in, replacing **<exotic>** with the directory of your choice. *Note the tilde symbol "~" is just a shortcut reference to the location of your home directory.


   ```python
   $ mkdir ~/<exotic>
   ```

   3.2 Switch to the downloads directory and get ready to copy the software – replacing **<downloads-directory>** with the full location of the directory in which your downloaded EXOTIC files are located. If you are using the Firefox browser in Ubuntu (and have not changed the default settings) then your downloaded file will be found in the Downloads directory within your home directory, i.e. "~/Downloads".

   ```python
   $ cd <downloads-directory>
   ```

   3.3 Start by copying the EXOTIC compressed files from where you downloaded them into the directory for running EXOTIC that you created above. This command will copy the software package, where "n.nn.n" is the number of your release of EXOTIC in the downloaded file name and **<exotic>** is the name of the directory that you created above.

      ```python
   $ cp EXOTIC-n.nn.n.zip ~/<exotic>
      ```

   3.4 Now you need to step into the directory where the EXOTIC files have just been copied and get ready to extract and decompress the software. Replace **<exotic>** with the directory that you created above.

      ```python
   $ cd ~/<exotic>
      ```

   3.5 Then the EXOTIC software needs to be extracted from the archive file into the directory where you are now located.

   ```python
   $ unzip EXOTIC-n.nn.n.zip
   ```

   3.6 Now you need to step into the directory where the EXOTIC files have just been placed before you can run the software. The actual directory name will reflect both the directory in which you originally chose to place your files and the version of EXOTIC that you have downloaded.

      ```python
   $ cd ~/<exotic>/EXOTIC-n.nn.n
      ```

   3.7 Run this script to verify that you have suitable versions of Python and pip and then install all the necessary packages to run EXOTIC. Enter the password you created above if prompted.

      ```python
   $ ./run_exotic_linux.sh
      ```

   3.8 Wait a while. This process will take up to around 15 minutes to download and install the complete EXOTIC system.

   A terminal window with a black background will then open up, showing you the progress of the installation for each of the files required to run EXOTIC. This process may take some time and you will see many lines of yellow and white text scroll by until the whole installation is finished. Do not worry when the process pauses on a line starting with "Building wheel for ...", since each of these steps will take some time. 

   3.9 The final stage, beginning "Installing collected packages:..." followed by a long list of packages, will take several minutes to run and the screen will not change at all during this time - please be patient and just wait for the installation to complete on its own.

   When the script completes, a new command line window in black appears that will show the progress of the data-reduction once EXOTIC is running. In front of this you will see the EXOTIC GUI window showing "Welcome to EXOTIC!" plus the version of EXOTIC you are running and other details. This GUI window is used to enter the values required to define how the data-reduction run will be performed.

   You can move each window around separately to see what is happening in each of them.

## Running EXOTIC next time

To run the EXOTIC program in the future, simply click locate and run the "run_exotic_linux.sh" file again. This will firstly check whether any EXOTIC files need updating and update them automatically if there has been another release of the software. In this case you will see the same installation process as outlined above, though it will run very much faster, as most files will not have changed between releases. 

After the installation of any new or updated EXOTIC components has completed (or is skipped through if nothing has changed) the program will run as before.

To learn more about how EXOTIC works, check our other guides on the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC/tree/main/Documentation). Once you have opened the document you wish to read, right click the "Download" button to save the PDF file, and then open it in your normal PDF reader, so that you can click the links in the document.

## Further Help and Support

If you have any questions or comments, please feel free to reach out to us on Slack or email at [exoplanetwatch@jpl.nasa.gov](exoplanetwatch@jpl.nasa.gov)



