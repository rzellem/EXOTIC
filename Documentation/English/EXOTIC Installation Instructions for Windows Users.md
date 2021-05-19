# EXOTIC Installation Instructions for Windows Users

~<!--Images of running programs and logos have yet to be restored-->~

## Install DS9 (Astronomical Image Viewing Software)

Before installing EXOTIC itself, you may wish to install the SAOImageDS9 (Astronomical Image Viewing Software) to view the “.FITS” images you obtain during observations. For more information on SAOImageDS9, check out the [User Guide](http://ds9.si.edu/doc/user/index.html)

- Follow [SAOImageDS9 download link](https://sites.google.com/cfa.harvard.edu/saoimageds9/download)
- Click on "Windows" then select "Windows 64 bit".
- Run the Windows installer in the usual way once the file has been downloaded.
- Follow the instructions in the installer to complete the installation.

## EXOTIC
There are two ways you can run EXOTIC: either through a web browser interface with a so-called Python Notebook using the Google Colab (or by running Jupyter Notebook) or locally on your own personal computer using the command line. 

Read the details of these two approaches below and select the best option for you:

### Google Colab - RECOMMENDED:

Google Colab is a free facility (available to anyone with a Google account) where you may write and run Python programs in your browser without installing anything on your local machine. Running EXOTIC in the Google Colab is the recommended option for all users as it offers a more user-friendly and interactive experience. It is especially recommended for new EXOTIC users, those who are unfamiliar with using the command line (where you have to type in commands - instead of using the mouse), students and those analysing data via remote observatories. 

To use Google Colab to run EXOTIC follow the links on the home page of the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC).

Once you are happy running running EXOTIC with the Colab and you want to run your own Python Notebook, please visit the instructions [‘How-to-Run-EXOTIC-with-the-Python-Notebook](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/How-to-Run-EXOTIC-with-the-Python-Notebook.pdf)’ in the [Documentation folder] (https://github.com/rzellem/EXOTIC/tree/main/Documentation) on GitHub.
~~*Please note: in order to be able to click on the links and select text in this PDF document, you must download it from GitHub. The GitHub preview simply shows you an image of the document, which does not allow for those functions.*~<!--I think ths note will not be applicable if the Notebook instructions are rewritten using Markdown - at which point I shall remove it-->~

### Command line:

If you are comfortable using the command line in Windows (with a windows terminal tool such as cmd.exe or powershell), instead of using your mouse to navigate around a web page, you may prefer to run EXOTIC on your own machine in order to:

- reduce large data sets
- learn more about how EXOTIC works and get access to the full range of options
- carry out a real time data-reduction while observing on the telescope

*Note: The instructions below tell you how to install WSL2 and Ubuntu on your Windows computer before installing the Exoplanet Watch software. EXOTIC can be run on your Windows computer without installing Ubuntu. However, we do not recommend this approach as the installation process is much more difficult and EXOTIC runs much slower on Windows. If you would still prefer running EXOTIC natively on Windows (instead of in Ubuntu), follow this guide instead: https://github.com/rzellem/EXOTIC/blob/master/Documentation/EXOTIC-Instructions.pdf*

1. #### What are the Windows Subsystem for Linux (WSL2) and Ubuntu?


Windows is one example of a Microsoft operating system, which is the software that manages your computer hardware and software. It supports interfaces for the user (such as keyboards and mice for input and graphical displays and speakers for output) as well as important utilities for managing the computer's' memory, processor and disk drives.

WLS2 is a new product from Microsoft that sits between Windows and your machine's hardware, so that you may then install another operating system in parallel with Windows that will run side by side. Ubuntu running on WSL2 behaves just like any other windows application. It runs in its own window, where you can execute the commands necessary to run EXOTIC. While Ubuntu is running on WSL2, you can at the same time still use your usual Windows web browser, the file manger or any other Windows application.

Ubuntu is one example of a Linux operating system. Linux underpins over 95% of the servers that power the world wide web and 100% of the world's top supercomputers. Installing Ubuntu with WSL2 will simply allow you to make use of this operating system when (and only when) you are using the app. It should not damage your current Windows setup (but you install it at your own risk) and it will make running EXOTIC a lot easier and faster!

INSERT IMAGE OF UBUNTU RUNNING ON WINDOWS.

2. #### Install Windows Subsystem for Linux (WSL2) and Ubuntu 20.04

Go to the Microsoft page: [Windows Subsystem for Linux Installation Guide for Windows 10.](https://docs.microsoft.com/en-us/windows/wsl/install-win10) Ignore the first section "Simplified Installation for Windows Insiders" (unless you have already joined this program and wish to experiment with new upcoming versions of the WSL2 installation process).  Instead, follow the first five steps in the "Manual Installation Steps" section, restarting your machine when instructed to do so. 
Finally, execute 'Step 6 - Install your Linux distribution of choice' by selecting Ubuntu version "20.04 LTS". This is a stable version of Ubuntu that was released a year ago and will be supported until 2025.

3. #### Start Ubuntu and open a terminal window 

- Click the Windows start button in the lower left-hand corner of your screen.
- Scroll down to “U” in the list of applications or type “Ubuntu”.
- Click on "Ubuntu" to start the application - this opens a so-called "terminal window" where you can type in commands.
- When opening Ubuntu for the first time, you will be prompted to create an administrator account. Make up and enter a new username and password for your Ubuntu user (not your Windows username/password) and take careful note of them. You will be prompted again to enter the password  later when executing certain commands in Ubuntu that start with the word "sudo".

- If you wish to leave Ubuntu, you can just click on any other running application in Windows, or start a new application in the usual way. If you wish to close the Ubuntu window, just enter exit and press return or click on the close window button in the ordinary way.


4. #### Install EOXTIC

~<!--This section has only been started and should not be read at this stage. See below for a new 30 second guide to using the Ubuntu terminal-->~

 ~<!--Execute the following commands in Ubuntu to update your system to the latest version and then install the EXOTIC program. After entering each command you will need to press the Enter (or Return) key.
Type “sudo apt update” – do not include the quotes.
Enter the password you just created when prompted. This will then show you all the packages that can be upgraded to their latest version.
Type “sudo apt upgrade” – do not include the quotes. You will see a list of the packages being upgraded, as well as a progress bar showing the progress of the overall upgrade.
Type “mk /PATH/” – replacing PATH with the folder you want to create to hold the EXOTIC folder and files to be downloaded.   For example, typing “mk Documents” will mean that the EXOTIC files can  be stored in "/mnt/c/Users/your_username/Documents/EXOTIC/".
Type “cd /PATH/” – replacing PATH with the directory you created above to switch to this dorectory.
This will also be the location in which you run the program.
  o Note: cd stands for “Change Directory”. In executing this command, you are navigating to your Downloads folder, just as you would by double-clicking on Downloads in File Explorer.
  • Install Python (or update Python to the latest version) on Ubuntu by using the command ‘sudo apt-get install python3’. Here you will be prompted for the password you created in the beginning (or already had if you’ve used Ubuntu before).
  Running these next few commands will take some time (probably a few minutes). Wait for them to complete.
  • Type ‘curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py’ and hit enter.
  • Type ‘python3 get-pip.py’ and hit enter. If this command fails, try replacing ‘python3’ with ‘python’.
  • Run the command ‘pip3 install exotic || pip install exotic’. This command will install EXOTIC and all necessary Python packages used to run EXOTIC.
  o Please Note: If this command fails, you can manually download EXOTIC off GitHub. To do this, follow the additional instructions below.
  • Wait a while. This process may take several minutes.
  • Type ‘exotic’ and hit enter. You should see the introductory header to EXOTIC as pictured at the bottom of this document, which tells you that it is all up and running!
  If you experienced issues with the starred command above, please follow these next steps:
  • Type “git clone https://github.com/rzellem/EXOTIC.git” and hit enter.
  • Type “chmod 755 exotic_installation_linux.sh” -- do not include the quotes. Hit Enter.
  o Note: this command alters the file you downloaded “exotic_installation_linux.sh” to be executable (i.e. you can now run it in your terminal).
  • Type “./exotic_installation_linux.sh” -- do not include the quotes.
  • Hit Enter.
  • Enter the password you created earlier (or already had) when prompted.
  o Note: this command runs the file you downloaded, which is called a script. A script is simply a list of commands to
  be executed in the Ubuntu. This script will download Python (unless you already have it) and install all the necessary packages to run EXOTIC. Finally, the script will run EXOTIC to test that it is functional.
  • Wait a while. This process may take several minutes.
  • Once the process has completed, you should see the introductory header to EXOTIC as pictured at the bottom of this document, which tells you that it is all up and running!
  If you are seeing this header in your Ubuntu terminal, EXOTIC is successfully installed, and you are ready to begin analyzing exoplanet data!
  And that’s it! You have successfully installed EXOTIC and can now use it at any time to reduce the data from your amazing transit observations!
  To learn how to run the code and how EXOTIC works, check our other guides on GitHub!
  → https://github.com/rzellem/EXOTIC/tree/main/Documentation
  end-->~

  

##### *30-second guide to using the Ubuntu terminal*

*You will see that when the Ubuntu terminal is waiting for you enter a command, it will show a prompt line  "ubuntu_username@machine _name: /mnt/c/Users/windows_username/”– where “windows_username” is the username for the account you are signed into on Windows and ubuntu_username is that for Ubuntu. This means that you are actually pointing to the windows filing system, in your normal user account space where you can create folders and save files that you can see from Windows itself. (Windows labels this C:\Users\your_username with back slashes not forward slashes but if you look there using Windows file manager you can see that these are in fact one and the same folder with the same content). Be aware that if you create or delete a file or folder here in Windows, it will be added or removed from Ubuntu and vice-versa.* 

*Note: Ubuntu and other flavours of Linux generally refer to "directories", while Windows uses the term "folders" but these terms mean exactly the same thing and can used interchangeably.*

*If you make a mistake when entering a command in Ubuntu, do not try to use the mouse to edit the command line. You can only use the left and right cursor keys to move through the line and the Delete and Back keys to remove characters. If you want to repeat a command then you can use the up and down cursor keys to move through your list of previous commands. Stop at the one you require and just press Enter.* 

*If you wish to copy some text (like a directory name) from Windows or from your Ubuntu terminal window, then highlight what you wish to copy with your mouse, then hold down the Ctrl key and press the "c" key to copy. Then go to the Ubuntu terminal window, move to where you want to paste the text and press the right button on your mouse (do not use Ctrl "v").*

*Note: Ubuntu allows you to perform a wide range of actions on your computer (run Python programs, install Linux applications, edit files, etc.) by typing in commands. If you are interested in learning more about the command line on Ubuntu and the wide range different commands you can use, follow this link: [The Linux command line for beginners](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview). This document refers to Ubuntu version 18.04 but there's no command that you can use in Ubuntu 18.04 that will not work in 20.04.*

## Further Help and Support

If you have any questions or comments, please feel free to reach out to us on Slack or email at exoplanetwatch@jpl.nasa.gov