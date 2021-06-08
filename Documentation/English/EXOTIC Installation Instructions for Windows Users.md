# EXOTIC Installation Instructions for Windows Users

## Install DS9 (Astronomical Image Viewing Software)

Before installing EXOTIC itself, you may wish to install the SAOImageDS9 (Astronomical Image Viewing Software) to view the “.FITS” images you obtain during observations. For more information on SAOImageDS9, check out the [User Guide](http://ds9.si.edu/doc/user/index.html)


- Follow [SAOImageDS9 download link](https://sites.google.com/cfa.harvard.edu/saoimageds9/download)
- Click on "***Windows***" then select "***Windows 64 bit***".
- Run the Windows installer in the usual way once the file has been downloaded.
- Follow the instructions in the installer to complete the installation.

## Running EXOTIC
There are two ways you can run EXOTIC: either through a web browser interface with a so-called Python Notebook using the Google Colab (or by running Jupyter Notebook) or locally on your own personal computer using the command line. 

Read the details of these two approaches below and select the best option for you:

### Google Colab:

Google Colab is a free facility (available to anyone with a Google account) where you may write and run Python programs in your browser without installing anything on your local machine. Running EXOTIC in the Google Colab is the recommended option for all users as it offers a more user-friendly and interactive experience. It is especially recommended for new EXOTIC users, those who are unfamiliar with using the command line (where you have to type in commands - instead of using the mouse), students and those analysing data via remote observatories. 

To use Google Colab to run EXOTIC follow the links on the home page of the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC).

Once you are happy running running EXOTIC with the Colab and you want to run your own Python Notebook, please visit the instructions [‘How-to-Run-EXOTIC-with-the-Python-Notebook](https://github.com/rzellem/EXOTIC/raw/main/Documentation/English/How-to-Run-EXOTIC-with-the-Python-Notebook.pdf)’ in the [Documentation folder] (https://github.com/rzellem/EXOTIC/tree/main/Documentation) on GitHub.

*Note: in order to be able to click on the links and select text in this PDF document, you must download it from GitHub. The GitHub preview simply shows you an image of the document, which does not allow for those functions.*

### Command line - for cmd.exe or PowerShell Users:

If you are comfortable using the command line in Windows (with a windows terminal tool such as cmd.exe or Powershell), instead of using your mouse to navigate around a web page, you may prefer to run EXOTIC on your own Windows machine in order to:

- reduce large data sets
- learn more about how EXOTIC works and get access to the full range of options
- carry out a real time data-reduction while observing on the telescope

To get a crash course in how to use the command line in Ubuntu and how it differs from the command line in Windows,  go to the ***60-second guide to using the Ubuntu terminal*** at the foot of this article. 

***Note: The instructions below tell you how to install WSL2 and Ubuntu on your Windows computer before installing the Exoplanet Watch software. EXOTIC can be run on your Windows computer without installing Ubuntu. We do not recommend this approach as the installation process is much more difficult and EXOTIC runs much slower on Windows.***

1. #### Windows Subsystem for Linux (WSL2) and Ubuntu


Windows is one example of a Microsoft operating system, which is the software that manages your computer hardware and software. It supports interfaces for the user (such as keyboards and mice for input and graphical displays and speakers for output) as well as important utilities for managing the computer's  memory, processor and disk drives.

WLS2 is a new product from Microsoft that sits between Windows and your machine's hardware, so that you may then install another operating system in parallel with Windows that will run side by side. Ubuntu running on WSL2 behaves just like any other windows application. It runs in its own window, where you can execute the commands necessary to run EXOTIC. While Ubuntu is running on WSL2, you can at the same time still use your usual Windows web browser, the file manger or any other Windows application.

Ubuntu is one example of a Linux operating system. Linux underpins over 95% of the servers that power the world wide web and 100% of the world's top supercomputers. Installing Ubuntu with WSL2 will simply allow you to make use of this operating system when (and only when) you are using the app. It should not damage your current Windows setup (but you install it at your own risk) and it will make running EXOTIC a lot easier and faster!

2. #### Install Windows Subsystem for Linux (WSL2) and Ubuntu 20.04

Go to the Microsoft page: [Windows Subsystem for Linux Installation Guide for Windows 10](https://docs.microsoft.com/en-us/windows/wsl/install-win10) where you will use the special Windows command line tool called PowerShell to install WSL2. Ignore the first section "***Simplified Installation for Windows Insiders***" (unless you have already joined this program and wish to experiment with new upcoming versions of the WSL2 installation process).  Instead, follow the first five steps in the "***Manual Installation Steps***" section, restarting your machine when instructed to do so. 
Finally, execute '***Step 6 - Install your Linux distribution of choice***' selecting "***Ubuntu 20.04 LTS***" from the choice of Linux distributions in the Microsoft Store. This is a stable version of Ubuntu that was released a year ago and will be supported until 2025.

3. #### Start Ubuntu and open a terminal window 

- Click the Windows start button in the lower left-hand corner of your screen.
- Scroll down to “***U***” in the list of applications or type “***Ubuntu***”.
- Click on "***Ubuntu***" to start the application - this opens a so-called "terminal window" where you can type in commands.
- When opening Ubuntu for the first time, you will be prompted to create an administrator account. Make up and enter a new username and password for your Ubuntu user (not your Windows username/password) and take careful note of them. You will be prompted again to enter the password  later when executing certain commands in Ubuntu that start with the word "sudo".

- If you wish to leave Ubuntu, you can just click on any other running application in Windows, or start a new application in the usual way. If you wish to close the Ubuntu window, just enter exit and press return or click on the close window button in the ordinary way.


4. #### Install EOXTIC

Execute the commands shown in the pink boxes below by typing them on the Ubuntu command line after the prompt symbol "$". After entering each command you will need to press the Enter (or Return) key and wait for it to complete, which is when the prompt reappears.

Start by updating your Ubuntu system to the latest version and then proceed install the EXOTIC program.

​    4.1. Update the repository of Ubuntu software (equivalent to Microsoft store).  Enter the password you created above when prompted. This command will show you all the packages that can be upgraded and any available essential new packages and then download them.

```python
  $ sudo apt update
```

​     4.2. Update each package  to its latest version. You will see a list of the packages being upgraded or installed, as well as a progress bar showing the progress of the overall upgrade to your Ubuntu system. Enter 'y' or 'Y' to continue with installation when prompted.

```python
   $ sudo apt upgrade
```
   4.3. Create the folder you want to locate the EXOTIC files in, replacing ***folder*** with the directory of your choice. If you wish to be able to access the results files from windows, you must choose to place this folder in the windows filing system. To do so use the prefix "/mnt/c/Users/your_windowsname" and replace “your_windowsname” with the username for the account you are signed into on Windows (not your Ubuntu username). 

For example, typing “mk /mnt/c/Users/your_windowsname/exoplanets” will mean that the EXOTIC files will be stored in "C:\Users\your_windowsname\exoplanets" in Windows.


```python
   $ mkdir /folder
```

   4.4. Switch to your chosen directory and get ready to install the software – replacing ***folder*** with the directory you created above (for example “cd /mnt/c/Users/your_windowsname/exoplanets”). This will also be the location from where you run the program.

```python
$ cd /folder
```
   4.5. Install Python 3 (or update Python 3 to the latest version) on Ubuntu. You may be prompted for the password you created above.

```python
$ sudo apt install python3
```
Running these next few commands will take some time (probably a few minutes). Wait for each one to complete (when you will see the prompt appear again) before proceeding to the next instruction.

   4.6.  Next install pip, which is a tool needed to install the rest of the EXOTIC software.

```python
$ sudo apt install python3-pip
```

   4.7. When the installation is complete, verify the installation by checking the pip version - enter:

```python
   $ pip3 --version
```
The version number may vary, but it will look something like this:
```output
$ pip 20.0.2 from /usr/lib/python3/dist-packages/pip (python 3.8)
```

   4.8. Install EXOTIC and all necessary Python packages used to run EXOTIC.

```python
   $ pip3 install exotic
```
   4.9. You will see a list of the various packages required for EXOTIC being collected and then installed. This process will take around 15 minutes or more and you will see many lines of output before the prompt reappears and you can enter further commands.

*Note: If this command fails, you can manually download EXOTIC from GitHub. To do this, follow the instructions below in the section entitled "**Alternative method for Installation of EXOTIC**".*

   4.10. Start up the EXOTIC program.

```python
$ exotic
```

As soon as it has started you should see the introductory header to EXOTIC as shown below, which confirms that it is all up and running! Do not worry if you see a different version of Python or of EXOTIC. 

```output

$ *************************************************************
$ Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)
$ Version 0.42.2
$ *************************************************************

```

If you are seeing this header in your Ubuntu terminal, you have successfully installed EXOTIC and can now use it at any time to reduce the data from your amazing transit observations! Do not worry if you see a different version of Python or of EXOTIC, since new versions are coming out all the time.

To learn more about how EXOTIC works, check our other guides on the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC/tree/main/Documentation). Once you have opened the document you wish to read, right click the "Download" button to save the PDF file, and then open it in your normal PDF reader so that you can click the links in the document.

5. #### Alternative method for Installation of EXOTIC

If you experienced issues with the downloading and installing EXOTIC above, please follow these steps instead. Note that they take you right back to the start of the process and run the steps above in one continuous set of commands. 

   5.1. Start by copying the EXOTIC project files from GitHub, the repository for EXOTIC. This command will install the software into a directory called "EXOTIC" within the directory you created above.

```python
   $ git clone https://github.com/rzellem/EXOTIC.git
```

   5.2 Now you need to step into the directory where the EXOTIC files have just been copied and get ready to install the software. Replacing ***folder*** with the directory you created above (for example “cd /mnt/c/Users/your_windowsname/exoplanets/EXOTIC”). This will also be the location from where you run the program.

```python
$ cd /path/EXOTIC
```

  5.3 One of the files you have downloaded, named “exotic_installation_linux.sh”, is called a "shell script". Such a script consists of a list of commands to be executed in the Ubuntu terminal. Modify this file to make it executable, meaning that you can run it in your terminal.

```python
   $ sudo chmod 755 exotic_installation_linux.sh
```

   5.4 Run this script to download or update Python (even if you already have it) and install all the necessary packages to run EXOTIC. Enter the password you created above if prompted.

Finally, the script will run EXOTIC for the first time to test that it is functional.

```python
   $ ./exotic_installation_linux.sh
```
   5.5 Wait a while. This process will take up to around 15 minutes to download and install the complete EXOTIC system. 

Once the process has completed, you should see the introductory header to EXOTIC as shown below, which confirms that it is all up and running! Do not worry if you see a different version of Python or of EXOTIC, since new versions are coming out all the time. 

```output

$ *************************************************************
$ Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)
$ Version 0.42.2
$ *************************************************************

```

To learn more about how EXOTIC works, check our other guides on the [EXOTIC Project on GitHub](https://github.com/rzellem/EXOTIC/tree/main/Documentation). Once you have opened the document you wish to read, right click the "Download" button to save the PDF file, and then open it in your normal PDF reader, so that you can click the links in the document.

To run the EXOTIC program next time around make sure that you have switched to the correct folder 

```python
$ cd /mnt/c/Users/your_windowsname/your_folder/EXOTIC
```

 Then to start reducing your data with EXOTIC just use the command:

```python
   $ exotic
```

##### *60-second guide to using the Ubuntu terminal*

*Many installations of Ubuntu come with a Graphical User Interface (GUI) just like Windows, where you can use the mouse as well as the keyboard to start up and interact with programs. But just using the command line on Ubuntu allows you to perform a wide range of actions on your computer (run Python programs, install Linux applications, edit files, etc.) by typing in commands. If you are interested in learning more about the command line on Ubuntu and the wide range different commands you can use, follow this link: [The Linux command line for beginners](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview). This document refers to Ubuntu version 18.04 but there's no command that you can use in Ubuntu 18.04 that will not work in 20.04.*

*If you make a mistake when entering a command in the Ubuntu terminal, do not try to use the mouse to help you edit the command line. You can only use the left and right cursor keys to move through the line and the Delete and Back keys to remove characters. Be aware that in all Linux systems, including Ubuntu, all commands (as well as file and directory names) are case sensitive. So typing "SUDO apt upgrade", instead of "sudo apt upgrade" to upgrade your Ubuntu installation, will result in an error message.*

*Everything that you type and the responses to your commands are stored in the terminal window, even though they may not be visible because they have scrolled off the top of the window. So if you want to repeat a command, then you can use the up and down cursor keys to move through your list of previous commands. Stop at the one you require and just press Enter. Another feature is that rather that entering a long directory/file name in full, if you start typing and then press the "Alt" key, Ubuntu will try to fill in the rest of the command for you. If two directories/files start with the characters you have entered, the completion mechanism will stop, to let you enter more characters to distinguish between them.*

*If you wish to stop a program running for any reason, (perhaps you realise that you made a mistake in entering your data), then hold down the CTRL key and press the "c" key. This will stop the program almost immediately and return you to the prompt ready for you to enter another command.*

*If you wish to copy some text (like a directory name) from Windows or from your Ubuntu terminal window, first highlight what you wish to copy with your mouse, then hold down the Ctrl key and press the "c" key to copy. Then go to the Ubuntu terminal window, move to where you want to paste the text using the cursor keys and press the right button on your mouse (Ctrl "v" will not work).*

*Ubuntu and other flavours of Linux generally refer to "directories", while Windows uses the term "folders" but these terms mean exactly the same thing and can used interchangeably. You can see the contents of your Ubuntu directory with the "ls" command (similar for the Windows "dir" command). Changing to another directory (or folder) is done with the "cd" command. In executing “cd ~/exoplanets” for example, you are navigating to your "exoplanets" folder, just as you would by double-clicking on "exoplanets" in Windows File Explorer.*

*A single dot prefix (e.g. ./filename) refers to a file in the directory you are in, while a double dot prefix (e.g. ../filename) refers to a file in the directory above. So typing "cd .." when you are in "/mnt/c/Users/windows_username/exoplanets” takes you to "/mnt/c/Users/windows_username".*

*You will see that when the Ubuntu terminal is waiting for you enter a command, it will show a prompt like this: "ubuntu_username@machine _name: directory". This is made up from the username for the account you are signed into on Windows, the name of your computer and an Ubuntu directory with @ sign and a colon between the three parts. The first part of the Ubuntu directory "/mnt/c/Users/windows_username/” represents your home folder in Windows and the remainder shows which directory you have moved into using the change directory command. This means that you are actually pointing to the windows filing system, in your normal user space so the folders you create and the files that you save in Ubuntu can be seen from Windows itself.* 

*Ubuntu labels your home folder in Windows as"/mnt/c/Users/windows_username/” (note the lower case "c" for the "C" drive) whereas Windows labels it "C:\Users\your_username" using back slashes not forward slashes. If you compare the contents using Windows file manager, you can see that these are in fact one and the same folder with the same content. In the Ubuntu filing system a directory named "EXOTIC" is quite separate from one name "exotic but in the Windows filing system you are in fact two references to the same folder. Be aware that if you create or delete a file or folder in your EXOTIC folders with Windows File Explorer, it will be added or removed from Ubuntu and vice-versa.* 

## Further Help and Support

If you have any questions or comments, please feel free to reach out to us on Slack or email at [exoplanetwatch@jpl.nasa.gov](exoplanetwatch@jpl.nasa.gov)

