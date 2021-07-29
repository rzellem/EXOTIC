          EXOTIC Installation Instructions for Mac Users

 There are two ways you can run EXOTIC, either with a Python Notebook in Jupyter or the Google Collab or on your own personal machine using the command line. Please review the details of the two processes below and select the best option for you.

 Google Collab: Running EXOTIC in the Google Collab is the recommended option for all users as it offers a more user-friendly and interactive experience. It is especially recommended for new users, those who are unfamiliar with using the command line, students, and those analyzing data via remote observatories. To run EXOTIC with the Python Notebook, please visit the instructions ‘How-to-Run-EXOTIC-with-the-Python-Notebook’ in the Documentation folder.

  Command line: If you feel comfortable using the command line, are reducing large data sets, or are looking to learn more about EXOTIC and how to run it on your own machine, follow the set of instructions on the following pages:

  Please note: in order to be able to click on the links and select text in this document, you must download it off GitHub. The GitHub preview simply shows you an image of the document, which does not allow for those functions.

I. Download DS9 (Astronomical Image Viewing Software) • Follow the link: https://sites.google.com/cfa.harvard.edu/saoimageds9?pli=1&authuser=1 • Download the version corresponding to the Mac operating system. • Run the installer once downloaded. SAOImageDS9 (@SAOImageDS9) | Twitter• Follow the instructions in the installer to complete the installation. 

  Note: This software will be used to view the “.FITS” images you obtain during observations. For more information on DS9, check out the User Guide: http://ds9.si.edu/doc/user/index.html Change your Mac Hostname via Terminal | OSXDaily

  II. Open the Terminal app • Select Launch Pad on the Task Bar • Double-click Terminal • When opening Terminal for the first time, you might be prompted to enter a password. Enter your new password and take note of it. You will be prompted again later to enter it. If you are not prompted to add one, your password will be the same as the password for the computer account you are currently working on.

  Note: The Terminal app allows you to perform actions on your computer (run python programs, install applications, edit files, etc.) by typing in commands. If you are interested in learning more about the terminal and the different commands you can use, follow this link: https://www.macworld.com/article/2042378/master-the-command-line-navigating-files-and-folders.html

III. Installing EXOTIC - Execute the following commands in your Terminal’s command line: 

   • Type “cd /PATH/” – replacing PATH with the directory you want the EXOTIC folder to be downloaded in. For example, typing “cd Documents” will mean that the EXOTIC files will be stored in Users/your\_username/Documents/EXOTIC/.

   • This will also be the location in which you run the program.

   o Note: cd stands for “Change Directory”. In executing this command, you are navigating to your Downloads folder, just as you would by double-clicking on Downloads in File Explorer.

   • Install Python (or update Python to the latest version) on your Mac by typing the command ‘/bin/bash -c "\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)’ and hitting enter. Do not include the single quotes.

   • Type “brew install python” and hit enter. If you already have Python installed, you can skip this step.

   • Type “curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py” and hit enter.

   • Type “python get-pip.py” and hit enter.

    • Run the command ‘pip3 install exotic || pip install exotic’. This command will install EXOTIC and all necessary Python packages used to run EXOTIC.

    Please Note: If this command fails, you can manually download EXOTIC off GitHub. To do this, follow the additional instructions below.

    • Wait a while. This process may take several minutes.

    • Type ‘exotic’ and hit enter. You should see the introductory header to EXOTIC as pictured at the bottom of this document, which tells you that it is all up and running!

    If you experienced issues with the starred command above, please follow these next steps:

    • Type “git clone https://github.com/rzellem/EXOTIC.git” and hit enter.

    • Type “chmod 755 exotic\_installation\_mac.sh” – do not include the quotes. Hit Enter.

    ▪ Note: this command alters the file you downloaded “exotic\_installation\_mac.sh” to be executable (i.e. you can now run it in your terminal).

    • Type “./exotic\_installation\_mac.sh” – do not include the quotes. Hit Enter.

    • Enter the password you created earlier (or already had) if prompted.

      ▪ Note: this command runs the file you downloaded, which is called a script. A script is simply a list of commands to be executed in the Terminal. This script will download Python (unless you already have it) and install all the necessary packages to run EXOTIC. Finally, the script will run EXOTIC to test that it is functional.
    • Wait a while. This process may take several minutes.

    • Once the process has completed, you should see the introductory header to EXOTIC as pictured at the bottom of this document, which tells you that it is all up and running! ----------------------------------------------------------------------------------------------------


    If you are seeing this header in your terminal, EXOTIC is successfully installed, and you are ready to begin analyzing exoplanet data!





     And that’s it! You have successfully installed EXOTIC and can now use it at any time to reduce the data from your amazing transit observations!

     To learn how to run the code and how EXOTIC works, check our other guides on GitHub!

      → https://github.com/rzellem/EXOTIC/tree/main/Documentation

      If you have any questions or comments, please feel free to reach out to us on Slack or email at exoplanetwatch@jpl.nasa.gov
