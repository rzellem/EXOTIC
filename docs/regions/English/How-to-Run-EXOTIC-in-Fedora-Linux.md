# Installation of EXOTIC in Fedora Linux

## 1 Introduction
The intent of these instructions is to guide the user in installing **EXOTIC** on ***Fedora 33*** (or later).  Two methods are provided, based on user preference.  These instructions assume a default ("*vanilla*") installation of ***Fedora*** Linux, yet may be modified as required for custom changes made by the user subsequently.

## 2 Installation of EXOTIC with virtualenv Tool
This method makes use of a virtual Python Environment builder to create an environment dedicated to **EXOTIC** and isolated from the native ***Fedora*** installation of Python.

### 2.1 Set Up Virtual Environment
1. Ensure that **Python 3.7** or later is installed

   `python3 --version`

2. Ensure that ***virtualenv*** and supporting tools are installed

   `sudo dnf install python3-virtualenvwrapper python3-virtualenv python3-virtualenv-clone`
   
   *Note: This command only performs the installation of these packages if not already installed (and may be ommited in this case)*
   
3. Edit ***.bashrc*** file

   `echo "export WORKON_HOME=~/env" >> ~/.bashrc`

   Run the modified script:
   
   `source ~/.bashrc`

4. Create directory for ***virtualenv***-compatible Python virtual environments

   `mkdir -p $WORKON_HOME` 

5. Initialize ***virtualenvwrapper***

   `/usr/bin/virtualenvwrapper.sh`
   
   *Note: This command creates various new directories within the newly created Python virtual environment directory. Optionally, the command may be added to the ***.bashrc*** file which executes each time a terminal window is launched:*
   
   `echo "source /usr/bin/virtualenvwrapper.sh" >> ~/.bashrc`

6. Create new Python virtual environment for **EXOTIC**

   `mkvirtualenv exotic`
   
   *Note: The virtual Python environment for* **EXOTIC** *is now created and can be confirmed by the* "***(exotic)***" *prefix in the terminal:*
   ```bash
   (exotic) user@localhost:~>
   ```
      
### 2.2 Install EXOTIC

7. Ensure that the **EXOTIC** virtual environment is active, verified by the "**(exotic)**" prefix in the terminal window

8. Install pre-requisite Python packages

    `pip install distlib numpy cython scikit-image`

9. Issue the installation command

    `pip install exotic`


### 2.3 Run EXOTIC 

10. If the **EXOTIC** virtual environment is ***not*** active, such as after a computer reboot or when launching a new terminal window, issue the following command:

    `workon exotic`

    *Note: If* ***virtualenvwrapper*** was not added to the ***.bashrc*** *file per Step 5, the above command must be preceded by the* `source /usr/bin/virtualenvwrapper.sh` *command*

11. Ensure that the **EXOTIC** virtual environment is active, verified by the "**(exotic)**" prefix in the terminal window

12. Run **EXOTIC**

    `exotic`

### 2.4 Exiting EXOTIC Virtual Environment

13. This final command is optional, and is to be invoked only if the user desires to exit the **EXOTIC** virtual environment completely:

    `deactivate`

    *Note: After this command is executed, perform Steps 10 - 12 to run* **EXOTIC**. 

## 3 Installation of EXOTIC with Miniconda 

This is an alternate method that makes use of ***Miniconda*** to create an environment dedicated to **EXOTIC** and isolated from the native ***Fedora*** installation of Python.

### 3.1 Install Miniconda

14. Download a ***Miniconda*** [Linux Installer](https://docs.conda.io/en/latest/miniconda.html#linux-installers) compatible with your computer

15. Install ***Miniconda*** with the following command and follow the script's installation instructions:

    `bash Miniconda3-latest-Linux-x86_64.sh`

    *Note: Replace "Miniconda3-latest-Linux-x86_64.sh" with actual downloaded filename.  Original installation instructions are found at* [conda.io](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### 3.2 Set Up Virtual Environment

16. Ensure that the *base* ***Miniconda*** virtual environment is active, verified by the "**(base)**" prefix in the terminal window

    ```bash
    (base) user@localhost:~>
    ```
    *Note: Depending on the user's selections during installation, the* ***Miniconda*** *virtual environment may activate by default in the terminal window or may require manual activation*

17. Create dedicated virtual environment for **EXOTIC**

    `conda create --clone base --name exotic`

18. Activate the newly created **EXOTIC** virtual environment in preparation for its installation

    `conda activate exotic`

19. Ensure that the *exotic* ***Miniconda*** virtual environment is active, verified by the "**(exotic)**" prefix in the terminal window

    ```bash
    (exotic) user@localhost:~>
    ```
### 3.3 Install EXOTIC

20. Install pre-requisite ***Fedora*** packages

    `sudo dnf install gcc gcc-c++`

21. Install pre-requisite Conda packages

    `conda install numpy cython scikit-image`

22. Issue the installation command

    `pip install exotic`

### 3.4 Run EXOTIC 

23. If the **EXOTIC** virtual environment is ***not*** active, such as after a computer reboot or when launching a new terminal window, issue the following command, otherwise skip to the next step:

    `eval "$(~/miniconda3/bin/conda shell.bash hook)"`

    `conda activate exotic`

24. Ensure that the **EXOTIC** virtual environment is active, verified by the "**(exotic)**" prefix in the terminal window

25. Run **EXOTIC**

    `exotic`

### 3.5 Exiting EXOTIC Virtual Environment

26. This final command is optional, and is to be invoked only if the user desires to exit the EXOTIC virtual environment completely:

    `conda deactivate`

    *Note: Deactivating the* "***(exotic)***" *virtual environment places the user in the* "***(base)***" *environment.  Running the command above a second time places the user in the default system environment.  After this command is executed, perform Steps 23 - 25 to run* **EXOTIC**. 
