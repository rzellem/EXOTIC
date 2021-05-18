# Installation of EXOTIC in openSUSE Linux

## 1 Introduction
The intent of these instructions is to guide the user in installing **EXOTIC** on ***openSUSE Tumbleweed*** or ***openSUSE RC13*** (or later).  Two methods are provided, based on user preference.  These instructions assume a default ("*vanilla*") installation of ***openSUSE*** Linux, yet may be modified as required for custom changes made by the user subsequently.

## 2 Installation of EXOTIC with virtualenv Tool
This method makes use of a virtual Python Environment builder to create an environment dedicated to **EXOTIC** and isolated from the native ***openSUSE*** installation of Python.

### 2.1 Set Up Virtual Environment
1. Ensure that **Python 3.7** or later is installed

   `python3 --version`

2. Note the installation location of **Python 3**

   `which python3`

3. Ensure that ***virtualenv*** and supporting tools are installed

   `sudo zypper install python3-virtualenvwrapper python38-virtualenv python38-virtualenv-clone`
   
   *Note: This command only performs the installation if the tool is not already installed (and may be ommited in this case)*

4. Create directory for Python virtual environments

   `mkdir ~/env`
   
   *Note: If a Python virtual environment directory already exists and is desired to be used, this step is omitted*

5. Edit ***.bashrc*** file

   `echo "export WORKON_HOME=~/env" >> ~/.bashrc`

   Run the modified script:
   
   `source ~/.bashrc`

6. Create Python virtual environment for **EXOTIC**

   `mkdir -p $WORKON_HOME` 

7. Initialize ***virtualenvwrapper***

   `virtualenvwrapper`
   
   *Note: This command creates various new directories within the newly created Python virtual environment directory. Optionally, the command may be added to the ***.bashrc*** file which executes each time a terminal window is launched:*
   
   `echo "source virtualenvwrapper" >> ~/.bashrc`

8. Create new Python virtual environment for **EXOTIC**

   `mkvirtualenv exotic`
   
   *Note: The virtual Python environment for ***EXOTIC*** is now created and can be confirmed by the "(exotic)" prefix in the terminal:*
   ```bash
   (exotic) user@localhost:~>
   ```
      
### 2.2 Install and Run EXOTIC

9. Ensure that the **EXOTIC** virtual environment is active (verified by the "*exotic*" prefix in the terminal window)

10. Issue the installation command

    `pip install exotic`

11. Run **EXOTIC**

    `exotic`
