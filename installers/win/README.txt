The following software packages should be installed in order: 

1-git-2.21.0-64-bit.zip
* Description: 
    Provides Git command line terminal with Unix-like tools.
* Download: 
    https://github.com/git-for-windows/git/releases/tag/v2.21.0.windows.1
* Notes: 
    Installation is complete when a Git Bash icon is available on the desktop.

2-python-3.7.3-amd64.zip
* Description: 
    Python programming interpreter.
* Download: 
    https://www.python.org/ftp/python/3.7.3/
* Notes: 
    Test install with `python3 --version`.

3-mingw-w64-install.zip
* Description: 
    Extended Unix-like tooling, including gcc compiler.
* Download: 
    https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download
    or 
    https://sourceforge.net/projects/mingw-w64/files/mingw-w64/mingw-w64-release/
* Notes: 
    Use a custom setup and install on top of Git Bash in the `<Git_install_directory>/mingw64` directory. This will make the latest tooling available in the Git command-line.

