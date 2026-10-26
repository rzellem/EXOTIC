# Install EXOTIC on Ubuntu

EXOTIC requires **Python 3.12 or newer**. The shared guide uses Python 3.14.
Use the [Linux installation guide](EXOTIC-Installation-Instructions-for-Linux.md)
to create an environment, install EXOTIC, and run its command-line or graphical
interface.

If your Ubuntu release supplies a supported default `python3`, its virtual
environment support can be installed with:

```bash
sudo apt install python3-venv
```

For the optional graphical interface, install `python3-tk`. These packages match
Ubuntu's default Python; if you install another Python version, use the matching
venv and Tk packages from that interpreter's provider instead.

Windows users can use the [native Windows instructions](<EXOTIC Installation Instructions for Windows 10.md>).
WSL also supports Linux GUI applications when configured as described in
[Microsoft's WSL GUI guide](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps).
The EXOTIC command-line interface does not require a graphical display.

See the [repository README](../../../README.md) for Google Colab tutorials,
sample data, and optional MPI support.
