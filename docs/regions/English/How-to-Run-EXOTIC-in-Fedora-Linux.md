# Install EXOTIC on Fedora Linux

EXOTIC requires Python 3.12 or newer. The commands here use Fedora's Python 3.13
package; the shared guide also covers Python 3.14 environments.
Fedora provides a [Python 3.13 package](https://packages.fedoraproject.org/pkgs/python3.13/python3.13/):

```bash
sudo dnf install python3.13
```

For EXOTIC's optional graphical interface, install the matching
[Tkinter package](https://packages.fedoraproject.org/pkgs/python3.13/python3.13-tkinter/):

```bash
sudo dnf install python3.13-tkinter
```

Continue with the shared [Linux installation instructions](EXOTIC-Installation-Instructions-for-Linux.md)
for environment creation, installation, and launching EXOTIC.
