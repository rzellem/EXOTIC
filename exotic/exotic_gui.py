# ########################################################################### #
# #    Copyright (c) 2019-2020, California Institute of Technology.
# #    All rights reserved.  Based on Government Sponsored Research under
# #    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
# #
# #    Redistribution and use in source and binary forms, with or without
# #    modification, are permitted provided that the following conditions
# #    are met:
# #      1. Redistributions of source code must retain the above copyright
# #         notice, this list of conditions and the following disclaimer.
# #      2. Redistributions in binary form must reproduce the above copyright
# #         notice, this list of conditions and the following disclaimer in
# #         the documentation and/or other materials provided with the
# #         distribution.
# #      3. Neither the name of the California Institute of
# #         Technology (Caltech), its operating division the Jet Propulsion
# #         Laboratory (JPL), the National Aeronautics and Space
# #         Administration (NASA), nor the names of its contributors may be
# #         used to endorse or promote products derived from this software
# #         without specific prior written permission.
# #
# #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# #    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# #    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# #    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
# #    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# #    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# #    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# #    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# #    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# #    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# #
# # ########################################################################### #
# #    EXOplanet Transit Interpretation Code (EXOTIC)
# #    # NOTE: See companion file version.py for version info.
# # ########################################################################### #

try:
    from animate import *
except ImportError:
    from .animate import *
animate_toggle(True)

import ast
from datetime import datetime
import json
import os
import platform
import python_version
import subprocess
import sys

try:
    import tkinter as tk
    from tkinter import filedialog
except ImportError:
    print(f"\nERROR: Incomplete Python 3 installation -- missing `GUI` tools. Please\n"
          "modify with support for 'python3-tk' or 'Tkinter', or reinstall. ... EXITING!\n")
    print("Press the <ENTER> key to continue. ...")
    sys.stdin.read(1)
    exit(79)  # cannot access a shared library


try:  # filters
    from .api.filters import fwhm as photometric_filters
except ImportError:  # package import
    from api.filters import fwhm as photometric_filters

try:
    from .api.nea import NASAExoplanetArchive
except ImportError:  # package import
    from api.nea import NASAExoplanetArchive

try:  # simple version
    from .version import __version__
except ImportError:  # package import
    from version import __version__

animate_toggle()


class FolderSelect(tk.Frame):
    def __init__(self, parent=None, folderDescription="", default_text="", **kw):
        tk.Frame.__init__(self, master=parent, **kw)
        self.folderPath = tk.StringVar()
        self.lblName = tk.Label(self, text=folderDescription, anchor=tk.W)
        self.lblName.grid(row=0, column=0)
        self.entPath = tk.Entry(self, textvariable=self.folderPath)
        self.entPath.insert(tk.END, default_text)
        self.entPath.grid(row=0, column=1)
        self.btnFind = tk.Button(self, text="Browse Folder", command=self.setFolderPath, anchor=tk.W)
        self.btnFind.grid(row=0, column=2)

    def setFolderPath(self):
        folder_selected = filedialog.askdirectory()
        self.folderPath.set(folder_selected)

    @property
    def folder_path(self):
        return self.folderPath.get()


class FileSelect(tk.Frame):
    def __init__(self, parent=None, folderDescription="", default_text="", **kw):
        tk.Frame.__init__(self, master=parent, **kw)
        self.filePath = tk.StringVar()
        self.lblName = tk.Label(self, text=folderDescription, anchor=tk.W)
        self.lblName.grid(row=0, column=0)
        self.entPath = tk.Entry(self, textvariable=self.filePath)
        self.entPath.insert(tk.END, default_text)
        self.entPath.grid(row=0, column=1)
        self.btnFind = tk.Button(self, text="Browse", command=self.setFilePath, anchor=tk.W)
        self.btnFind.grid(row=0, column=2)

    def setFilePath(self):
        file_selected = filedialog.askopenfilename()
        self.filePath.set(file_selected)

    @property
    def file_path(self):
        return self.filePath.get()


def main():
    try:
        python_version.check(min=(3, 8, 0), max=(4, 0, 0))
    except Exception as e:
        print(f"\nERROR: EXOTIC {str(e)}. ... EXITING!\n")
        print("Press the <ENTER> key to continue. ...")
        sys.stdin.read(1)
        exit(65)  # package is not installed
    else:
        print(f"\nSUCCESS: Valid Python version {platform.python_version()} detected!\n")

    root=tk.Tk() 
    root.protocol("WM_DELETE_WINDOW", exit)

    root.title(f"EXOTIC v{__version__}")

    tk.Label(root,
             text="Welcome to EXOTIC!",
             justify=tk.LEFT,
             font=("Arial", 36),
             padx=20).pack(anchor=tk.CENTER)
    tk.Label(root,
             text="the EXOplanet Transit Interpretation Code",
             justify=tk.LEFT,
             font=("Arial", 22),
             padx=20).pack(anchor=tk.CENTER)
    tk.Label(root,
             text=f"Version {__version__}",
             justify=tk.LEFT,
             padx=20).pack(anchor=tk.CENTER)

    tk.Label(root,
             text="Copyright (c) 2021, California Institute of Technology. "
                  "All rights reserved.\nBased on Government Sponsored Research "
                  "under contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.",
             font=("Arial", 12),
             justify=tk.CENTER,
             padx=20).pack(anchor=tk.CENTER)

    # Button for closing
    exit_button = tk.Button(root, text="Next", command=root.destroy)
    exit_button.pack(anchor=tk.CENTER)

    root.mainloop()

    root=tk.Tk() 
    root.protocol("WM_DELETE_WINDOW", exit)
    root.title(f"EXOTIC v{__version__}")

    reduction_opt = tk.IntVar()
    tk.Label(root,
             text="\nHow would you like to run EXOTIC?",
             justify=tk.LEFT,
             font="Helvetica 14 bold",
             padx=20).pack(anchor=tk.W)

    tk.Radiobutton(root,
                   text="Real Time Reduction\n(for quickly analyzing your data while simultaneously observing)",
                   justify=tk.LEFT,
                   padx=20,
                   variable=reduction_opt,
                   value=1).pack(anchor=tk.W)

    tk.Radiobutton(root,
                   text="Complete Reduction\n(for analyzing your data after an observing run)",
                   justify=tk.LEFT,
                   padx=20,
                   variable=reduction_opt,
                   value=2).pack(anchor=tk.W)
    reduction_opt.set(2)

    # Button for closing
    exit_button = tk.Button(root, text="Next", command=root.destroy)
    exit_button.pack(anchor=tk.E)

    root.mainloop()

    if reduction_opt.get() == 1:
        # First ask user how they want to enter the observing information
        root=tk.Tk() 
        root.protocol("WM_DELETE_WINDOW", exit)
        obsinfo = tk.StringVar()

        input_data = {}

        root.title(f"EXOTIC v{__version__}")

        tk.Label(root,
                 text="""How would you like to input your observing information?""",
                 justify=tk.LEFT,
                 font="Helvetica 14 bold",
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually (Recommended for first-time users)",
                       padx=20,
                       variable=obsinfo,
                       value="manual").pack(anchor=tk.W)
        obsinfo.set("manual")

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g. inits.json) - for advanced users",
                       padx=20,
                       variable=obsinfo,
                       value='inits').pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20, anchor=tk.E)

        root.mainloop()

        if obsinfo.get() == "manual":
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            window_label = tk.Label(root,
                                    text="""Please enter the following information about your observation:""",
                                    font="Helvetica 14 bold",
                                    justify=tk.CENTER,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1
            j = 0

            # "Directory with FITS files": "sample-data/HatP32Dec202017",
            FITS_dir = FolderSelect(root, "Folder with your FITS files")
            FITS_dir.grid(row=i)
            i += 1

            #         "Planet Name": "HAT-P-32 b",
            planet_label = tk.Label(root, text="Planet Name", justify=tk.LEFT)
            planet_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            planet_entry.insert(tk.END, "HAT-P-32 b")
            planet_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            planet_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            targetpos_label = tk.Label(root, text="Target Star X & Y Pixel Position", justify=tk.LEFT)
            targetpos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            targetpos_entry.insert(tk.END, "[x, y]")
            targetpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            comppos_label = tk.Label(root, text="Comparison Star(s) X & Y Pixel Position(s)\n    "
                                                "(Note: You can use the AAVSO's VSP to help you find\n    "
                                                "good comparison stars: https://apps.aavso.org/vsp/)",
                                     justify=tk.LEFT)
            comppos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            comppos_entry.insert(tk.END, "[x, y]")
            comppos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comppos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)

            def save_input():
                input_data['comppos'] = ast.literal_eval(comppos_entry.get())
                input_data['targetpos'] = ast.literal_eval(targetpos_entry.get())
                input_data['pName'] = planet_entry.get()
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()
        else:
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            inits_dir = FileSelect(root, "Please select your initialization file")
            inits_dir.grid(row=0)

            def save_input():
                input_data['inits_dir'] = inits_dir.file_path
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Run EXOTIC", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        if obsinfo.get() == "manual":
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            window_label = tk.Label(root,
                                    text="To make analyzing data easier in the future, EXOTIC will create an initialization file for you.",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            initssave_dir = FolderSelect(root, "Folder to save your initialization file:")
            initssave_dir.grid(row=i)
            i += 1

            def save_input():
                input_data['initssave_dir'] = initssave_dir.folderPath.get()
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

            new_inits = {'inits_guide': {}, 'user_info': {}, 'planetary_parameters': {}, 'optional_info': {}}
            new_inits['inits_guide'] = {
                "Title": "EXOTIC's Initialization File",
                "Comment": "Please answer all the following requirements below by following the format of the given",
                "Comment1": "sample dataset HAT-P-32 b. Edit this file as needed to match the data wanting to be reduced.",
                "Comment2": "Do not delete areas where there are quotation marks, commas, and brackets.",
                "Comment3": "The inits_guide dictionary (these lines of text) does not have to be edited",
                "Comment4": "and is only here to serve as a guide. Will be updated per user's advice.",
                "Image Calibrations Directory Guide": "Enter in the path to image calibrations or enter in null for none.",
                "Planetary Parameters Guide": "For planetary parameters that are not filled in, enter in null.",
                "Comparison Star(s) Guide": "Up to 10 comparison stars can be added following the format given below.",
                "Obs. Latitude Guide": "Indicate the sign (+ North, - South) before the degrees. Needs to be in decimal or HH:MM:SS format.",
                "Obs. Longitude Guide": "Indicate the sign (+ East, - West) before the degrees. Needs to be in decimal or HH:MM:SS format.",
                "Plate Solution": "For your image to be given a plate solution, type y.",
                "Plate Solution Disclaimer": "One of your imaging files will be publicly viewable on nova.astrometry.net.",
                "Standard Filter": "To use EXOTIC standard filters, type only the filter name.",
                "Custom Filter": "To use a custom filter, enter in the FWHM in optional_info.",
                "Target Star RA": "Must be in HH:MM:SS sexagesimal format.",
                "Target Star DEC": "Must be in +/-DD:MM:SS sexagesimal format with correct sign at the beginning (+ or -).",
                "Demosaic Format": "Optional control for handling Bayer pattern color images - to use, provide Bayer color patttern of your camera (RGGB, BGGR, GRBG, GBRG) - null (no color processing) is default",
                "Demosaic Output": "Select how to process color data (gray for grayscale, red or green or blue for single color channel, blueblock for grayscale without blue, [ R, G, B ] for custom weights for mixing colors.  green is default",
                "Formatting of null": "Due to the file being a .json, null is case sensitive and must be spelled as shown.",
                "Decimal Format": "Leading zero must be included when appropriate (Ex: 0.32, .32 or 00.32 causes errors.)."
            }
            new_inits['user_info'] = {
                "Directory with FITS files": FITS_dir.folder_path,

                "Target Star X & Y Pixel": input_data['targetpos'],
                "Comparison Star(s) X & Y Pixel": [input_data['comppos']],
                "Demosaic Format": null, # TODO add GUI input for these
                "Demosaic Output": null
            }
            new_inits['planetary_parameters'] = {
                "Planet Name": input_data['pName'],
            }

            now = datetime.now()
            dt_string = now.strftime("%d_%m_%Y__%H_%M_%S")
            fname = os.path.join(input_data['initssave_dir'], "inits_" + dt_string + ".json")
            with open(fname, "w") as initsf:
                json.dump(new_inits, initsf, indent=4)
            print(f"{fname} saved!")

            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            window_label = tk.Label(root,
                                    text="Your initialization file has been created!",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            window_label = tk.Label(root,
                                    text=f"{fname}",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=1, column=0, sticky=tk.N, pady=6)

            # Button for closing
            exit_button = tk.Button(root, text="Run EXOTIC", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=2, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        run_mthd = '--realtime'

        fname = fname if obsinfo.get() == "manual" else inits_dir.file_path

        try:
            try:
                subprocess.run(['exotic', run_mthd, fname], check=True)
            except:
                try:
                    subprocess.run(['python3', 'exotic.py', run_mthd, fname], check=True)
                except:
                    subprocess.run(['python3', 'exotic/exotic.py', run_mthd, fname], check=True)
        except:
            print("\n\n################################################")
            print("ERROR: Please contact the Exoplanet Watch Team for help "
                  "on our Slack Workspace in the #data-reduction channel!")
            print("You can sign up for our free Slack Workspace here: "
                  "https://join.slack.com/t/uol-ets/shared_invite/zt-mvb4ljbo-LRBgpk3uMmUokbs4ge2JlA")
            print("################################################\n\n")
            pass
    else:
        root=tk.Tk() 
        root.protocol("WM_DELETE_WINDOW", exit)
        root.title(f"EXOTIC v{__version__}")

        fitsortext = tk.IntVar()

        tk.Label(root,
                 text="""How do you want to run EXOTIC?""",
                 font="Helvetica 14 bold",
                 justify=tk.LEFT,
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Analyze your image files",
                       padx=20,
                       variable=fitsortext,
                       value=1).pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="Start with pre-reduced data in a .txt format",
                       padx=20,
                       variable=fitsortext,
                       value=2).pack(anchor=tk.W)
        fitsortext.set(2)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20, anchor=tk.E)
        input_data = {}  # for saving entries

        root.mainloop()

        if fitsortext.get() == 2:
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")
            i=0; j=0

            prered_file = FileSelect(root, "Pre-reduced Data File")
            prered_file.grid(row=i)
            i += 1

            #            # "Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)": "BJD_TDB",
            pretime_label = tk.Label(root, text="Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)",
                                     justify=tk.LEFT)
            pretime_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            pretime_entry.insert(tk.END, "")
            pretime_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            pretime_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #            # "Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)": "flux",
            preunit_label = tk.Label(root, text="Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)",
                                     justify=tk.LEFT)
            preunit_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            preunit_entry.insert(tk.END, "")
            preunit_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            preunit_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #              "Exposure Time (s)": 60.0
            exp_label = tk.Label(root, text="Exposure Time (s)", justify=tk.LEFT)
            exp_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            exp_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            exp_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            comp_star_label = tk.Label(root, text="Comparison Star used in Photometry (leave blank if none):",
                                       justify=tk.LEFT)
            comp_star_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            i += 1

            comp_star_ra_label = tk.Label(root, text="Comparison Star RA", justify=tk.LEFT)
            comp_star_ra_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            comp_star_ra_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comp_star_ra_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            comp_star_dec_label = tk.Label(root, text="Comparison Star DEC", justify=tk.LEFT)
            comp_star_dec_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            comp_star_dec_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comp_star_dec_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            comp_star_x_label = tk.Label(root, text="Comparison Star X Pixel Coordinate", justify=tk.LEFT)
            comp_star_x_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            comp_star_x_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comp_star_x_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            comp_star_y_label = tk.Label(root, text="Comparison Star Y Pixel Coordinate", justify=tk.LEFT)
            comp_star_y_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            comp_star_y_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comp_star_y_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            def save_input():
                input_data['file_time'] = pretime_entry.get()
                input_data['file_units'] = preunit_entry.get()
                input_data['exp'] = float(exp_entry.get())
                input_data['phot_comp_star'] = {'ra': comp_star_ra_entry.get(), 'dec': comp_star_dec_entry.get(),
                                                'x': comp_star_x_entry.get(), 'y': comp_star_y_entry.get()}
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            root.mainloop()

        # First ask user how they want to enter the observing information
        root=tk.Tk() 
        root.protocol("WM_DELETE_WINDOW", exit)
        obsinfo = tk.StringVar()

        root.title(f"EXOTIC v{__version__}")

        tk.Label(root,
                 text="""How would you like to input your observing information?""",
                 justify=tk.LEFT,
                 font="Helvetica 14 bold",
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually (Recommended for first-time users)",
                       padx=20,
                       variable=obsinfo,
                       value="manual").pack(anchor=tk.W)
        obsinfo.set("manual")

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g. inits.json) - for advanced users",
                       padx=20,
                       variable=obsinfo,
                       value='inits').pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20, anchor=tk.E)

        root.mainloop()


        if obsinfo.get() == "manual":
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="""Please enter the following information about your observation:""",
                                    font="Helvetica 14 bold",
                                    justify=tk.CENTER,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1
            j = 0

            folderPath = tk.StringVar()

            # "Directory with FITS files": "sample-data/HatP32Dec202017",
            if fitsortext.get() == 1:
                FITS_dir = FolderSelect(root, "Folder with your FITS files")
                FITS_dir.grid(row=i)
                i += 1

            # "Directory to Save Plots": "sample-data/",
            save_dir = FolderSelect(root, "Folder to save EXOTIC output")
            save_dir.grid(row=i)
            i += 1

            window_label = tk.Label(root,
                                    text="""Please enter your calibration file information
                                    \n(Note: each calibration type must be in its OWN, SEPARATE folder)""",
                                    # font=("Helvetica 14 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=i, column=0, sticky=tk.N, pady=6)
            i += 1

            # "Directory of Flats": null,
            if fitsortext.get() == 1:
                flats_dir = FolderSelect(root, "Folder with your flats (null if none)", "null")
                flats_dir.grid(row=i)
                i += 1

            # "Directory of Darks": null,
            if fitsortext.get() == 1:
                darks_dir = FolderSelect(root, "Folder with your darks (null if none)", "null")
                darks_dir.grid(row=i)
                i += 1

            # "Directory of Biases": null,
            if fitsortext.get() == 1:
                biases_dir = FolderSelect(root, "Folder with your biases (null if none)", "null")
                biases_dir.grid(row=i)
                i += 1

            # "Directory of Biases": null,
            # biases_dir = FolderSelect(root, "Folder with your biases (null if none)", "null")
            # biases_dir.grid(row=i)
            # i += 1

            #
            # obscode2_label = TextInput(root, "AAVSO Obs Code", "RTZ")
            # obscode2_label.grid(row=1)
            # i += 1

            # if fitsortext.get() == 2:
            #     prered_file = FileSelect(root, "Pre-reduced Data File Path")
            #     prered_file.grid(row=i)
            #     i += 1
            #
            #     #            # "Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)": "BJD_TDB",
            #     pretime_label = tk.Label(root, text="Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)", justify=tk.LEFT)
            #     pretime_entry  = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            #     pretime_entry.insert(tk.END, "")
            #     pretime_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     pretime_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            #     i += 1
            #
            #     #            # "Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)": "flux",
            #     preunit_label = tk.Label(root, text="Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)", justify=tk.LEFT)
            #     preunit_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            #     preunit_entry.insert(tk.END, "")
            #     preunit_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     preunit_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            #     i += 1

            #             # "AAVSO Observer Code (blank if none)": "RTZ",
            obscode_label = tk.Label(root, text="AAVSO Observer Code (leave blank if none)", justify=tk.LEFT)
            obscode_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            obscode_entry.insert(tk.END, "")
            obscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            #
            #             # "Secondary Observer Codes (blank if none)": "",
            secondobscode_label = tk.Label(root, text="Secondary Observer Codes (leave blank if none)", justify=tk.LEFT)
            secondobscode_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            secondobscode_entry.insert(tk.END, "")
            secondobscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            secondobscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            #             # "Observation date": "17-December-2017",
            obsdate_label = tk.Label(root, text="Observation date (e.g. DAY-MONTH-YEAR)", justify=tk.LEFT)
            obsdate_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            obsdate_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obsdate_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            # #             "Obs. Latitude": "+32.41638889",
            lat_label = tk.Label(root, text="Obs. Latitude (+ = North; - = South; e.g. +32.41)", justify=tk.LEFT)
            lat_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            lat_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            lat_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            # #             "Obs. Longitude": "-110.73444444",
            long_label = tk.Label(root, text="Obs. Longitude (+ = East; - = West; e.g. -110.74) ", justify=tk.LEFT)
            long_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            long_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            long_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            # #             "Obs. Elevation (meters)": 2616,
            elevation_label = tk.Label(root, text="Obs. Elevation [meters]", justify=tk.LEFT)
            elevation_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            elevation_entry.insert(tk.END, "0")
            elevation_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            elevation_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            # #             "Camera Type (CCD or DSLR)": "CCD",
            cameratype_label = tk.Label(root,
                                        text="Camera Type (e.g., CCD or DSLR;\n"
                                             "Note: if you are using a CMOS, please enter CCD here and\n"
                                             "then note your actual camera type under \"Observing Notes\" below)",
                                        justify=tk.LEFT)
            cameratype_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            cameratype_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            cameratype_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
            #
            # #             "Pixel Binning": "1x1",
            pixbin_label = tk.Label(root, text="Pixel Binning  (e.g 1x1)", justify=tk.LEFT)
            pixbin_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            pixbin_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            pixbin_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #     "Pixel Scale (Ex: 5.21 arcsecs/pixel)": null,
            pixscale_label = tk.Label(root, text="Pixel Scale (e.g., 5.21 arcsecs/pixel; null if unknown)", justify=tk.LEFT)
            pixscale_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            # pixscale_entry.insert(tk.END, "")
            pixscale_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            pixscale_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            # # #             "Filter Name (aavso.org/filters)": "V",
            # filtername_label = tk.Label(root, text="Filter Name (please see aavso.org/filters)", justify=tk.LEFT)
            # filtername_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            # filtername_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            # filtername_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            # i += 1

            # choices = ['one', 'two', 'three']
            choices = [item for item in photometric_filters.keys()] + ["N/A"]
            choices = sorted(set(choices))  # sort and list unique values
            filteroptions = tk.StringVar(root)
            filteroptions.set(choices[0])  # default value

            l3 = tk.Label(root, text='Filter (use N/A for custom)', justify=tk.LEFT)
            l3.grid(row=i, column=j, sticky=tk.W, pady=2)

            om1 = tk.OptionMenu(root, filteroptions, *choices)
            om1.grid(row=i, column=1)
            # my_w.mainloop()  # Keep the window open
            i += 1

            # if fitsortext.get() == 2:
            #     #              "Exposure Time (s)": 60.0
            #     exp_label = tk.Label(root, text="Exposure Time (s)", justify=tk.LEFT)
            #     exp_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            #     exp_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     exp_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            #     i += 1

            #              "Observing Notes": "Weather, seeing was nice.",
            obsnotes_label = tk.Label(root, text="Observing Notes", justify=tk.LEFT)
            obsnotes_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            obsnotes_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obsnotes_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            # root=tk.Tk() 
            # root.protocol("WM_DELETE_WINDOW", exit)
            # root.title(f"EXOTIC v{__version__}")
            #
            # window_label = tk.Label(root,
            #                         text="""Please enter the following information:""",
            #                         font=("Helvetica 15 bold"),
            #                         justify=tk.LEFT,
            #                         padx=20)  # .pack()
            # window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # i, j = 1, 0
            if fitsortext.get() == 1:
                platesolve = tk.BooleanVar()
                platesolve_check = tk.Checkbutton(root, text='Plate solve my first image'
                                                             '\n(Note: Do not select if all images are'
                                                             '\nplate solved and you trust it)',
                                                  variable=platesolve, onvalue=True, offvalue=False, justify=tk.LEFT)
                platesolve_check.grid(row=i, column=j, sticky=tk.W, pady=2)
                i += 1

                aavsocomp = tk.BooleanVar()
                aavsocomp_check = tk.Checkbutton(root, text="Add Comparison Stars Automatically"
                                                            "\nfrom AAVSO's VSP API for Stellar"
                                                            "\nVariability Reduction \n"
                                                            "(must have plate solution)",
                                                  variable=aavsocomp, onvalue=True, offvalue=False, justify=tk.LEFT)
                aavsocomp_check.grid(row=i, column=j, sticky=tk.W, pady=2)
                i += 1

                targetpos_label = tk.Label(root, text="Target Star X & Y Pixel Position", justify=tk.LEFT)
                targetpos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
                targetpos_entry.insert(tk.END, "[x, y]")
                targetpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                targetpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                i += 1

                comppos_label = tk.Label(root, text="Comparison Star(s) X & Y Pixel Position(s)\n    "
                                                    "(Note: You can use the AAVSO's VSP to help you find\n    "
                                                    "good comparison stars: https://apps.aavso.org/vsp/)",
                                         justify=tk.LEFT)
                comppos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
                comppos_entry.insert(tk.END, "[x1, y1], [x2, y2]")
                comppos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                comppos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                i += 1

                # # Button for closing
                # exit_button = tk.Button(root, text="Next", command=root.destroy)
                # exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)
                # # exit_button.pack(pady=20)
                # root.update()
                # root.mainloop()

            def save_input():
                input_data['obsnotes'] = obsnotes_entry.get()
                if filteroptions.get() == "N/A":
                    input_data['obsfilter'] = "N/A"
                else:
                    input_data['obsfilter'] = photometric_filters[filteroptions.get()]["name"]
                input_data['pixbin'] = pixbin_entry.get()
                input_data['cameratype'] = cameratype_entry.get()
                input_data['obscode'] = obscode_entry.get()
                input_data['secondobscode'] = secondobscode_entry.get()
                input_data['obsdate'] = obsdate_entry.get()
                input_data['lat'] = lat_entry.get()
                input_data['long'] = long_entry.get()
                input_data['elevation'] = float(elevation_entry.get())
                input_data['pixscale'] = pixscale_entry.get()
                if fitsortext.get() == 1:
                    input_data['comppos'] = str(list(ast.literal_eval(comppos_entry.get())))
                    input_data['targetpos'] = str(list(ast.literal_eval(targetpos_entry.get())))

                    if platesolve.get():
                        input_data['platesolve'] = 'y'
                    else:
                        input_data['platesolve'] = 'n'

                    if aavsocomp.get():
                        input_data['aavso_comp'] = 'y'
                    else:
                        input_data['aavso_comp'] = 'n'
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()
        else:
            # raise(Exception("feature not supported yet"))
            pass

        try:
            if filteroptions.get() == "N/A":
                root=tk.Tk() 
                root.protocol("WM_DELETE_WINDOW", exit)
                root.title(f"EXOTIC v{__version__}")

                i = 0
                j = 0

                tk.Label(root,
                         text="Please enter the following information for your custom filter:",
                         font="Helvetica 15 bold",
                         justify=tk.LEFT,
                         padx=20).grid(row=i, column=j, sticky=tk.W, pady=2)
                i += 1

                # #
                # "Filter Minimum Wavelength (nm)": null,
                filtermin_label = tk.Label(root, text="Filter Minimum Wavelength (nm)", justify=tk.LEFT)
                filtermin_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
                filtermin_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                filtermin_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                i += 1

                # "Filter Maximum Wavelength (nm)": null
                filtermax_label = tk.Label(root, text="Filter Maximum Wavelength (nm)", justify=tk.LEFT)
                filtermax_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
                filtermax_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                filtermax_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                i += 1

                def save_input():
                    input_data['filtermax'] = float(filtermax_entry.get())
                    input_data['filtermin'] = float(filtermin_entry.get())
                    root.destroy()

                # Button for closing
                exit_button = tk.Button(root, text="Next", command=save_input)
                # exit_button.pack(pady=20)
                exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

                root.mainloop()
        except:
            pass

        # then ask user how they want to enter the planetary information
        root=tk.Tk() 
        root.protocol("WM_DELETE_WINDOW", exit)
        root.title(f"EXOTIC v{__version__}")

        planetparams = tk.StringVar()

        tk.Label(root,
                 text="""How would you like to input your planetary system parameters?""",
                 font="Helvetica 15 bold",
                 justify=tk.LEFT,
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually (Recommended for first-time users)",
                       padx=20,
                       variable=planetparams,
                       value="manual").pack(anchor=tk.W)
        planetparams.set("manual")

        tk.Radiobutton(root,
                       text="Automatically adopt all planetary parameters\nfrom the NASA Exoplanet Archive",
                       padx=20,
                       variable=planetparams,
                       value="nea").pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g., inits.json) - for advanced users",
                       padx=20,
                       variable=planetparams,
                       value="inits").pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20, anchor=tk.E)

        root.mainloop()

        if planetparams.get() == "manual":
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="""Please enter the following information:""",
                                    font="Helvetica 15 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1
            j = 0

            #         "Planet Name": "HAT-P-32 b",
            planet_label = tk.Label(root, text="Planet Name", justify=tk.LEFT)
            planet_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            planet_entry.insert(tk.END, "HAT-P-32 b")
            planet_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            planet_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Host Star Name": "HAT-P-32",
            star_label = tk.Label(root, text="Host Star Name", justify=tk.LEFT)
            star_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            star_entry.insert(tk.END, "HAT-P-32")
            star_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            star_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Target Star RA": "02:04:10",
            targetRA_label = tk.Label(root, text="Host Star Right Ascension", justify=tk.LEFT)
            targetRA_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            targetRA_entry.insert(tk.END, "02:04:10")
            targetRA_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetRA_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Target Star Dec": "+46:41:23",
            targetDEC_label = tk.Label(root, text="Host Star Declination", justify=tk.LEFT)
            targetDEC_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            targetDEC_entry.insert(tk.END, "+46:41:23")
            targetDEC_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetDEC_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Orbital Period (days)": 2.1500082,
            period_label = tk.Label(root, text="Planet's orbital period (days)", justify=tk.LEFT)
            period_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            period_entry.insert(tk.END, "2.1500082")
            period_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            period_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Orbital Period Uncertainty": 1.3e-07,
            perioderr_label = tk.Label(root, text="Planet's orbital period uncertainty (days)", justify=tk.LEFT)
            perioderr_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            perioderr_entry.insert(tk.END, "1.3e-07")
            perioderr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            perioderr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Published Mid-Transit Time (BJD-UTC)": 2455867.402743,
            Tmid_label = tk.Label(root, text="Published Mid-Transit Time (BJD-UTC)", justify=tk.LEFT)
            Tmid_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            Tmid_entry.insert(tk.END, "2455867.402743")
            Tmid_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tmid_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Mid-Transit Time Uncertainty": 4.9e-05,
            Tmiderr_label = tk.Label(root, text="Mid-Transit Time Uncertainty", justify=tk.LEFT)
            Tmiderr_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            Tmiderr_entry.insert(tk.END, "4.9e-05")
            Tmiderr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tmiderr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.14886235252742716,
            rprs_label = tk.Label(root, text="Ratio of Planet to Stellar Radius (Rp/Rs)", justify=tk.LEFT)
            rprs_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            rprs_entry.insert(tk.END, "0.14886235252742716")
            rprs_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            rprs_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.0005539487393037134,
            rprserr_label = tk.Label(root, text="Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty",
                                     justify=tk.LEFT)
            rprserr_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            rprserr_entry.insert(tk.END, "0.0005539487393037134")
            rprserr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            rprserr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Ratio of Distance to Stellar Radius (a/Rs)": 5.344,
            aRs_label = tk.Label(root, text="Ratio of Distance to Stellar Radius (a/Rs)", justify=tk.LEFT)
            aRs_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            aRs_entry.insert(tk.END, "5.344")
            aRs_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            aRs_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.039496835316262996,
            aRserr_label = tk.Label(root, text="Ratio of Distance to Stellar Radius (a/Rs) Uncertainty",
                                    justify=tk.LEFT)
            aRserr_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            aRserr_entry.insert(tk.END, "0.039496835316262996")
            aRserr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            aRserr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Orbital Inclination (deg)": 88.98,
            inc_label = tk.Label(root, text="Orbital Inclination (degrees)", justify=tk.LEFT)
            inc_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            inc_entry.insert(tk.END, "88.98")
            inc_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            inc_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Orbital Inclination (deg) Uncertainty": 0.7602631123499285,
            incerr_label = tk.Label(root, text="Orbital Inclination (degrees) Uncertainty", justify=tk.LEFT)
            incerr_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            incerr_entry.insert(tk.END, "0.7602631123499285")
            incerr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            incerr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Orbital Eccentricity (0 if null)": 0.159,
            ecc_label = tk.Label(root, text="Orbital Eccentricity (0 if null)", justify=tk.LEFT)
            ecc_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            ecc_entry.insert(tk.END, "0.159")
            ecc_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            ecc_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            # "Argument of Periastron (deg)": 50,
            omega_label = tk.Label(root, text="Argument of Periastron (degrees)", justify=tk.LEFT)
            omega_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            omega_entry.insert(tk.END, "0.")
            omega_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            omega_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Effective Temperature (K)": 6001.0,
            Teff_label = tk.Label(root, text="Star Effective Temperature (K)", justify=tk.LEFT)
            Teff_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            Teff_entry.insert(tk.END, "6001.0")
            Teff_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Teff_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Effective Temperature (+) Uncertainty": 88.0,
            Tefferrpos_label = tk.Label(root, text="Star Effective Temperature Positive (+) Uncertainty",
                                        justify=tk.LEFT)
            Tefferrpos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            Tefferrpos_entry.insert(tk.END, "+88.0")
            Tefferrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tefferrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Effective Temperature (-) Uncertainty": -88.0,
            Tefferrneg_label = tk.Label(root, text="Star Effective Temperature Negative (-) Uncertainty",
                                        justify=tk.LEFT)
            Tefferrneg_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            Tefferrneg_entry.insert(tk.END, "-88.0")
            Tefferrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tefferrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Metallicity ([FE/H])": -0.16,
            FeH_label = tk.Label(root, text="Star Metallicity ([Fe/H])", justify=tk.LEFT)
            FeH_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            FeH_entry.insert(tk.END, "-0.16")
            FeH_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeH_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Metallicity (+) Uncertainty": 0.08,
            FeHerrpos_label = tk.Label(root, text="Star Metallicity Positive (+) Uncertainty", justify=tk.LEFT)
            FeHerrpos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            FeHerrpos_entry.insert(tk.END, "0.08")
            FeHerrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeHerrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Metallicity (-) Uncertainty": -0.08,
            FeHerrneg_label = tk.Label(root, text="Star Metallicity Negative (-) Uncertainty", justify=tk.LEFT)
            FeHerrneg_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            FeHerrneg_entry.insert(tk.END, "-0.08")
            FeHerrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeHerrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Surface Gravity (log(g))": 4.22,
            logg_label = tk.Label(root, text="Star Surface Gravity (log(g))", justify=tk.LEFT)
            logg_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            logg_entry.insert(tk.END, "4.22")
            logg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            logg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Surface Gravity (+) Uncertainty": 0.04,
            loggerrpos_label = tk.Label(root, text="Star Surface Gravity Positive (+) Uncertainty", justify=tk.LEFT)
            loggerrpos_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            loggerrpos_entry.insert(tk.END, "0.04")
            loggerrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            loggerrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Star Surface Gravity (-) Uncertainty": -0.04
            loggerrneg_label = tk.Label(root, text="Star Surface Gravity Negative (-) Uncertainty", justify=tk.LEFT)
            loggerrneg_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            loggerrneg_entry.insert(tk.END, "-0.04")
            loggerrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            loggerrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            def save_input():
                input_data['ra'] = targetRA_entry.get()
                input_data['dec'] = targetDEC_entry.get()
                input_data['pName'] = planet_entry.get()
                input_data['sName'] = star_entry.get()
                input_data['pPer'] = period_entry.get()
                input_data['pPerUnc'] = perioderr_entry.get()
                input_data['midT'] = Tmid_entry.get()
                input_data['midTUnc'] = Tmiderr_entry.get()
                input_data['rprs'] = rprs_entry.get()
                input_data['rprsUnc'] = rprserr_entry.get()
                input_data['aRs'] = aRs_entry.get()
                input_data['aRsUnc'] = aRserr_entry.get()
                input_data['inc'] = inc_entry.get()
                input_data['incUnc'] = incerr_entry.get()
                input_data['ecc'] = ecc_entry.get()
                input_data['omega'] = omega_entry.get()
                input_data['teff'] = Teff_entry.get()
                input_data['teffUncPos'] = Tefferrpos_entry.get()
                input_data['teffUncNeg'] = Tefferrneg_entry.get()
                input_data['met'] = FeH_entry.get()
                input_data['metUncPos'] = FeHerrpos_entry.get()
                input_data['metUncNeg'] = FeHerrneg_entry.get()
                input_data['logg'] = logg_entry.get()
                input_data['loggUncPos'] = loggerrpos_entry.get()
                input_data['loggUncNeg'] = loggerrneg_entry.get()
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        elif planetparams.get() == 'nea':
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="""Please enter the following information:""",
                                    font="Helvetica 15 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1
            j = 0

            #         "Planet Name": "HAT-P-32 b",
            planet_label = tk.Label(root, text="Planet Name", justify=tk.LEFT)
            planet_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            planet_entry.insert(tk.END, "HAT-P-32 b")
            planet_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            planet_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            #         "Host Star Name": "HAT-P-32",
            star_label = tk.Label(root, text="Host Star Name", justify=tk.LEFT)
            star_entry = tk.Entry(root, font="Helvetica 12", justify=tk.LEFT)
            star_entry.insert(tk.END, "HAT-P-32")
            star_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            star_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1

            def save_input():
                input_data['sName'] = star_entry.get()
                input_data['pName'] = planet_entry.get()
                nea_obj = NASAExoplanetArchive(planet=input_data['pName'])
                if not os.path.exists('pl_names.json') or time.time() - os.path.getmtime('pl_names.json') > 2592000:
                    nea_obj.planet_names(filename="pl_names.json")
                if os.path.exists('pl_names.json'):
                    with open("pl_names.json", "r") as f:
                        planets = json.load(f)
                        for key, value in planets.items():
                            if input_data['pName'].lower().replace(' ', '').replace('-', '') == key:
                                input_data['pName'] = value
                                break
                input_data['pName'], CandidatePlanetBool, pDict = nea_obj.planet_info()
                for key in pDict:
                    input_data[key] = pDict[key]
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        if (planetparams.get() == "inits") or (obsinfo.get() == "inits"):
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            inits_dir = FileSelect(root, "Please select your initialization file")
            inits_dir.grid(row=0)

            def save_inputs():
                input_data['inits_dir'] = inits_dir.file_path
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_inputs)
            # exit_button.pack(pady=20)
            exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

            # if (obsinfo.get() == "manual"):
            #     root=tk.Tk() 
            # root.protocol("WM_DELETE_WINDOW", exit)
            #     root.title(f"EXOTIC v{__version__}")
            #
            #     window_label = tk.Label(root,
            #                             text="""Please enter the following information:""",
            #                             font=("Helvetica 15 bold"),
            #                             justify=tk.LEFT,
            #                             padx=20)  # .pack()
            #     window_label.grid(row=0, column=0, sticky=tk.N, pady=6)
            #
            #     i, j = 1, 0
            #     platesolve = tk.BooleanVar()
            #     platesolve_check = tk.Checkbutton(root, text='Plate solve my images', variable=platesolve, onvalue=True,
            #                                       offvalue=False)
            #     platesolve_check.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     i += 1
            #
            #     targetpos_label = tk.Label(root, text="Target Star X & Y Pixel Position", justify=tk.LEFT)
            #     targetpos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            #     targetpos_entry.insert(tk.END, "[x, y]")
            #     targetpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     targetpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            #     targetpos = targetpos_entry.get()
            #     i += 1
            #
            #     comppos_label = tk.Label(root, text="Comparison Star(s) X & Y Pixel Position", justify=tk.LEFT)
            #     comppos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            #     comppos_entry.insert(tk.END, "[x1, y1], [x2, y2]")
            #     comppos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            #     comppos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            #     comppos = comppos_entry.get()
            #     i += 1
            #
            #     # Button for closing
            #     exit_button = tk.Button(root, text="Next", command=root.destroy)
            #     exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)
            #     # exit_button.pack(pady=20)
            #     root.update()
            #     root.mainloop()
            # root.mainloop()

        # Create an initialization file if it does not already exist
        if (planetparams.get() != "inits") or (obsinfo.get() != "inits"):
            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="To make analyzing data easier in the future, EXOTIC will create an initialization file for you.",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1
            j = 0

            folderPath = tk.StringVar()

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            initssave_dir = FolderSelect(root, "Folder to save your initialization file:")
            initssave_dir.grid(row=i)
            i += 1

            def save_input():
                input_data['initssave_dir'] = initssave_dir.folderPath.get()
                root.destroy()

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=save_input)
            # exit_button.pack(pady=20)
            exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

            # create the inits file here
            new_inits = {'inits_guide': {}, 'user_info': {}, 'planetary_parameters': {}, 'optional_info': {}}
            new_inits['inits_guide'] = {
                "Title": "EXOTIC's Initialization File",
                "Comment": "Please answer all the following requirements below by following the format of the given",
                "Comment1": "sample dataset HAT-P-32 b. Edit this file as needed to match the data wanting to be reduced.",
                "Comment2": "Do not delete areas where there are quotation marks, commas, and brackets.",
                "Comment3": "The inits_guide dictionary (these lines of text) does not have to be edited",
                "Comment4": "and is only here to serve as a guide. Will be updated per user's advice.",
                "Image Calibrations Directory Guide": "Enter in the path to image calibrations or enter in null for none.",
                "Planetary Parameters Guide": "For planetary parameters that are not filled in, enter in null.",
                "Comparison Star(s) Guide": "Up to 10 comparison stars can be added following the format given below.",
                "Obs. Latitude Guide": "Indicate the sign (+ North, - South) before the degrees. Needs to be in decimal or HH:MM:SS format.",
                "Obs. Longitude Guide": "Indicate the sign (+ East, - West) before the degrees. Needs to be in decimal or HH:MM:SS format.",
                "Plate Solution": "For your image to be given a plate solution, type y.",
                "Plate Solution Disclaimer": "One of your imaging files will be publicly viewable on nova.astrometry.net.",
                "Standard Filter": "To use EXOTIC standard filters, type only the filter name.",
                "Custom Filter": "To use a custom filter, enter in the FWHM in optional_info.",
                "Target Star RA": "Must be in HH:MM:SS sexagesimal format.",
                "Target Star DEC": "Must be in +/-DD:MM:SS sexagesimal format with correct sign at the beginning (+ or -).",
                "Demosaic Format": "Optional control for handling Bayer pattern color images - to use, provide Bayer color patttern of your camera (RGGB, BGGR, GRBG, GBRG) - null (no color processing) is default",
                "Demosaic Output": "Select how to process color data (gray for grayscale, red or green or blue for single color channel, blueblock for grayscale without blue, [ R, G, B ] for custom weights for mixing colors.  green is default",
                "Formatting of null": "Due to the file being a .json, null is case sensitive and must be spelled as shown.",
                "Decimal Format": "Leading zero must be included when appropriate (Ex: 0.32, .32 or 00.32 causes errors.)."
            }

            null = None

            if fitsortext.get() == 1:
                if obsinfo.get() == 'manual':
                    new_inits['user_info'] = {
                        "Directory with FITS files": FITS_dir.folder_path,
                        "Directory to Save Plots": save_dir.folder_path,
                        # "Directory of Flats": flats_dir.folder_path,
                        # "Directory of Darks": darks_dir.folder_path,
                        # "Directory of Biases": biases_dir.folder_path,

                        "AAVSO Observer Code (blank if none)": input_data['obscode'],
                        "Secondary Observer Codes (blank if none)": input_data['secondobscode'],

                        "Observation date": input_data['obsdate'],
                        "Obs. Latitude": input_data['lat'],
                        "Obs. Longitude": input_data['long'],
                        "Obs. Elevation (meters; Note: leave blank if unknown)": float(input_data.get('elevation')),
                        "Camera Type (CCD or DSLR)": input_data['cameratype'],
                        "Pixel Binning": input_data['pixbin'],
                        "Filter Name (aavso.org/filters)": input_data['obsfilter'],
                        "Observing Notes": input_data['obsnotes'],

                        "Plate Solution? (y/n)": input_data['platesolve'],
                        "Add Comparison Stars from AAVSO? (y/n)": input_data['aavso_comp'],

                        "Target Star X & Y Pixel": (input_data['targetpos']),
                        "Comparison Star(s) X & Y Pixel": (input_data['comppos']),
                        
                        "Demosaic Format": null, # TODO add GUI input for these
                        "Demosaic Output": null
                    }

                    if flats_dir.folder_path == 'null':
                        new_inits['user_info']["Directory of Flats"] = null
                    else:
                        new_inits['user_info']["Directory of Flats"] = flats_dir.folder_path

                    if darks_dir.folder_path == 'null':
                        new_inits['user_info']["Directory of Darks"] = null
                    else:
                        new_inits['user_info']["Directory of Darks"] = darks_dir.folder_path

                    if biases_dir.folder_path == 'null':
                        new_inits['user_info']["Directory of Biases"] = null
                    else:
                        new_inits['user_info']["Directory of Biases"] = biases_dir.folder_path

                elif obsinfo.get() == 'inits':
                    with open(input_data['inits_dir'], "r") as confirmed:
                        original_inits = json.load(confirmed)

                    new_inits['user_info'] = original_inits['user_info']

                new_inits['optional_info'] = {
                    "Filter Minimum Wavelength (nm)": input_data.get('filtermin', null),
                    "Filter Maximum Wavelength (nm)": input_data.get('filtermax', null)
                }

                if 'pixscale' not in input_data.keys():
                    new_inits['optional_info']['Pixel Scale (Ex: 5.21 arcsecs/pixel)'] = null
                elif input_data['pixscale'] == 'null':
                    new_inits['optional_info']['Pixel Scale (Ex: 5.21 arcsecs/pixel)'] = null
                else:
                    new_inits['optional_info']['Pixel Scale (Ex: 5.21 arcsecs/pixel)'] = input_data['pixscale']

            elif fitsortext.get() == 2:

                if obsinfo.get() == 'manual':
                    new_inits['user_info'] = {
                        "Directory to Save Plots": save_dir.folder_path,

                        "AAVSO Observer Code (blank if none)": input_data['obscode'],
                        "Secondary Observer Codes (blank if none)": input_data['secondobscode'],

                        "Observation date": input_data['obsdate'],
                        "Obs. Latitude": input_data['lat'],
                        "Obs. Longitude": input_data['long'],
                        "Obs. Elevation (meters)": float(input_data.get('elevation', 0)),
                        "Camera Type (CCD or DSLR)": input_data['cameratype'],
                        "Pixel Binning": input_data['pixbin'],
                        "Filter Name (aavso.org/filters)": input_data['obsfilter'],
                        "Observing Notes": input_data['obsnotes'],
                    }

                elif obsinfo.get() == 'inits':
                    with open(input_data['inits_dir'], "r") as confirmed:
                        original_inits = json.load(confirmed)

                    new_inits['user_info'] = original_inits['user_info']

                new_inits['optional_info'] = {
                    "Pre-reduced File:": prered_file.file_path,
                    "Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)": input_data['file_time'],
                    "Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)": input_data['file_units'],
                    "Comparison Star used in Photometry (blank if none)": input_data['phot_comp_star'],
                    "Exposure Time (s)": input_data['exp']
                }

            if planetparams.get() in ["manual", "nea"]:
                new_inits['planetary_parameters'] = {
                    "Target Star RA": input_data['ra'],
                    "Target Star Dec": input_data['dec'],
                    "Planet Name": input_data['pName'],
                    "Host Star Name": input_data['sName'],
                    "Orbital Period (days)": float(input_data['pPer']),
                    "Orbital Period Uncertainty": float(input_data['pPerUnc']),
                    "Published Mid-Transit Time (BJD-UTC)": float(input_data['midT']),
                    "Mid-Transit Time Uncertainty": float(input_data['midTUnc']),
                    "Ratio of Planet to Stellar Radius (Rp/Rs)": float(input_data['rprs']),
                    "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": float(input_data['rprsUnc']),
                    "Ratio of Distance to Stellar Radius (a/Rs)": float(input_data['aRs']),
                    "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": float(input_data['aRsUnc']),
                    "Orbital Inclination (deg)": float(input_data['inc']),
                    "Orbital Inclination (deg) Uncertainty": float(input_data['incUnc']),
                    "Orbital Eccentricity (0 if null)": float(input_data['ecc']),
                    "Argument of Periastron (deg)": float(input_data['omega']),
                    "Star Effective Temperature (K)": float(input_data['teff']),
                    "Star Effective Temperature (+) Uncertainty": float(input_data['teffUncPos']),
                    "Star Effective Temperature (-) Uncertainty": float(input_data['teffUncNeg']),
                    "Star Metallicity ([FE/H])": float(input_data['met']),
                    "Star Metallicity (+) Uncertainty": float(input_data['metUncPos']),
                    "Star Metallicity (-) Uncertainty": float(input_data['metUncNeg']),
                    "Star Surface Gravity (log(g))": float(input_data['logg']),
                    "Star Surface Gravity (+) Uncertainty": float(input_data['loggUncPos']),
                    "Star Surface Gravity (-) Uncertainty": float(input_data['loggUncNeg'])
                }

            elif planetparams.get() == "inits":
                with open(input_data['inits_dir'], "r") as confirmed:
                    original_inits = json.load(confirmed)

                new_inits['planetary_parameters'] = original_inits['planetary_parameters']

            now = datetime.now()
            dt_string = now.strftime("%d_%m_%Y__%H_%M_%S")
            fname = os.path.join(input_data['initssave_dir'], "inits_" + dt_string + ".json")
            with open(fname, "w") as initsf:
                json.dump(new_inits, initsf, indent=4)
            print(f"{fname} saved!")

            root=tk.Tk() 
            root.protocol("WM_DELETE_WINDOW", exit)
            root.title(f"EXOTIC v{__version__}")

            window_label = tk.Label(root,
                                    text="Your initialization file has been created!",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            window_label = tk.Label(root,
                                    text=f"{fname}",
                                    font="Helvetica 14 bold",
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=1, column=0, sticky=tk.N, pady=6)

            # Button for closing
            exit_button = tk.Button(root, text="Run EXOTIC", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=2, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        run_mthd = None

        if fitsortext.get() == 1:
            run_mthd = '--reduce'
        elif fitsortext.get() == 2:
            run_mthd = '--prereduced'

        #         If the user already has an inits file, then go for it
        try:
            if planetparams.get() == 'inits':
                try:
                    subprocess.run(['exotic', run_mthd, inits_dir.file_path, '-ov'], check=True)
                except:
                    try:
                        subprocess.run(['python3', 'exotic.py', run_mthd, inits_dir.file_path, '-ov'], check=True)
                    except:
                        subprocess.run(['python3', 'exotic/exotic.py', run_mthd, inits_dir.file_path, '-ov'],
                                       check=True)

            elif planetparams.get() == 'nea':
                try:
                    subprocess.run(['exotic', run_mthd, fname, '-ov'], check=True)
                except:
                    try:
                        subprocess.run(['python3', 'exotic.py', run_mthd, fname, '-ov'], check=True)
                    except:
                        subprocess.run(['python3', 'exotic/exotic.py', run_mthd, fname, '-ov'], check=True)

            else:
                try:
                    subprocess.run(['exotic', run_mthd, fname], check=True)
                except:
                    try:
                        subprocess.run(['python3', 'exotic.py', run_mthd, fname], check=True)
                    except:
                        subprocess.run(['python3', 'exotic/exotic.py', run_mthd, fname], check=True)
        except:
            print("\n\n################################################")
            print("ERROR: Please contact the Exoplanet Watch Team for help "
                  "on our Slack Workspace in the #data-reduction channel!")
            print("You can sign up for our free Slack Workspace here: "
                  "https://join.slack.com/t/uol-ets/shared_invite/zt-mvb4ljbo-LRBgpk3uMmUokbs4ge2JlA")
            print("################################################\n\n")
            pass


if __name__ == "__main__":
    main()
