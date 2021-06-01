# from exotic import NASAExoplanetArchive
from tkinter import filedialog
import tkinter as tk
import exotic
import os
import json

try:  # filters
    from .api.filters import fwhm as photometric_filters
except ImportError:  # package import
    from api.filters import fwhm as photometric_filters

try:  # simple version
    from .version import __version__
except ImportError:  # package import
    from version import __version__

class FolderSelect(tk.Frame):
    def __init__(self,parent=None,folderDescription="",default_text="",**kw):
        tk.Frame.__init__(self,master=parent,**kw)
        self.folderPath = tk.StringVar()
        self.lblName = tk.Label(self, text=folderDescription, anchor=tk.W)
        self.lblName.grid(row=0,column=0)
        self.entPath = tk.Entry(self, textvariable=self.folderPath)
        self.entPath.insert(tk.END, default_text)
        self.entPath.grid(row=0,column=1)
        self.btnFind = tk.Button(self, text="Browse Folder",command=self.setFolderPath, anchor=tk.W)
        self.btnFind.grid(row=0,column=2)
    def setFolderPath(self):
        folder_selected = filedialog.askdirectory()
        self.folderPath.set(folder_selected)
    @property
    def folder_path(self):
        return self.folderPath.get()

class FileSelect(tk.Frame):
    def __init__(self,parent=None,folderDescription="",default_text="",**kw):
        tk.Frame.__init__(self,master=parent,**kw)
        self.filePath = tk.StringVar()
        self.lblName = tk.Label(self, text=folderDescription, anchor=tk.W)
        self.lblName.grid(row=0,column=0)
        self.entPath = tk.Entry(self, textvariable=self.filePath)
        self.entPath.insert(tk.END, default_text)
        self.entPath.grid(row=0,column=1)
        self.btnFind = tk.Button(self, text="Browse",command=self.setFilePath, anchor=tk.W)
        self.btnFind.grid(row=0,column=2)
    def setFilePath(self):
        file_selected = filedialog.askopenfilename()
        self.filePath.set(file_selected)
    @property
    def file_path(self):
        return self.filePath.get()

# def doStuff():
#     folder1 = directory1Select.folder_path
#     folder2 = directory2Select.folder_path
#     folder3 = directory3Select.folder_path
#     print("Doing stuff with folder", folder1, folder2, folder3)

root = tk.Tk()
root.title(f"EXOTIC v{__version__}")

tk.Label(root,
        text="Welcome to EXOTIC!",
        justify = tk.LEFT,
        font=("Arial", 36),
        padx = 20).pack(anchor=tk.CENTER)
tk.Label(root,
        text="the EXOplanet Transit Interpretation Code",
        justify = tk.LEFT,
        font=("Arial", 22),
        padx = 20).pack(anchor=tk.CENTER)
tk.Label(root,
        text=f"Version {__version__}",
        justify = tk.LEFT,
        padx = 20).pack(anchor=tk.CENTER)

tk.Label(root,
        text="Copyright (c) 2021, California Institute of Technology. All rights reserved.\nBased on Government Sponsored Research under contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.",
        font=("Arial", 12),
        justify = tk.CENTER,
        padx = 20).pack(anchor=tk.CENTER)

# Button for closing
exit_button = tk.Button(root, text="Next", command=root.destroy)
exit_button.pack(anchor=tk.CENTER)

root.mainloop()

root = tk.Tk()
root.title(f"EXOTIC v{__version__}")

reduction_opt = tk.IntVar()
tk.Label(root,
        text="\nHow would you like to run EXOTIC?",
        justify = tk.LEFT,
        font=("Helvetica 14 bold"),
        padx = 20).pack(anchor=tk.W)

tk.Radiobutton(root,
               text="Real Time Reduction\n(for analyzing your data while simultaneously observing)",
               justify = tk.LEFT,
               padx = 20,
               variable=reduction_opt,
               state=tk.DISABLED,
               value=1).pack(anchor=tk.W)

tk.Radiobutton(root,
               text="Complete Reduction\n(for analyzing your data after an observing run)",
               justify = tk.LEFT,
               padx = 20,
               variable=reduction_opt,
               value=2).pack(anchor=tk.W)

# Button for closing
exit_button = tk.Button(root, text="Next", command=root.destroy)
exit_button.pack(anchor=tk.E)

root.mainloop()


if reduction_opt.get() == 2:
    root = tk.Tk()
    root.title(f"EXOTIC v{__version__}")

    fitsortext = tk.IntVar()

    tk.Label(root,
             text="""How do you want to run EXOTIC?""",
             font=("Helvetica 14 bold"),
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
                   state=tk.DISABLED,
                   value=2).pack(anchor=tk.W)

    # Button for closing
    exit_button = tk.Button(root, text="Next", command=root.destroy)
    exit_button.pack(pady=20,anchor=tk.E)

    root.mainloop()

    if fitsortext.get() == 1:
        # First ask user how they want to enter the observing information
        root = tk.Tk()
        root.title(f"EXOTIC v{__version__}")

        obsinfo = tk.StringVar()

        tk.Label(root,
                 text="""How would you like to input your observing information?""",
                 justify=tk.LEFT,
                 font=("Helvetica 14 bold"),
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually",
                       padx=20,
                       variable=obsinfo,
                       value="manual").pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g., inits.json)",
                       padx=20,
                       variable=obsinfo,
                       value='inits').pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20,anchor=tk.E)

        root.mainloop()



        if obsinfo.get() == "manual":
            root = tk.Tk()
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="""Please enter the following information about your observation:""",
                                    font=("Helvetica 14 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1;
            j = 0

            folderPath = tk.StringVar()

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            FITS_dir = FolderSelect(root, "Folder with your FITS files")
            FITS_dir.grid(row=i)
            i += 1

            # "Directory to Save Plots": "sample-data/",
            save_dir = FolderSelect(root, "Folder to save EXOTIC output")
            save_dir.grid(row=i)
            i += 1

            # "Directory of Flats": null,
            flats_dir = FolderSelect(root, "Folder with your flats (null if none)", "null")
            flats_dir.grid(row=i)
            i += 1

            # "Directory of Darks": null,
            darks_dir = FolderSelect(root, "Folder with your darks (null if none)", "null")
            darks_dir.grid(row=i)
            i += 1

            # "Directory of Biases": null,
            biases_dir = FolderSelect(root, "Folder with your biases (null if none)", "null")
            biases_dir.grid(row=i)
            i += 1

            # "Directory of Biases": null,
            biases_dir = FolderSelect(root, "Folder with your biases (null if none)", "null")
            biases_dir.grid(row=i)
            i += 1

            #
            # obscode2_label = TextInput(root, "AAVSO Obs Code", "RTZ")
            # obscode2_label.grid(row=1)
            # i += 1


            #             # "AAVSO Observer Code (N/A if none)": "RTZ",
            obscode_label = tk.Label(root, text="AAVSO Observer Code (N/A if none)", justify=tk.LEFT)
            obscode_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            obscode_entry.insert(tk.END, "N/A")
            obscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            obscode = obscode_entry.get()
            i += 1
            #
            #
            #             # "Secondary Observer Codes (N/A if none)": "N/A",
            secondobscode_label = tk.Label(root, text="Secondary Observer Codes (N/A if none)", justify=tk.LEFT)
            secondobscode_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            secondobscode_entry.insert(tk.END, "N/A")
            secondobscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            secondobscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            secondobscode = secondobscode_entry.get()
            i += 1
            #
            #             # "Observation date": "17-December-2017",
            obsdate_label = tk.Label(root, text="Observation date", justify=tk.LEFT)
            obsdate_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            obsdate_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obsdate_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            obsdate = obsdate_entry.get()
            i += 1
            #
            # #             "Obs. Latitude": "+32.41638889",
            lat_label = tk.Label(root, text="Obs. Latitude", justify=tk.LEFT)
            lat_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            lat_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            lat_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            lat = lat_entry.get()
            i += 1
            #
            # #             "Obs. Longitude": "-110.73444444",
            long_label = tk.Label(root, text="Obs. Longitude", justify=tk.LEFT)
            long_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            long_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            long_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            long = long_entry.get()
            i += 1
            #
            # #             "Obs. Elevation (meters)": 2616,
            elevation_label = tk.Label(root, text="Obs. Elevation", justify=tk.LEFT)
            elevation_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            elevation_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            elevation_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            elevation = elevation_entry.get()
            i += 1
            #
            # #             "Camera Type (CCD or DSLR)": "CCD",
            cameratype_label = tk.Label(root, text="Camera Type (CCD or DSLR)", justify=tk.LEFT)
            cameratype_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            cameratype_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            cameratype_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            cameratype = cameratype_entry.get()
            i += 1
            #
            # #             "Pixel Binning": "1x1",
            pixbin_label = tk.Label(root, text="Pixel Binning", justify=tk.LEFT)
            pixbin_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            pixbin_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            pixbin_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            pixbin = pixbin_entry.get()
            i += 1

            #     "Pixel Scale (Ex: 5.21 arcsecs/pixel)": null,
            pixscale_label = tk.Label(root, text="Pixel Scale (null if unknown)", justify=tk.LEFT)
            pixscale_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            pixscale_entry.insert(tk.END, "5.21 arcsecs/pixel")
            pixscale_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            pixscale_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            pixscale = pixscale_entry.get()
            i += 1

            # # #             "Filter Name (aavso.org/filters)": "V",
            # filtername_label = tk.Label(root, text="Filter Name (please see aavso.org/filters)", justify=tk.LEFT)
            # filtername_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            # filtername_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            # filtername_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            # i += 1

            # choices = ['one', 'two', 'three']
            choices = [item for sublist in photometric_filters for item in sublist]

            filteroptions = tk.StringVar(root)
            filteroptions.set(choices[0])  # default value

            l3 = tk.Label(root, text='Filter', width=15)
            l3.grid(row=i, column=0)

            om1 = tk.OptionMenu(root, filteroptions, *choices)
            om1.grid(row=i, column=1)
            # my_w.mainloop()  # Keep the window open
            obsfilter = filteroptions.get()
            i += 1

            # #             "Observing Notes": "Weather, seeing was nice.",
            obsnotes_label = tk.Label(root, text="Observing Notes", justify=tk.LEFT)
            obsnotes_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            obsnotes_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obsnotes_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            obsnotes = obsnotes_entry.get()
            i += 1


            # root = tk.Tk()
            # root.title(f"EXOTIC v{__version__}")
            #
            # window_label = tk.Label(root,
            #                         text="""Please enter the following information:""",
            #                         font=("Helvetica 15 bold"),
            #                         justify=tk.LEFT,
            #                         padx=20)  # .pack()
            # window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # i, j = 1, 0
            platesolve = tk.BooleanVar()
            platesolve_check = tk.Checkbutton(root, text='Plate solve my images', variable=platesolve, onvalue=True,
                                              offvalue=False)
            platesolve_check.grid(row=i, column=j, sticky=tk.W, pady=2)
            i += 1

            alignment = tk.BooleanVar()
            alignment_check = tk.Checkbutton(root, text='Align my images', variable=alignment, onvalue=1, offvalue=0)
            alignment_check.grid(row=i, column=j, sticky=tk.W, pady=2)
            i += 1

            targetpos_label = tk.Label(root, text="Target Star X & Y Pixel Position", justify=tk.LEFT)
            targetpos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            targetpos_entry.insert(tk.END, "[x, y]")
            targetpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            targetpos = targetpos_entry.get()
            i += 1

            comppos_label = tk.Label(root, text="Comparison Star(s) X & Y Pixel Position", justify=tk.LEFT)
            comppos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            comppos_entry.insert(tk.END, "[x1, y1], [x2, y2]")
            comppos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            comppos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            comppos = comppos_entry.get()
            i += 1

            # # Button for closing
            # exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)
            # # exit_button.pack(pady=20)
            # root.update()
            # root.mainloop()


            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        # else:
        #     root = tk.Tk()
        #     root.title(f"EXOTIC v{__version__}")
        #
        #     filePath = tk.StringVar()
        #
        #     # # "Directory with FITS files": "sample-data/HatP32Dec202017",
        #     inits_dir = FileSelect(root, "Please select your initialization file")
        #     inits_dir.grid(row=0)
        #
        #     # Button for closing
        #     exit_button = tk.Button(root, text="Next", command=root.destroy)
        #     # exit_button.pack(pady=20)
        #     exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)
        #
        #     tk.mainloop()

        try:
            if filteroptions.get() == "N/A":
                root = tk.Tk()
                root.title(f"EXOTIC v{__version__}")

                i = 0; j = 0

                tk.Label(root,
                         text="Please enter the following information for your custom filter:",
                         font=("Helvetica 15 bold"),
                         justify=tk.LEFT,
                         padx=20).grid(row=i, column=j, sticky=tk.W, pady=2)
                i += 1

                # #
                # "Filter Minimum Wavelength (nm)": null,
                filtermin_label = tk.Label(root, text="Filter Minimum Wavelength (nm)", justify=tk.LEFT)
                filtermin_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
                filtermin_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                filtermin_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                filtermin = filtermin_entry.get()
                i += 1

                # "Filter Maximum Wavelength (nm)": null
                filtermax_label = tk.Label(root, text="Filter Maximum Wavelength (nm)", justify=tk.LEFT)
                filtermax_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
                filtermax_label.grid(row=i, column=j, sticky=tk.W, pady=2)
                filtermax_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
                filtermax = filtermax_entry.get()
                i += 1

                # Button for closing
                exit_button = tk.Button(root, text="Next", command=root.destroy)
                # exit_button.pack(pady=20)
                exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

                root.mainloop()
        except:
            pass

        # then ask user how they want to enter the planetary information
        root = tk.Tk()
        root.title(f"EXOTIC v{__version__}")

        planetparams = tk.StringVar()

        tk.Label(root,
                 text="""How would you like to input your planetary system parameters?""",
                 font=("Helvetica 15 bold"),
                 justify=tk.LEFT,
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually",
                       padx=20,
                       variable=planetparams,
                       value="manual").pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="Automatically adopt all planetary parameters\nfrom the NASA Exoplanet Archive",
                       padx=20,
                       variable=planetparams,
                       value="nea").pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g., inits.json)",
                       padx=20,
                       variable=planetparams,
                       value="inits").pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20, anchor=tk.E)

        root.mainloop()

        if planetparams.get() == "manual":
            root = tk.Tk()
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="""Please enter the following information:""",
                                    font=("Helvetica 15 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1;
            j = 0

            #         "Planet Name": "HAT-P-32 b",
            planet_label = tk.Label(root, text="Planet Name", justify=tk.LEFT)
            planet_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            planet_entry.insert(tk.END, "HAT-P-32 b")
            planet_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            planet_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            planet = planet_entry.get()
            i += 1

            #         "Host Star Name": "HAT-P-32",
            star_label = tk.Label(root, text="Host Star Name", justify=tk.LEFT)
            star_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            star_entry.insert(tk.END, "HAT-P-32")
            star_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            star_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            star = star_entry.get()
            i += 1

            #         "Target Star RA": "02:04:10",
            targetRA_label = tk.Label(root, text="Host Star Right Ascension", justify=tk.LEFT)
            targetRA_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            targetRA_entry.insert(tk.END, "02:04:10")
            targetRA_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetRA_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            targetRA = targetRA_entry.get()
            i += 1

            #         "Target Star Dec": "+46:41:23",
            targetDEC_label = tk.Label(root, text="Host Star Declination", justify=tk.LEFT)
            targetDEC_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            targetDEC_entry.insert(tk.END, "+46:41:23")
            targetDEC_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            targetDEC_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            targetDEC = targetDEC_entry.get()
            i += 1

            #         "Orbital Period (days)": 2.1500082,
            period_label = tk.Label(root, text="Planet's orbital period (days)", justify=tk.LEFT)
            period_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            period_entry.insert(tk.END, "2.1500082")
            period_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            period_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            period = period_entry.get()
            i += 1

            #         "Orbital Period Uncertainty": 1.3e-07,
            perioderr_label = tk.Label(root, text="Planet's orbital period uncertainty (days)", justify=tk.LEFT)
            perioderr_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            perioderr_entry.insert(tk.END, "1.3e-07")
            perioderr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            perioderr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            perioderr = perioderr_entry.get()
            i += 1

            #         "Published Mid-Transit Time (BJD-UTC)": 2455867.402743,
            Tmid_label = tk.Label(root, text="Published Mid-Transit Time (BJD-UTC)", justify=tk.LEFT)
            Tmid_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            Tmid_entry.insert(tk.END, "2455867.402743")
            Tmid_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tmid_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            Tmid = Tmid_entry.get()
            i += 1

            #         "Mid-Transit Time Uncertainty": 4.9e-05,
            Tmiderr_label = tk.Label(root, text="Mid-Transit Time Uncertainty", justify=tk.LEFT)
            Tmiderr_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            Tmiderr_entry.insert(tk.END, "4.9e-05")
            Tmiderr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tmiderr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            Tmiderr = Tmiderr_entry.get()
            i += 1

            #         "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.14886235252742716,
            rprs_label = tk.Label(root, text="Ratio of Planet to Stellar Radius (Rp/Rs)", justify=tk.LEFT)
            rprs_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            rprs_entry.insert(tk.END, "0.14886235252742716")
            rprs_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            rprs_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            rprs = rprs_entry.get()
            i += 1

            #         "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.0005539487393037134,
            rprserr_label = tk.Label(root, text="Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty", justify=tk.LEFT)
            rprserr_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            rprserr_entry.insert(tk.END, "0.0005539487393037134")
            rprserr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            rprserr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            rprserr = rprserr_entry.get()
            i += 1

            #         "Ratio of Distance to Stellar Radius (a/Rs)": 5.344,
            aRs_label = tk.Label(root, text="Ratio of Distance to Stellar Radius (a/Rs)", justify=tk.LEFT)
            aRs_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            aRs_entry.insert(tk.END, "5.344")
            aRs_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            aRs_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            aRs = aRs_entry.get()
            i += 1

            #         "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.039496835316262996,
            aRserr_label = tk.Label(root, text="Ratio of Distance to Stellar Radius (a/Rs) Uncertainty", justify=tk.LEFT)
            aRserr_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            aRserr_entry.insert(tk.END, "0.039496835316262996")
            aRserr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            aRserr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            aRserr = aRserr_entry.get()
            i += 1

            #         "Orbital Inclination (deg)": 88.98,
            inc_label = tk.Label(root, text="Orbital Inclination (degrees)", justify=tk.LEFT)
            inc_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            inc_entry.insert(tk.END, "88.98")
            inc_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            inc_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            inc = inc_entry.get()
            i += 1

            #         "Orbital Inclination (deg) Uncertainty": 0.7602631123499285,
            incerr_label = tk.Label(root, text="Orbital Inclination (degrees) Uncertainty", justify=tk.LEFT)
            incerr_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            incerr_entry.insert(tk.END, "0.7602631123499285")
            incerr_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            incerr_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            incerr = incerr_entry.get()
            i += 1

            #         "Orbital Eccentricity (0 if null)": 0.159,
            ecc_label = tk.Label(root, text="Orbital Eccentricity (0 if null)", justify=tk.LEFT)
            ecc_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            ecc_entry.insert(tk.END, "0.159")
            ecc_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            ecc_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            ecc = ecc_entry.get()
            i += 1

            #         "Star Effective Temperature (K)": 6001.0,
            Teff_label = tk.Label(root, text="Star Effective Temperature (K)", justify=tk.LEFT)
            Teff_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            Teff_entry.insert(tk.END, "6001.0")
            Teff_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Teff_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            Teff = Teff_entry.get()
            i += 1

            #         "Star Effective Temperature (+) Uncertainty": 88.0,
            Tefferrpos_label = tk.Label(root, text="Star Effective Temperature Positive (+) Uncertainty", justify=tk.LEFT)
            Tefferrpos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            Tefferrpos_entry.insert(tk.END, "+88.0")
            Tefferrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tefferrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            Tefferrpos = Tefferrpos_entry.get()
            i += 1

            #         "Star Effective Temperature (-) Uncertainty": -88.0,
            Tefferrneg_label = tk.Label(root, text="Star Effective Temperature Negative (-) Uncertainty", justify=tk.LEFT)
            Tefferrneg_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            Tefferrneg_entry.insert(tk.END, "-88.0")
            Tefferrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            Tefferrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            Tefferrneg = Tefferrneg_entry.get()
            i += 1

            #         "Star Metallicity ([FE/H])": -0.16,
            FeH_label = tk.Label(root, text="Star Metallicity ([Fe/H])", justify=tk.LEFT)
            FeH_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            FeH_entry.insert(tk.END, "-0.1")
            FeH_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeH_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            FeH = FeH_entry.get()
            i += 1

            #         "Star Metallicity (+) Uncertainty": 0.08,
            FeHerrpos_label = tk.Label(root, text="Star Metallicity Positive (+) Uncertainty", justify=tk.LEFT)
            FeHerrpos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            FeHerrpos_entry.insert(tk.END, "0.08")
            FeHerrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeHerrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            FeHerrpos = FeHerrpos_entry.get()
            i += 1

            #         "Star Metallicity (-) Uncertainty": -0.08,
            FeHerrneg_label = tk.Label(root, text="Star Metallicity Negative (-) Uncertainty", justify=tk.LEFT)
            FeHerrneg_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            FeHerrneg_entry.insert(tk.END, "-0.08")
            FeHerrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            FeHerrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            FeHerrneg = FeHerrneg_entry.get()
            i += 1

            #         "Star Surface Gravity (log(g))": 4.22,
            logg_label = tk.Label(root, text="Star Surface Gravity (log(g))", justify=tk.LEFT)
            logg_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            logg_entry.insert(tk.END, "4.22")
            logg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            logg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            logg = logg_entry.get()
            i += 1

            #         "Star Surface Gravity (+) Uncertainty": 0.04,
            loggerrpos_label = tk.Label(root, text="Star Surface Gravity Positive (+) Uncertainty", justify=tk.LEFT)
            loggerrpos_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            loggerrpos_entry.insert(tk.END, "4.22")
            loggerrpos_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            loggerrpos_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            loggerrpos = loggerrpos_entry.get()
            i += 1

            #         "Star Surface Gravity (-) Uncertainty": -0.04
            loggerrneg_label = tk.Label(root, text="Star Surface Gravity Negative (-) Uncertainty", justify=tk.LEFT)
            loggerrneg_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            loggerrneg_entry.insert(tk.END, "4.22")
            loggerrneg_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            loggerrneg_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            loggerrneg = loggerrneg_entry.get()
            i += 1

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        if (planetparams.get() == "inits") or (obsinfo.get() == "inits"):
            root = tk.Tk()
            root.title(f"EXOTIC v{__version__}")

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            inits_dir = FileSelect(root, "Please select your initialization file")
            inits_dir.grid(row=0)

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=1, column=3, sticky=tk.W, pady=10)

            tk.mainloop()

        # if (obsinfo.get() == "manual"):
        #     root = tk.Tk()
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
        #     alignment = tk.BooleanVar()
        #     alignment_check = tk.Checkbutton(root, text='Align my images', variable=alignment, onvalue=1, offvalue=0)
        #     alignment_check.grid(row=i, column=j, sticky=tk.W, pady=2)
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
            root = tk.Tk()
            root.title(f"EXOTIC v{__version__}")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                                    text="To make analyzing these data in the future more easy, EXOTIC will create an initialization file for you.",
                                    font=("Helvetica 14 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i = 1;
            j = 0

            folderPath = tk.StringVar()

            # # "Directory with FITS files": "sample-data/HatP32Dec202017",
            initssave_dir = FolderSelect(root, "Folder to save your initialization file:")
            initssave_dir.grid(row=i)
            i += 1

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=i, column=3, sticky=tk.W, pady=10)

            tk.mainloop()


            #TODO: create the inits file here
            new_inits = {'inits_guide':{}, 'user_info':{}, 'planetary_parameters':{}, 'optional_info':{}}
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
                "Formatting of null": "Due to the file being a .json, null is case sensitive and must be spelled as shown.",
                "Decimal Format": "Leading zero must be included when appropriate (Ex: 0.32, .32 or 00.32 causes errors.)."
            }

            if obsinfo.get() == 'manual':
                new_inits['user_info'] = {
                    "Directory with FITS files": FITS_dir.folder_path,
                    "Directory to Save Plots": save_dir.folder_path,
                    "Directory of Flats": flats_dir.folder_path,
                    "Directory of Darks": darks_dir.folder_path,
                    "Directory of Biases": biases_dir.folder_path,

                    "AAVSO Observer Code (N/A if none)": obscode,
                    "Secondary Observer Codes (N/A if none)": secondobscode,

                    "Observation date": obsdate,
                    "Obs. Latitude": lat,
                    "Obs. Longitude": long,
                    "Obs. Elevation (meters)": float(elevation),
                    "Camera Type (CCD or DSLR)": cameratype,
                    "Pixel Binning": pixbin,
                    "Filter Name (aavso.org/filters)": obsfilter,
                    "Observing Notes": obsnotes,

                    "Plate Solution? (y/n)": platesolve,
                    "Align Images? (y/n)": alignment,

                    "Target Star X & Y Pixel": list(targetpos),
                    "Comparison Star(s) X & Y Pixel": list(comppos)
                }
            elif obsinfo.get() == 'inits':
                with open(inits_dir.file_path, "r") as confirmed:
                    original_inits = json.load(confirmed)

            if (planetparams.get() == "manual"):
                new_inits['planetary_parameters'] = {
                    "Target Star RA": targetRA,
                    "Target Star Dec": targetDEC,
                    "Planet Name": planet,
                    "Host Star Name": star,
                    "Orbital Period (days)": float(period),
                    "Orbital Period Uncertainty": float(perioderr),
                    "Published Mid-Transit Time (BJD-UTC)": float(Tmid),
                    "Mid-Transit Time Uncertainty": float(Tmiderr),
                    "Ratio of Planet to Stellar Radius (Rp/Rs)": float(rprs),
                    "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": float(rprserr),
                    "Ratio of Distance to Stellar Radius (a/Rs)": float(aRs),
                    "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": float(aRserr),
                    "Orbital Inclination (deg)": float(inc),
                    "Orbital Inclination (deg) Uncertainty": float(incerr),
                    "Orbital Eccentricity (0 if null)": float(ecc),
                    "Star Effective Temperature (K)": float(Teff),
                    "Star Effective Temperature (+) Uncertainty": float(Tefferrpos),
                    "Star Effective Temperature (-) Uncertainty": float(Tefferrneg),
                    "Star Metallicity ([FE/H])": float(FeH),
                    "Star Metallicity (+) Uncertainty": float(FeHerrpos),
                    "Star Metallicity (-) Uncertainty": float(FeHerrneg),
                    "Star Surface Gravity (log(g))": float(logg),
                    "Star Surface Gravity (+) Uncertainty": float(loggerrpos),
                    "Star Surface Gravity (-) Uncertainty": float(loggerrneg)
                }
            elif (planetparams.get() == "nea"):
                # Just put in dummy values as they will be overwritten by the NEA later
                new_inits['planetary_parameters'] = {
                    "Target Star RA": "00:00:00",
                    "Target Star Dec": "+00:00:00",
                    "Planet Name": planet,
                    "Host Star Name": star,
                    "Orbital Period (days)": 0.,
                    "Orbital Period Uncertainty": 0.,
                    "Published Mid-Transit Time (BJD-UTC)": 0.,
                    "Mid-Transit Time Uncertainty": 0.,
                    "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.,
                    "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.,
                    "Ratio of Distance to Stellar Radius (a/Rs)": 0.,
                    "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.,
                    "Orbital Inclination (deg)": 0.,
                    "Orbital Inclination (deg) Uncertainty": 0.,
                    "Orbital Eccentricity (0 if null)": 0.,
                    "Star Effective Temperature (K)": 0.,
                    "Star Effective Temperature (+) Uncertainty": 0.,
                    "Star Effective Temperature (-) Uncertainty": 0.,
                    "Star Metallicity ([FE/H])": 0.,
                    "Star Metallicity (+) Uncertainty": 0.,
                    "Star Metallicity (-) Uncertainty": 0.,
                    "Star Surface Gravity (log(g))": 0.,
                    "Star Surface Gravity (+) Uncertainty": 0.,
                    "Star Surface Gravity (-) Uncertainty": 0.
                }

            new_inits['optional_info'] = {
                "Pixel Scale (Ex: 5.21 arcsecs/pixel)": pixscale,
                "Filter Minimum Wavelength (nm)": filtermin,
                "Filter Maximum Wavelength (nm)": filtermax
            }


            with open(inits_dir.file_path, "w") as confirmed:
                data = json.dump(new_inits)

            # with open(inits_dir.file_path, "r") as confirmed:
            #     data = json.load(confirmed)

            root = tk.Tk()
            root.title(f"EXOTIC v{__version__}")

            window_label = tk.Label(root,
                                    text="Your initialization file has been created!",
                                    font=("Helvetica 14 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            window_label = tk.Label(root,
                                    text=f"{initsfile}",
                                    font=("Helvetica 14 bold"),
                                    justify=tk.LEFT,
                                    padx=20)  # .pack()
            window_label.grid(row=1, column=0, sticky=tk.N, pady=6)

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            # exit_button.pack(pady=20)
            exit_button.grid(row=2, column=3, sticky=tk.W, pady=10)

            tk.mainloop()


        #         If the user already has an inits file, then go for it
        try:
            inits_dir.file_path

            if planetparams.get() == 'inits':
                os.system(f'python3 exotic/exotic.py --reduce {inits_dir.file_path} -ov')
            elif planetparams.get() == 'nea':
                os.system(f'python3 exotic/exotic.py --reduce {inits_dir.file_path} -nea')
            else:
                os.system(f'python3 exotic/exotic.py --reduce {inits_dir.file_path}')


        except:
            pass

# import pdb; pdb.set_trace()