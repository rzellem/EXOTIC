import tkinter as tk
# from exotic import NASAExoplanetArchive

try:  # simple version
    from .version import __version__
except ImportError:  # package import
    from version import __version__

import tkinter as tk

root = tk.Tk()
root.title("EXOTIC")

reduction_opt = tk.IntVar()

tk.Label(root,
        text="""How do you want to run EXOTIC?""",
        justify = tk.LEFT,
        padx = 20).pack()

tk.Radiobutton(root,
               text="Real Time Reduction\n(for analyzing your data while observing) ",
               padx = 20,
               variable=reduction_opt,
               value=1).pack(anchor=tk.CENTER)

tk.Radiobutton(root,
               text="Complete Reduction\n(for analyzing your data after an observing run).",
               padx = 20,
               variable=reduction_opt,
               value=2).pack(anchor=tk.CENTER)

# Button for closing
exit_button = tk.Button(root, text="Next", command=root.destroy)
exit_button.pack(pady=20)

root.mainloop()


if reduction_opt.get() == 2:
    root = tk.Tk()
    root.title("EXOTIC")

    fitsortext = tk.IntVar()

    tk.Label(root,
             text="""How do you want to run EXOTIC?""",
             justify=tk.LEFT,
             padx=20).pack()

    tk.Radiobutton(root,
                   text="Perform aperture photometry on image files",
                   padx=20,
                   variable=fitsortext,
                   value=1).pack(anchor=tk.W)

    tk.Radiobutton(root,
                   text="Start with pre-reduced data in a .txt format",
                   padx=20,
                   variable=fitsortext,
                   value=2).pack(anchor=tk.W)

    # Button for closing
    exit_button = tk.Button(root, text="Next", command=root.destroy)
    exit_button.pack(pady=20)

    root.mainloop()

    if fitsortext.get() == 1:
        root = tk.Tk()
        root.title("EXOTIC")

        initparams = tk.IntVar()

        tk.Label(root,
                 text="""How would you like to input your initial parameters?""",
                 justify=tk.LEFT,
                 padx=20).pack()

        tk.Radiobutton(root,
                       text="Manually",
                       padx=20,
                       variable=initparams,
                       value=1).pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="Automatically adopt all planetary parameters\nfrom the NASA Exoplanet Archive",
                       padx=20,
                       variable=initparams,
                       value=2).pack(anchor=tk.W)

        tk.Radiobutton(root,
                       text="From a pre-existing input file (e.g., inits.json)",
                       padx=20,
                       variable=initparams,
                       value=3).pack(anchor=tk.W)

        # Button for closing
        exit_button = tk.Button(root, text="Next", command=root.destroy)
        exit_button.pack(pady=20)

        root.mainloop()

        if initparams.get() == 1:
            root = tk.Tk()
            root.title("EXOTIC")

            initparams = tk.IntVar()

            window_label = tk.Label(root,
                     text="""Please enter the following information:""",
                     font = ("Helvetica 15 bold"),
                     justify=tk.LEFT,
                     padx=20)#.pack()
            window_label.grid(row=0, column=0, sticky=tk.N, pady=6)

            # Set up rows + columns
            i=1; j=0
            # "Directory with FITS files": "sample-data/HatP32Dec202017",
            inputdir_label = tk.Label(root, text="Directory with FITS files", font=("Helvetica 12 bold"), justify=tk.LEFT)
            inputdir_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            inputdir_label.grid(row=i, column= j, sticky=tk.W, pady=2)
            inputdir_entry.grid(row=i, column=j+1, sticky=tk.W, pady=2)
            i+=1

#             # "Directory to Save Plots": "sample-data/",
            savedir_label = tk.Label(root, text="Directory to Save EXOTIC Output", font=("Helvetica 12 bold"), justify=tk.LEFT)
            savedir_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            savedir_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            savedir_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i+=1

#             # "Directory of Flats": null,
            flatdir_label = tk.Label(root, text="Directory of Flat(s) (null if none)", font=("Helvetica 12 bold"), justify=tk.LEFT)
            flatdir_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            flatdir_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            flatdir_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i+=1
#
#             # "Directory of Darks": null,
            darkdir_label = tk.Label(root, text="Directory of Dark(s) (null if none)", font=("Helvetica 12 bold"), justify=tk.LEFT)
            darkdir_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            darkdir_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            darkdir_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
#             # "Directory of Biases": null,
            biasdir_label = tk.Label(root, text="Directory of Bias(es) (null if none)", font=("Helvetica 12 bold"), justify=tk.LEFT)
            biasdir_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            biasdir_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            biasdir_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
#             # "AAVSO Observer Code (N/A if none)": "RTZ",
            obscode_label = tk.Label(root, text="AAVSO Observer Code (N/A if none)", font=("Helvetica 12 bold"), justify=tk.LEFT)
            obscode_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            obscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
#
#             # "Secondary Observer Codes (N/A if none)": "N/A",
            secondobscode_label = tk.Label(root, text="Secondary Observer Codes (N/A if none)", font=("Helvetica 12 bold"), justify=tk.LEFT)
            secondobscode_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            secondobscode_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            secondobscode_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
#             # "Observation date": "17-December-2017",
            obsdate_label = tk.Label(root, text="Observation date", font=("Helvetica 12 bold"), justify=tk.LEFT)
            obsdate_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            obsdate_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            obsdate_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
# #             "Obs. Latitude": "+32.41638889",
            lat_label = tk.Label(root, text="Obs. Latitude", font=("Helvetica 12 bold"), justify=tk.LEFT)
            lat_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            lat_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            lat_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
# #             "Obs. Longitude": "-110.73444444",
            long_label = tk.Label(root, text="Obs. Longitude", font=("Helvetica 12 bold"), justify=tk.LEFT)
            long_entry = tk.Entry(root, font=("Helvetica 12"), justify=tk.LEFT)
            long_label.grid(row=i, column=j, sticky=tk.W, pady=2)
            long_entry.grid(row=i, column=j + 1, sticky=tk.W, pady=2)
            i += 1
#
# #             "Obs. Elevation (meters)": 2616,
#             elevation_frame = tk.Frame(root)
#             tk.Label(elevation_frame, text="Obs. Elevation (meters)", font=("Helvetica 15 bold")).pack(side=tk.LEFT)
#             word = tk.Entry(elevation_frame, font=("Helvetica 15 bold"))
#             word.pack()
#             elevation_frame.pack(pady=10)
#
# #             "Camera Type (CCD or DSLR)": "CCD",
#             cameratype_frame = tk.Frame(root)
#             tk.Label(cameratype_frame, text="Camera Type (CCD or DSLR)", font=("Helvetica 15 bold")).pack(side=tk.LEFT)
#             word = tk.Entry(cameratype_frame, font=("Helvetica 15 bold"))
#             word.pack()
#             cameratype_frame.pack(pady=10)
#
# #             "Pixel Binning": "1x1",
#             pixbinning_frame = tk.Frame(root)
#             tk.Label(pixbinning_frame, text="Pixel Binning", font=("Helvetica 15 bold")).pack(side=tk.LEFT)
#             word = tk.Entry(pixbinning_frame, font=("Helvetica 15 bold"))
#             word.pack()
#             pixbinning_frame.pack(pady=10)
#
# #             "Filter Name (aavso.org/filters)": "V",
#             filter_frame = tk.Frame(root)
#             tk.Label(filter_frame, text="Filter Name (aavso.org/filters)", font=("Helvetica 15 bold")).pack(side=tk.LEFT)
#             word = tk.Entry(filter_frame, font=("Helvetica 15 bold"))
#             word.pack()
#             filter_frame.pack(pady=10)
#
# #             "Observing Notes": "Weather, seeing was nice.",
#             obsnotes_frame = tk.Frame(root)
#             tk.Label(obsnotes_frame, text="Observing Notes", font=("Helvetica 15 bold")).pack(side=tk.LEFT)
#             word = tk.Entry(obsnotes_frame, font=("Helvetica 15 bold"))
#             word.pack()
#             obsnotes_frame.pack(pady=10)
#
#             "Plate Solution? (y/n)": "y",
#             "Align Images? (y/n)": "y",
#
#             "Target Star X & Y Pixel": [424, 286],
#             "Comparison Star(s) X & Y Pixel": [[465, 183], [512, 263]]
#         },
#         "planetary_parameters": {
#         "Target Star RA": "02:04:10",
#         "Target Star Dec": "+46:41:23",
#         "Planet Name": "HAT-P-32 b",
#         "Host Star Name": "HAT-P-32",
#         "Orbital Period (days)": 2.1500082,
#         "Orbital Period Uncertainty": 1.3e-07,
#         "Published Mid-Transit Time (BJD-UTC)": 2455867.402743,
#         "Mid-Transit Time Uncertainty": 4.9e-05,
#         "Ratio of Planet to Stellar Radius (Rp/Rs)": 0.14886235252742716,
#         "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": 0.0005539487393037134,
#         "Ratio of Distance to Stellar Radius (a/Rs)": 5.344,
#         "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": 0.039496835316262996,
#         "Orbital Inclination (deg)": 88.98,
#         "Orbital Inclination (deg) Uncertainty": 0.7602631123499285,
#         "Orbital Eccentricity (0 if null)": 0.159,
#         "Star Effective Temperature (K)": 6001.0,
#         "Star Effective Temperature (+) Uncertainty": 88.0,
#         "Star Effective Temperature (-) Uncertainty": -88.0,
#         "Star Metallicity ([FE/H])": -0.16,
#         "Star Metallicity (+) Uncertainty": 0.08,
#         "Star Metallicity (-) Uncertainty": -0.08,
#         "Star Surface Gravity (log(g))": 4.22,
#         "Star Surface Gravity (+) Uncertainty": 0.04,
#         "Star Surface Gravity (-) Uncertainty": -0.04
#     },
#     "optional_info": {
#     "Pixel Scale (Ex: 5.21 arcsecs/pixel)": null,
#     "Filter Minimum Wavelength (nm)": null,
#     "Filter Maximum Wavelength (nm)": null
# }

            tk.mainloop()

        elif initparams.get() == 2:
            root = tk.Tk()
            root.title("EXOTIC")

            initparams = tk.IntVar()

            tk.Label(root,
                     text="""How would you like to input your parameters?""",
                     justify=tk.LEFT,
                     padx=20).pack()

            tk.Radiobutton(root,
                           text="Manually",
                           padx=20,
                           variable=initparams,
                           value=1).pack(anchor=tk.W)


            tk.Radiobutton(root,
                           text="From a pre-existing input file (e.g., inits.json)",
                           padx=20,
                           variable=initparams,
                           value=2).pack(anchor=tk.W)

            # Button for closing
            exit_button = tk.Button(root, text="Next", command=root.destroy)
            exit_button.pack(pady=20)

            root.mainloop()

        else:
            print('Ok')


import pdb; pdb.set_trace()