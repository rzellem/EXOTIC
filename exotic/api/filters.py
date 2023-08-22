# Sources for FWHM band wavelengths are referenced below, respectively 
# note that all below units default to nm
__fwhm = {
        # AAVSO, Source(s): AAVSO International Database; https://www.aavso.org/filters
        # Johnson
        "Johnson U": {"name": "U", "fwhm": ("333.8", "398.8")},
        "Johnson B": {"name": "B", "fwhm": ("391.6", "480.6")},
        "Johnson V": {"name": "V", "fwhm": ("502.8", "586.8")},
        "Johnson R": {"name": "RJ", "fwhm": ("590.0", "810.0")},
        "Johnson I": {"name": "IJ", "fwhm": ("780.0", "1020.0")},

        # Cousins
        "Cousins R": {"name": "R", "fwhm": ("561.7", "719.7")},
        "Cousins I": {"name": "I", "fwhm": ("721.0", "875.0")},

        # Near-Infrared
        "Near-Infrared J": {"name": "J", "fwhm": ("1040.0", "1360.0")},
        "Near-Infrared H": {"name": "H", "fwhm": ("1420.0", "1780.0")},
        "Near-Infrared K": {"name": "K", "fwhm": ("2015.0", "2385.0")},

        # Sloan
        "Sloan u": {"name": "SU", "fwhm": ("321.8", "386.8")},
        "Sloan g": {"name": "SG", "fwhm": ("402.5", "551.5")},
        "Sloan r": {"name": "SR", "fwhm": ("553.1", "693.1")},
        "Sloan i": {"name": "SI", "fwhm": ("697.5", "827.5")},
        "Sloan z": {"name": "SZ", "fwhm": ("841.2", "978.2")},

        # Stromgren
        "Stromgren u": {"name": "STU", "fwhm": ("336.3", "367.7")},
        "Stromgren v": {"name": "STV", "fwhm": ("401.5", "418.5")},
        "Stromgren b": {"name": "STB", "fwhm": ("459.55", "478.05")},
        "Stromgren y": {"name": "STY", "fwhm": ("536.7", "559.3")},
        "Stromgren Hbw": {"name": "STHBW", "fwhm": ("481.5", "496.5")},
        "Stromgren Hbn": {"name": "STHBN", "fwhm": ("487.5", "484.5")},

        # Optec
        "Optec Wing A": {"name": "MA", "fwhm": ("706.5", "717.5")},
        "Optec Wing B": {"name": "MB", "fwhm": ("748.5", "759.5")},
        "Optec Wing C": {"name": "MI", "fwhm": ("1003.0", "1045.0")},

        # PanSTARRS
        "PanSTARRS z-short": {"name": "ZS", "fwhm": ("826.0", "920.0")},
        "PanSTARRS Y": {"name": "Y", "fwhm": ("946.4", "1054.4")},
        "PanSTARRS w": {"name": "N/A", "fwhm": ("404.2", "845.8")},

        # MObs Clear Filter; Source(s): Martin Fowler
        "MObs CV": {"name": "CV", "fwhm": ("350.0", "850.0")},

        # Astrodon CBB; Source(s): George Silvis; https://astrodon.com/products/astrodon-exo-planet-filter/
        "Astrodon ExoPlanet-BB": {"name": "CBB", "fwhm": ("500.0", "1000.0")},
}
# expose as fwhm and for convenience set 'desc' field equal to key
fwhm = {k: v for k, v in __fwhm.items() if (v.update(desc=k),)}

# aliases to back-reference naming standard updates
fwhm_alias = {
        "LCO Bessell B": "Johnson B",
        "LCO Bessell V": "Johnson V",
        "LCO Bessell R": "Cousins R",
        "LCO Bessell I": "Cousins I",

        "J NIR 1.2 micron": "Near-Infrared J",
        "H NIR 1.6 micron": "Near-Infrared H",
        "K NIR 2.2 micron": "Near-Infrared K",

        "LCO SDSS u'": "Sloan u",
        "LCO SDSS g'": "Sloan g",
        "LCO SDSS r'": "Sloan r",
        "LCO SDSS i'": "Sloan i",

        "LCO Pan-STARRS zs": "PanSTARRS z-short",
        "LCO Pan-STARRS Y": "PanSTARRS Y",
        "LCO Pan-STARRS w": "PanSTARRS w",

        "Clear (unfiltered) reduced to V sequence": "MObs CV",
        "Clear (unfiltered) reduced to R sequence": "Cousins R",

        "Clear with blue-blocking": "Astrodon ExoPlanet-BB",
        "Exop": "Astrodon ExoPlanet-BB",
}

# standard filters w/o precisely defined FWHM values
fwhm_names_nonspecific = {
        'CR': "Clear (unfiltered) reduced to R sequence",
        'CBB': "Clear with blue-blocking",
        'CV': "Clear (unfiltered) reduced to V sequence",
        'TB': "DSLR Blue",
        'TG': "DSLR Green",
        'TR': "DSLR Red",
}
