# Sources for FWHM band wavelengths are referenced below, respectively 
# note that all below units default to nm
fwhm = {
        # AAVSO, Source(s): AAVSO International Database; https://www.aavso.org/filters
        # Johnson
        "Johnson U": {"name": "U", "desc": "Johnson U", "fwhm": ("333.8", "398.8")},
        "Johnson B": {"name": "B", "desc": "Johnson B", "fwhm": ("391.6", "480.6")},
        "Johnson V": {"name": "V", "desc": "Johnson V", "fwhm": ("502.8", "586.8")},
        "Johnson R": {"name": "RJ", "desc": "Johnson R", "fwhm": ("590.0", "810.0")},
        "Johnson I": {"name": "IJ", "desc": "Johnson I", "fwhm": ("780.0", "1020.0")},

        # Cousins
        "Cousins R": {"name": "R", "desc": "Cousins R", "fwhm": ("561.7", "719.7")},
        "Cousins I": {"name": "I", "desc": "Cousins I", "fwhm": ("721.0", "875.0")},

        # Near-Infrared
        "Near-Infrared J": {"name": "J", "desc": "Near-Infrared J", "fwhm": ("1040.0", "1360.0")},
        "Near-Infrared H": {"name": "H", "desc": "Near-Infrared H", "fwhm": ("1420.0", "1780.0")},
        "Near-Infrared K": {"name": "K", "desc": "Near-Infrared K", "fwhm": ("2015.0", "2385.0")},

        # Sloan
        "Sloan u": {"name": "SU", "desc": "Sloan u", "fwhm": ("321.8", "386.8")},
        "Sloan g": {"name": "SG", "desc": "Sloan g", "fwhm": ("402.5", "551.5")},
        "Sloan r": {"name": "SR", "desc": "Sloan r", "fwhm": ("553.1", "693.1")},
        "Sloan i": {"name": "SI", "desc": "Sloan i", "fwhm": ("697.5", "827.5")},
        "Sloan z": {"name": "SZ", "desc": "Sloan z", "fwhm": ("841.2", "978.2")},

        # Stromgren
        "Stromgren u": {"name": "STU", "desc": "Stromgren u", "fwhm": ("336.3", "367.7")},
        "Stromgren v": {"name": "STV", "desc": "Stromgren v", "fwhm": ("401.5", "418.5")},
        "Stromgren b": {"name": "STB", "desc": "Stromgren b", "fwhm": ("459.55", "478.05")},
        "Stromgren y": {"name": "STY", "desc": "Stromgren y", "fwhm": ("536.7", "559.3")},
        "Stromgren Hbw": {"name": "STHBW", "desc": "Stromgren Hbw", "fwhm": ("481.5", "496.5")},
        "Stromgren Hbn": {"name": "STHBN", "desc": "Stromgren Hbn", "fwhm": ("487.5", "484.5")},

        # Optec
        "Optec Wing A": {"name": "MA", "desc": "Optec Wing A", "fwhm": ("706.5", "717.5")},
        "Optec Wing B": {"name": "MB", "desc": "Optec Wing B", "fwhm": ("748.5", "759.5")},
        "Optec Wing C": {"name": "MI", "desc": "Optec Wing C", "fwhm": ("1003.0", "1045.0")},

        # PanSTARRS
        "PanSTARRS z-short": {"name": "ZS", "desc": "PanSTARRS z-short", "fwhm": ("826.0", "920.0")},
        "PanSTARRS Y": {"name": "Y", "desc": "PanSTARRS Y", "fwhm": ("946.4", "1054.4")},

        # MObs Clear Filter; Source(s): Martin Fowler
        "MObs CV": {"name": "CV", "desc": "MObs CV", "fwhm": ("350.0", "850.0")},

        # Astrodon CBB; Source(s): George Silvis; https://astrodon.com/products/astrodon-exo-planet-filter/
        "Astrodon ExoPlanet-BB": {"name": "CBB", "desc": "Astrodon ExoPlanet-BB", "fwhm": ("500.0", "1000.0")},
}

# aliases to back-reference naming standard updates
fwhm_alias = {
        "J NIR 1.2 micron": "Near-Infrared J",
        "J NIR 1.2micron": "Near-Infrared J",
        "H NIR 1.6 micron": "Near-Infrared H",
        "H NIR 1.6micron": "Near-Infrared H",
        "K NIR 2.2 micron": "Near-Infrared K",
        "K NIR 2.2micron": "Near-Infrared K",

        "LCO Bessell B": "Johnson B",
        "LCO Bessell V": "Johnson V",
        "LCO Bessell R": "Cousins R",
        "LCO Bessell I": "Cousins I",

        "LCO SDSS u'": "Sloan u",
        "LCO SDSS g'": "Sloan g",
        "LCO SDSS r'": "Sloan r",
        "LCO SDSS i'": "Sloan i",

        "LCO Pan-STARRS Y": "PanSTARRS Y",
        "LCO Pan-STARRS zs": "PanSTARRS z-short",

        "Exop": "Astrodon ExoPlanet-BB",

        "Clear (unfiltered) reduced to V sequence": "MObs CV",
}
