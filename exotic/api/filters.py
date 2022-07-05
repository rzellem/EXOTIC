# Source for FWHM band wavelengths (units: nm): https://www.aavso.org/filters
# Near-Infrared
fwhm = {
        "J NIR 1.2 micron": {"name": "J", "desc": "J NIR 1.2micron", "fwhm": ("1040.0", "1360.0")},
        "H NIR 1.6 micron": {"name": "H", "desc": "H NIR 1.6micron", "fwhm": ("1420.0", "1780.0")},
        "K NIR 2.2 micron": {"name": "K", "desc": "K NIR 2.2micron", "fwhm": ("2015.0", "2385.0")},

        # Sloan
        "Sloan u": {"name": "SU", "desc": "Sloan u", "fwhm": ("321.8", "386.8")},
        "Sloan g": {"name": "SG", "desc": "Sloan g", "fwhm": ("402.5", "551.5")},
        "Sloan r": {"name": "SR", "desc": "Sloan r", "fwhm": ("553.1", "693.1")},
        "Sloan i": {"name": "SI", "desc": "Sloan i", "fwhm": ("697.5", "827.5")},
        "Sloan z": {"name": "SZ", "desc": "Sloan z", "fwhm": ("841.2", "978.2")},

        # Stromgren
        "Stromgren b": {"name": "STB", "desc": "Stromgren b", "fwhm": ("459.55", "478.05")},
        "Stromgren y": {"name": "STY", "desc": "Stromgren y", "fwhm": ("536.7", "559.3")},
        "Stromgren Hbw": {"name": "STHBW", "desc": "Stromgren Hbw", "fwhm": ("481.5", "496.5")},
        "Stromgren Hbn": {"name": "STHBN", "desc": "Stromgren Hbn", "fwhm": ("487.5", "484.5")},

        # Johnson
        "Johnson U": {"name": "U", "desc": "Johnson U", "fwhm": ("333.8", "398.8")},
        "Johnson B": {"name": "B", "desc": "Johnson B", "fwhm": ("391.6", "480.6")},
        "Johnson V": {"name": "V", "desc": "Johnson V", "fwhm": ("502.8", "586.8")},
        "Johnson R": {"name": "RJ", "desc": "Johnson R", "fwhm": ("590.0", "810.0")},
        "Johnson I": {"name": "IJ", "desc": "Johnson I", "fwhm": ("780.0", "1020.0")},

        # Cousins
        "Cousins R": {"name": "R", "desc": "Cousins R", "fwhm": ("561.7", "719.7")},
        "Cousins I": {"name": "I", "desc": "Cousins I", "fwhm": ("721.0", "875.0")},

        # MObs Clear Filter, Source: Martin Fowler
        "MObs CV": {"name": "CV", "desc": "MObs CV", "fwhm": ("350.0", "850.0")},

        # Astrodon CBB: George Silvis: https://astrodon.com/products/astrodon-exo-planet-filter/
        "Astrodon ExoPlanet-BB": {"name": "CBB", "desc": "Astrodon ExoPlanet-BB", "fwhm": ("500.0", "1000.0")},

        # LCO, Source: Kalee Tock & Michael Fitzgerald, https://lco.global/observatory/instruments/filters/
        "LCO Bessell B": {"name": "N/A", "desc": "LCO Bessell B", "fwhm": ("391.6", "480.6")},
        "LCO Bessell V": {"name": "N/A", "desc": "LCO Bessell V", "fwhm": ("502.8", "586.8")},
        "LCO Pan-STARRS w": {"name": "N/A", "desc": "LCO Pan-STARRS w", "fwhm": ("404.2", "845.8")},
        "LCO Pan-STARRS zs": {"name": "N/A", "desc": "LCO Pan-STARRS zs", "fwhm": ("818.0", "922.0")},
        "LCO SDSS g'": {"name": "N/A", "desc": "LCO SDSS g'", "fwhm": ("402.0", "552.0")},
        "LCO SDSS r'": {"name": "N/A", "desc": "LCO SDSS r'", "fwhm": ("552.0", "691.0")},
        "LCO SDSS i'": {"name": "N/A", "desc": "LCO SDSS i'", "fwhm": ("690.0", "819.0")}
}
