# Source for FWHM band wavelengths (units: nm): https://www.aavso.org/filters
# Near-Infrared
fwhm = {('J NIR 1.2micron', 'J'): (1040.00, 1360.00), ('H NIR 1.6micron', 'H'): (1420.00, 1780.00),
        ('K NIR 2.2micron', 'K'): (2015.00, 2385.00),

        # Sloan
        ('Sloan u', 'SU'): (321.80, 386.80), ('Sloan g', 'SG'): (402.50, 551.50),
        ('Sloan r', 'SR'): (553.10, 693.10), ('Sloan i', 'SI'): (697.50, 827.50),
        ('Sloan z', 'SZ'): (841.20, 978.20),

        # Stromgren
        ('Stromgren b', 'STB'): (459.55, 478.05), ('Stromgren y', 'STY'): (536.70, 559.30),
        ('Stromgren Hbw', 'STHBW'): (481.50, 496.50), ('Stromgren Hbn', 'STHBN'): (487.50, 484.50),

        # Johnson
        ('Johnson U', 'U'): (333.80, 398.80), ('Johnson B', 'B'): (391.60, 480.60),
        ('Johnson V', 'V'): (502.80, 586.80), ('Johnson R', 'RJ'): (590.00, 810.00),
        ('Johnson I', 'IJ'): (780.00, 1020.00),

        # Cousins
        ('Cousins R', 'R'): (561.70, 719.70), ('Cousins I', 'I'): (721.00, 875.00),

        # MObs Clear Filter, Source: Martin Fowler
        ('MObs CV', 'CV'): (350.00, 850.00),

        # Astrodon CBB: George Silvis: https://astrodon.com/products/astrodon-exo-planet-filter/
        ('Astrodon ExoPlanet-BB', 'CBB'): (500.00, 1000.00),

        # LCO, Source: Kalee Tock & Michael Fitzgerald, https://lco.global/observatory/instruments/filters/
        ('LCO Bessell B', 'N/A'): (391.60, 480.60), ('LCO Bessell V', 'N/A'): (502.80, 586.80),
        ('LCO Bessell U', 'N/A'): (325.0, 375.00), ('Harris B', 'N/A'): (377.75, 484.65),
        ('LCO Pan-STARRS w', 'N/A'): (404.20, 845.80), ('LCO Pan-STARRS w', 'N/A'): (404.20, 845.80),
        ('LCO Pan-STARRS zs', 'N/A'): (818.00, 922.00), ('LCO SDSS g\'', 'N/A'): (402.00, 552.00),
        ('LCO SDSS r\'', 'N/A'): (552.00, 691.00), ('LCO SDSS i\'', 'N/A'): (690.00, 819.00)}