from numpy import pi
import astropy.constants as const
import astropy.units as u

AU = const.au                                                       # m
PI = pi
R_SUN = const.R_sun                                                 # m
R_JUP = const.R_jup                                                 # m

# CALCULATED VALUES
G = const.G.to(AU**3 / (const.M_sun * u.day**2))                    # AU^3 /(msun * day^2)
SA = lambda m, p: (G * m * p ** 2. / (4. * PI ** 2.)) ** (1. / 3.)  # Keplerian semi-major axis (au)
