from numpy import pi
import astropy.constants as const
import astropy.units as u

PI = pi
AU = const.au                                       # m
R_SUN = const.R_sun                                 # m
R_JUP = const.R_jup                                 # m
G = const.G.to(AU**3 / (const.M_sun*u.day**2))      # AU^3 /(msun * day^2)

# SHARED LAMBDAS
# Keplerian semi-major axis (au)
SA = lambda m, P: (G * m * P ** 2 / (4 * PI ** 2)) ** (1. / 3)
