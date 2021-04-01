class NEA:

    def __init__(self, init_file, init_opt):
        self.init_file = init_file
        self.init_opt = init_opt
        self.planet_dict = {
            'ra': None, 'dec': None, 'p_name': None, 's_name': None, 'per': None, 'per_unc': None,
            'mid_t': None, 'mid_t_unc': None, 'rprs': None, 'rprs_unc': None,
            'ars': None, 'ars_unc': None, 'inc': None, 'inc_unc': None, 'ecc': None,
            'teff': None, 'teff_unc_pos': None, 'teff_unc_neg': None,
            'met': None, 'met_unc_pos': None, 'met_unc_neg': None,
            'logg': None, 'logg_unc_pos': None, 'logg_unc_neg': None
        }
        comp_planet = {
            'ra': 'Target Star RA', 'dec': 'Target Star Dec', 'p_name': "Planet Name", 's_name': "Host Star Name",
            'per': 'Orbital Period (days)', 'per_unc': 'Orbital Period Uncertainty',
            'mid_t': 'Published Mid-Transit Time (BJD-UTC)', 'mid_t_unc': 'Mid-Transit Time Uncertainty',
            'rprs': 'Ratio of Planet to Stellar Radius (Rp/Rs)',
            'rprs_unc': 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
            'ars': 'Ratio of Distance to Stellar Radius (a/Rs)',
            'ars_unc': 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
            'inc': 'Orbital Inclination (deg)', 'inc_unc': 'Orbital Inclination (deg) Uncertainity',
            'ecc': 'Orbital Eccentricity (0 if null)', 'teff': 'Star Effective Temperature (K)',
            'teff_unc_pos': 'Star Effective Temperature (+) Uncertainty',
            'teff_unc_neg': 'Star Effective Temperature (-) Uncertainty', 'met': 'Star Metallicity ([FE/H])',
            'met_unc_pos': 'Star Metallicity (+) Uncertainty', 'met_unc_neg': 'Star Metallicity (-) Uncertainty',
            'logg': 'Star Surface Gravity (log(g))',
            'logg_unc_pos': 'Star Surface Gravity (+) Uncertainty',
            'logg_unc_neg': 'Star Surface Gravity (-) Uncertainty'
        }

        self.planet_dict = init_params(comp_planet, self.planet_dict, data['planetary_parameters'])
