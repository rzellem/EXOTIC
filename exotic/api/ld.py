import numpy as np

try:
    from .gael_ld import createldgrid
except ImportError:
    from gael_ld import createldgrid
try:
    from .filters import fwhm, fwhm_alias
except ImportError:
    from filters import fwhm, fwhm_alias


class LimbDarkening:
    ld0 = ld1 = ld2 = ld3 = None
    filter_desc = filter_type = None
    wl_min = wl_max = None
    fwhm = fwhm

    def __init__(self, stellar):
        self.priors = {
            'T*': stellar['teff'],
            'T*_uperr': stellar['teffUncPos'],
            'T*_lowerr': stellar['teffUncNeg'],
            'FEH*': stellar['met'],
            'FEH*_uperr': stellar['metUncPos'],
            'FEH*_lowerr': stellar['metUncNeg'],
            'LOGG*': stellar['logg'],
            'LOGG*_uperr': stellar['loggUncPos'],
            'LOGG*_lowerr': stellar['loggUncNeg']
        }

    def standard_list(self):
        print("\n\n***************************")
        print("Limb Darkening Coefficients")
        print("***************************")
        print("\nThe standard bands that are available for limb darkening parameters (https://www.aavso.org/filters)"
              "\nas well as filters for MObs and LCO (0.4m telescope) datasets:\n")
        for val in self.fwhm.values():
            match = next((f" (or {k})" for k, v in fwhm_alias.items() if v.lower() == val['desc'].lower()), '')
            print(f"\u2022 {val['desc']}{match}:\n\t-Abbreviation: {val['name']}"
                  f"\n\t-FWHM: ({val['fwhm'][0]}-{val['fwhm'][1]}) nm")

    def check_standard(self, filter_):
        if filter_['filter'] in fwhm_alias:
            filter_['filter'] = fwhm_alias[filter_['filter']]
        for val in self.fwhm.values():
            if filter_['filter'].lower() in (val['desc'].lower(), val['name'].lower()):
                self.filter_type = val['desc']
                self.set_filter(self.fwhm[self.filter_type]['name'], self.filter_type,
                                float(self.fwhm[self.filter_type]['fwhm'][0]),
                                float(self.fwhm[self.filter_type]['fwhm'][1]))
                return True
        else:
            return False

    def set_filter(self, filter_type, filter_desc, wl_min, wl_max):
        self.filter_type = filter_type
        self.filter_desc = filter_desc
        self.wl_min = wl_min
        self.wl_max = wl_max

    def calculate_ld(self):
        ld_params = createldgrid(np.array([self.wl_min / 1000]), np.array([self.wl_max / 1000]), self.priors)
        self.set_ld((ld_params['LD'][0][0], ld_params['ERR'][0][0]),
                    (ld_params['LD'][1][0], ld_params['ERR'][1][0]),
                    (ld_params['LD'][2][0], ld_params['ERR'][2][0]),
                    (ld_params['LD'][3][0], ld_params['ERR'][3][0]))

    def set_ld(self, ld0, ld1, ld2, ld3):
        self.ld0 = ld0
        self.ld1 = ld1
        self.ld2 = ld2
        self.ld3 = ld3

    def output_ld(self):
        print("\nEXOTIC-calculated nonlinear limb-darkening coefficients: ")
        print(f"{self.ld0[0]:5f} +/- + {self.ld0[1]:5f}")
        print(f"{self.ld1[0]:5f} +/- + {self.ld1[1]:5f}")
        print(f"{self.ld2[0]:5f} +/- + {self.ld2[1]:5f}")
        print(f"{self.ld3[0]:5f} +/- + {self.ld3[1]:5f}")


def test_ld(ld_obj_, filter_):
    try:
        ld_obj_.check_standard(filter_)
        ld_obj_.calculate_ld()
    except (KeyError, AttributeError):
        filter_['filter'] = "Custom"

        if filter_['wl_min'] and filter_['wl_max']:
            ld_obj_.set_filter('N/A', filter_['filter'], float(filter_['wl_min']), float(filter_['wl_max']))
            ld_obj_.calculate_ld()
        else:
            ld_ = [(filter_[key]["value"], filter_[key]["uncertainty"]) for key in filter_.keys()
                   if key in ['u0', 'u1', 'u2', 'u3']]
            ld_obj_.set_filter('N/A', filter_['filter'], filter_['wl_min'], filter_['wl_max'])
            ld_obj_.set_ld(ld_[0], ld_[1], ld_[2], ld_[3])


if __name__ == "__main__":
    stellar_params = {
        'teff': 6001.0,
        'teffUncPos': 88.0,
        'teffUncNeg': -88.0,
        'met': -0.16,
        'metUncPos': 0.08,
        'metUncNeg': -0.08,
        'logg': 4.22,
        'loggUncPos': 0.04,
        'loggUncNeg': -0.04
    }

    # Test existing filter
    filter_info1 = {
        'filter': "CV",
        'wl_min': None,
        'wl_max': None,
        'u0': {"value": None, "uncertainty": None},
        'u1': {"value": None, "uncertainty": None},
        'u2': {"value": None, "uncertainty": None},
        'u3': {"value": None, "uncertainty": None}
    }

    # Test alias filter
    filter_info2 = {
        'filter': "LCO SDSS u'",
        'wl_min': None,
        'wl_max': None,
        'u0': {"value": None, "uncertainty": None},
        'u1': {"value": None, "uncertainty": None},
        'u2': {"value": None, "uncertainty": None},
        'u3': {"value": None, "uncertainty": None}
    }

    # Test given only FWHM
    filter_info3 = {
        'filter': None,
        'wl_min': "350.0",
        'wl_max': "850.0",
        'u0': {"value": None, "uncertainty": None},
        'u1': {"value": None, "uncertainty": None},
        'u2': {"value": None, "uncertainty": None},
        'u3': {"value": None, "uncertainty": None}
    }

    # Test custom-entered ld coefficients
    filter_info4 = {
        'filter': None,
        'wl_min': None,
        'wl_max': None,
        'u0': {"value": 2.118, "uncertainty": 0.051},
        'u1': {"value": -3.88, "uncertainty": 0.21},
        'u2': {"value": 4.39, "uncertainty": 0.27},
        'u3': {"value": -1.63, "uncertainty": 0.12}
    }

    ld_obj = LimbDarkening(stellar_params)
    ld_obj.standard_list()

    test_ld(ld_obj, filter_info1)
    # test_ld(ld_obj, filter_info2)
    # test_ld(ld_obj, filter_info3)
    # test_ld(ld_obj, filter_info4)

    ld = [ld_obj.ld0[0], ld_obj.ld1[0], ld_obj.ld2[0], ld_obj.ld3[0]]
    ld_unc = [ld_obj.ld0[1], ld_obj.ld1[1], ld_obj.ld2[1], ld_obj.ld3[1]]
    ld_obj.output_ld()
