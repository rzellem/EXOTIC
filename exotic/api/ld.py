# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #
import logging
import numpy as np
import re

try:
    from gael_ld import createldgrid
except ImportError:
    from .gael_ld import createldgrid
try:
    from filters import fwhm, fwhm_alias, fwhm_names_nonspecific
except ImportError:
    from .filters import fwhm, fwhm_alias, fwhm_names_nonspecific

log = logging.getLogger(__name__)
ld_re_punct = r'[\W]'
ld_re_punct_p = re.compile(ld_re_punct)


class LimbDarkening:
    # implies custom filter, always
    filter_names_undefined = {'N/A', 'NA', 'NONE', }

    # include filter maps as references from this class
    fwhm = fwhm
    fwhm_alias = fwhm_alias
    fwhm_names_nonspecific = fwhm_names_nonspecific

    # lookup table: fwhm_lookup references filters irrespective of spacing and punctuation
    # 1 - combine optimized str lookups in lookup table
    fwhm_lookup = {k.strip().replace(' ', '').lower(): k for k in fwhm.keys()}
    fwhm_lookup.update({k.strip().replace(' ', '').lower(): v for k, v in fwhm_alias.items()})
    # 2 - ignore punctuation in lookup table
    fwhm_lookup = {re.sub(ld_re_punct_p, '', k): v for k, v in fwhm_lookup.items()}
    # lookup set: filter_desc_nonspecific_lookup_set references descriptions that do not represent a specific filter
    #     irrespective of spacing and punctuation
    # 1 - combine optimized str lookups in lookup set
    # 2 - ignore punctuation in lookup set
    filter_desc_nonspecific_lookup_set = tuple([re.sub(ld_re_punct_p, '', k.strip().replace(' ', '')).lower()
                                                for k in fwhm_names_nonspecific.values()])

    def __init__(self, stellar):
        self.filter_name = self.filter_desc = None
        self.ld0 = self.ld1 = self.ld2 = self.ld3 = [0, 0]
        self.priors = {
            'T*': stellar.get('teff'),
            'T*_uperr': stellar.get('teffUncPos'),
            'T*_lowerr': stellar.get('teffUncNeg'),
            'FEH*': stellar.get('met'),
            'FEH*_uperr': stellar.get('metUncPos'),
            'FEH*_lowerr': stellar.get('metUncNeg'),
            'LOGG*': stellar.get('logg'),
            'LOGG*_uperr': stellar.get('loggUncPos'),
            'LOGG*_lowerr': stellar.get('loggUncNeg')
        }
        self.wl_max = self.wl_min = None
        super()

    @staticmethod
    def standard_list():
        print("\n\n***************************")
        print("Limb Darkening Coefficients")
        print("***************************\n")
        print("The standard bands that are available for limb darkening parameters (https://www.aavso.org/filters)\n"
              "as well as filters for MObs and LCO (0.4m telescope) datasets:\n")
        for val in LimbDarkening.fwhm.values():
            match = next((f" (or {k})" for k, v in LimbDarkening.fwhm_alias.items() if
                          v.lower() == val['desc'].lower()), '')
            print(f"\u2022 {val['desc']}{match}:\n\t-Abbreviation: {val['name']}"
                  f"\n\t-FWHM: ({val['fwhm'][0]}-{val['fwhm'][1]}) nm")
        return

    def check_standard(self, filter_: dict = None, loose: bool = False, loose_len: int = 5) -> bool:
        """
        Utility method to detect filter from full or partial values and to inject
            results into LimbDarkening object. Detection algorithm is loose and
            finds alias from unique partial matches on the `filter` key. The input
            dict must contain one or more of the required keys to facilitate
            matching.
        Order of matching operations proceeds as follows:
            1 - 'filter': loose match against name (irrespective of spacing, caps, punct)
            2 - Both supplied min/max wavelength values: precise match
            3 - 'filter': if it precisely matches any filter key (same as 'desc')
            4 - Abbrev. name if it is in the 'filter' or 'name' fields: precise match
            5 - Any portion of filter name that is unique to one filter
        @param filter_: A dictionary containing any combination of ('filter', 'name',
            'wl_min', 'wl_max') keys implying a certain visibility filter used
             during telescope observations.
        @type filter_: dict
        @param loose: (default False) Flag to use loose failover matching, if all else fails.
        @type loose: bool
        @param loose_len: (default 5) Length of optimized matcher to use for loose matches.
        @type loose_len: int
        @return: True if filter parameters match known aliases and false if any do not.
        @rtype: bool
        """
        filter_alias = filter_matcher = None
        try:
            if not isinstance(filter_, dict):  # set sane failure
                log.error("Filter not defined according to spec -- parsing failure.")
                return False
            for k in ('filter', 'name', 'wl_min', 'wl_max'):  # set missing values to nonetype
                filter_[k] = filter_.get(k)
                if isinstance(filter_[k], str):
                    # eliminate errant spaces on edges
                    filter_[k] = filter_[k].strip()
                    if k == 'name' and filter_[k]:  # format 'name' (if exists) to uppercase, no spaces
                        filter_[k] = filter_[k].upper().replace(' ', '')
            if filter_['filter']:  # make matcher by removing spaces, remove punctuation and lowercase
                filter_matcher = filter_['filter'].lower().replace(' ', '')
                filter_matcher = re.sub(ld_re_punct_p, '', filter_matcher)
            # names that do not represent a specific filter combined into one tuple
            filter_names_nonspecific = set(LimbDarkening.fwhm_names_nonspecific.keys())
            filter_names_nonspecific.update(LimbDarkening.filter_names_undefined)
            # identify defined filters via optimized lookup table
            if (filter_matcher and filter_matcher in LimbDarkening.fwhm_lookup and
                    filter_matcher not in LimbDarkening.filter_desc_nonspecific_lookup_set):
                filter_['filter'] = LimbDarkening.fwhm_lookup[filter_matcher]  # sets to actual filter reference key
            for f in LimbDarkening.fwhm.values():
                # match to wavelength values (strict)
                if (filter_['wl_min'] and filter_['wl_min'] == f['fwhm'][0] and
                        filter_['wl_max'] and filter_['wl_max'] == f['fwhm'][1]):
                    filter_alias = f
                    break
                # match 'filter' (strict) to actual reference key, e.g. 'desc'
                elif filter_['filter'] == f['desc']:
                    filter_alias = f
                    break
                # match 'name' or 'filter' (strict) to actual name abbreviation, e.g. 'name'
                elif filter_['name'] == f['name'].strip().upper() or \
                        (filter_['filter'] and filter_['filter'][:5].upper() == f['name'].strip().upper()):
                    # exclude unknown vals for 'name'
                    if filter_['name'] in filter_names_nonspecific or \
                            (filter_['filter'] and filter_['filter'][:5].upper() in filter_names_nonspecific):
                        pass
                    else:
                        filter_alias = f
                        break
            # match 'filter' (loose) to any portion of actual reference key, e.g. 'desc' -- if possible
            if loose and (filter_matcher and len(filter_matcher) > loose_len) \
                    and not filter_alias:  # filter not identified so try another algorithm
                f_count = 0
                for f in LimbDarkening.fwhm_lookup:
                    if filter_matcher in f:
                        f_count += 1
                        filter_alias = LimbDarkening.fwhm[LimbDarkening.fwhm_lookup[f]]
                    if f_count > 1:  # no unique match, reset
                        filter_alias = None
                        break
            if filter_alias:
                filter_['name'] = filter_alias['name']
                filter_['filter'] = filter_alias['desc']
                filter_['wl_min'] = filter_alias['fwhm'][0]
                filter_['wl_max'] = filter_alias['fwhm'][1]
                self.set_filter(filter_['name'], filter_['filter'], float(filter_['wl_min']), float(filter_['wl_max']))
        except BaseException as be:
            log.error(f"Filter matching failed -- {be}")
        return filter_alias is not None

    @staticmethod
    def check_fwhm(filter_: dict = None) -> bool:
        """
        Validates wavelength values in a filter dict to verify they are in the correct
            order (e.g. min, max properly ordered) and are not out of bounds, between
            numeric values of 300nm and 4000nm. This mutates min/max value strings to
            ensure they end with '.0' and are ordered. NOTE: ASSUMES NANOMETER INPUTS.
        @param filter_: A dictionary containing full-width, half-maximum (fwhm)
            wavelength values ('wl_min' and 'wl_max') keys implying a certain visibility
            filter used during telescope observations.
        @type filter_: dict
        @return: True if wavelength parameters meet requirements and false if they do not.
            Note that the input dict is modified to correct malformed values, including
            popping detected bad keys.
        @rtype: bool
        """
        fwhm_tuple = None
        try:
            if not isinstance(filter_, dict):  # set sane failure
                log.error("Filter not defined according to spec (dict required) -- parsing failure.")
                return False
            for k in ('wl_min', 'wl_max'):  # clean inputs
                filter_[k] = filter_.get(k)
                filter_[k] = str(filter_[k]).strip().replace(' ', '').rstrip('.') if filter_[k] else filter_[k]
                if not 200. <= float(filter_[k]) <= 4000.:  # also fails if nan
                    raise ValueError(f"FWHM '{k}' is outside of bounds (200., 4000.). ...")
                else:  # add .0 to end of str to aid literal matching
                    if not filter_[k].count('.'):
                        filter_[k] += '.0'
            if float(filter_['wl_min']) > float(filter_['wl_max']):  # reorder if min > max
                fwhm_tuple = (filter_['wl_max'], filter_['wl_min'])
                filter_['wl_min'] = fwhm_tuple[0]
                filter_['wl_max'] = fwhm_tuple[1]
            else:
                fwhm_tuple = (filter_['wl_min'], filter_['wl_max'])
        except BaseException as be:
            log.error(f"FWHM matching failed [{filter_.get('wl_min')}, {filter_.get('wl_max')}]-- {be}")
        return fwhm_tuple is not None

    def set_filter(self, filter_name, filter_desc, wl_min, wl_max):
        self.filter_name = filter_name
        self.filter_desc = filter_desc
        self.wl_min = wl_min
        self.wl_max = wl_max
        return

    def calculate_ld(self):
        ld_params = createldgrid(np.array([self.wl_min / 1000.]), np.array([self.wl_max / 1000.]), self.priors)
        self.set_ld((ld_params['LD'][0][0], ld_params['ERR'][0][0]),
                    (ld_params['LD'][1][0], ld_params['ERR'][1][0]),
                    (ld_params['LD'][2][0], ld_params['ERR'][2][0]),
                    (ld_params['LD'][3][0], ld_params['ERR'][3][0]))
        return

    def set_ld(self, ld0, ld1, ld2, ld3):
        self.ld0 = ld0
        self.ld1 = ld1
        self.ld2 = ld2
        self.ld3 = ld3
        return

    def __str__(self):
        return (f"\nFilter Name: {self.filter_name}"
                f"\nFilter Description: {self.filter_desc}"
                f"\nMinimum Wavelength (nm): {self.wl_min}"
                f"\nMaxiumum Wavelength (nm):: {self.wl_max}"
                "\nEXOTIC-calculated nonlinear limb-darkening coefficients: "
                f"\n\t- {self.ld0[0]:5f} +/- + {self.ld0[1]:5f}"
                f"\n\t- {self.ld1[0]:5f} +/- + {self.ld1[1]:5f}"
                f"\n\t- {self.ld2[0]:5f} +/- + {self.ld2[1]:5f}"
                f"\n\t- {self.ld3[0]:5f} +/- + {self.ld3[1]:5f}")
