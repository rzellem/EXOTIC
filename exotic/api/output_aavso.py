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
import hashlib
from json import dump, dumps
from numpy import mean, median, std
import os
from pathlib import Path
import re
from tkinter import NONE

try:
    from .utils import round_to_2
except ImportError:
    from utils import round_to_2
try:
    from .version import __version__
except ImportError:
    from version import __version__


class OutputFiles:
    def __init__(self, fit, p_dict, i_dict, planetdir):
        self.fit = fit
        self.p_dict = p_dict
        self.plname = p_dict['pl_name'] 
        self.p_dict['pl_name'] = re.sub(r"[-.]([a-kA-K]|\d{1,2})$", r" \1", self.p_dict['pl_name'], count=1)
        self.i_dict = i_dict
        self.durs = fit.duration_measured
        self.dir = Path(planetdir)
        
    def final_lightcurve(self, phase):
        params_file = self.dir / f"FinalLightCurve_{self.plname}_TESS.csv"

        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {self.p_dict['pl_name']}\n")
            f.write("# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n")

            for bjd, phase, flux, fluxerr, model, am in zip(self.fit.time, phase, self.fit.detrended,
                                                            self.fit.dataerr / self.fit.airmass_model,
                                                            self.fit.transit, self.fit.airmass_model):
                f.write(f"{bjd}, {phase}, {flux}, {fluxerr}, {model}, {am}\n")

    def final_planetary_params(self, phot_opt, comp_star=None, comp_coords=None, min_aper=None, min_annul=None):
        params_file = self.dir / f"FinalParams_{self.plname}_TESS.json"

        params_num = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- "
                                       f"{round_to_2(self.fit.errors['tmid'])} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/Rs)": f"{round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- "
                                                         f"{round_to_2(self.fit.errors['rprs'])}",
            "Transit depth (Rp/Rs)^2": f"{round_to_2(100. * (self.fit.parameters['rprs'] ** 2.))} +/- "
                                       f"{round_to_2(100. * 2. * self.fit.parameters['rprs'] * self.fit.errors['rprs'])} [%]",
            "Semi Major Axis/Star Radius (a/Rs)": f"{round_to_2(self.fit.parameters['ars'], self.fit.errors['ars'])} +/- "
                                                  f"{round_to_2(self.fit.errors['ars'])} ",
            "Airmass coefficient 1 (a1)": f"{round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- "
                                          f"{round_to_2(self.fit.errors['a1'])}",
            "Airmass coefficient 2 (a2)": f"{round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- "
                                          f"{round_to_2(self.fit.errors['a2'])}",
            "Scatter in the residuals of the lightcurve fit is": f"{round_to_2(100. * std(self.fit.residuals / median(self.fit.data)))} %",
        }

        if phot_opt:
            phot_ext = {"Best Comparison Star": f"#{comp_star} - {comp_coords}" if min_aper >= 0 else str(comp_star)}
            if min_aper == 0:
                phot_ext["Optimal Method"] = "PSF photometry"
            else:
                phot_ext["Optimal Aperture"] = f"{abs(min_aper)}"
                phot_ext["Optimal Annulus"] = f"{min_annul}"
            params_num.update(phot_ext)

        params_num["Transit Duration (day)"] = (f"{round_to_2(mean(self.durs), std(self.durs))} +/- "
                                                f"{round_to_2(std(self.durs))}")
        final_params = {'FINAL PLANETARY PARAMETERS': params_num}

        with params_file.open('w') as f:
            dump(final_params, f, indent=4)

    def aavso(self, airmasses, ld0, ld1, ld2, ld3, tmidstr):
        priors_dict, filter_dict, results_dict = aavso_dicts(self.p_dict, self.fit, self.i_dict, self.durs,
                                                             ld0, ld1, ld2, ld3)

        # compute 32 character hash of the results_dict
        hash_object = hashlib.sha256(dumps(results_dict).encode())
        hash_id = hash_object.hexdigest()[:32]

        #params_file = self.dir / f"TESS_{hash_id}_{self.plname}_{tmidstr}_AAVSO.txt"
        params_file = self.dir / f"{tmidstr}_{hash_id}_{self.plname}_AAVSO.txt"
        # 2459642_61_5164d266e1755aead98dbec0f26e7b7c_gj436b_AAVSO

        with params_file.open('w') as f:
            f.write("#TYPE=EXOPLANET\n"  # fixed
                    f"#OBSCODE=TESS\n"  # UI
                    f"#SECONDARY_OBSCODES=\n"  # UI 
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSTYPE=CCD\n"
                    f"#STAR_NAME={self.p_dict['hostname']}\n"  # code yields
                    f"#EXOPLANET_NAME={self.p_dict['pl_name']}\n"  # code yields
                    f"#BINNING=1x1\n"  # uhhh i just put One. 
                    f"#EXPOSURE_TIME={self.i_dict.get('exposure', -1)}\n"  # UI 
                    f"#COMP_STAR-XC=null\n"
                    f"#NOTES=TESS Data\n"
                    "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"  # fixed
                    "#MEASUREMENT_TYPE=Rnflux\n"  # fixed
                    f"#FILTER=I\n" 
                    f"#FILTER-XC={dumps(filter_dict)}\n"
                    f"#PRIORS=Period={round_to_2(self.p_dict['pl_orbper'], self.p_dict['pl_orbpererr1'])} +/- {round_to_2(self.p_dict['pl_orbpererr1'])}"
                    f",a/R*={round_to_2(self.p_dict['pl_ratdor'], self.p_dict['pl_ratdorerr1'])} +/- {round_to_2(self.p_dict['pl_ratdorerr1'])}"
                    f",inc={round_to_2(self.p_dict['pl_orbincl'], self.p_dict['pl_orbinclerr1'])} +/- {round_to_2(self.p_dict['pl_orbinclerr1'])}"
                    f",ecc={round_to_2(self.p_dict['pl_orbeccen'])}"
                    f",u0={round_to_2(ld0)}"
                    f",u1={round_to_2(ld1)}"
                    f",u2={round_to_2(ld2)}"
                    f",u3={round_to_2(ld3)}\n"
                    f"#PRIORS-XC={dumps(priors_dict)}\n"  # code yields
                    f"#RESULTS=Tc={round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- {round_to_2(self.fit.errors['tmid'])}"
                    f",Rp/R*={round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- {round_to_2(self.fit.errors['rprs'])}"
                    f",Am1=0"
                    f",Am2=0\n"
                    f"#RESULTS-XC={dumps(results_dict)}\n")  # code yields

            f.write(
                "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                "This work is supported by NASA under award number NNX16AC65A to the "
                "Space Telescope Science Institute.\n"
                "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                "aavso.org/data-usage-guidelines\n")

            f.write("#DATE,DIFF,ERR,DETREND_1,DETREND_2\n")
            for aavsoC in range(0, len(self.fit.time)):
                # f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC] / self.fit.parameters['a1'], 7)},"
                #         f"{round(self.fit.dataerr[aavsoC] / self.fit.parameters['a1'], 7)},{round(airmasses[aavsoC], 7)},"
                #         f"{round(self.fit.airmass_model[aavsoC] / self.fit.parameters['a1'], 7)}\n")
                f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC], 7)},"
                        f"{round(self.fit.dataerr[aavsoC], 7)},{0.0},"
                        f"{1.000}\n")

# *********** FOR CSV ************
    def aavso_csv(self, airmasses, ld0, ld1, ld2, ld3,tmidstr):
            priors_dict, filter_dict, results_dict = aavso_dicts(self.p_dict, self.fit, self.i_dict, self.durs,
                                                                ld0, ld1, ld2, ld3)

            params_file = self.dir / f"TESS_{tmidstr}_{self.p_dict['pl_name']}_lightcurve.csv"

            with params_file.open('w') as f:
                f.write("#TYPE=EXOPLANET\n"  # fixed
                        f"#OBSCODE=TESS\n"  # UI
                        f"#SECONDARY_OBSCODES=\n"  # UI 
                        f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                        "#DELIM=,\n"  # fixed
                        "#DATE_TYPE=BJD_TDB\n"  # fixed
                        f"#OBSTYPE=CCD\n"
                        f"#STAR_NAME={self.p_dict['hostname']}\n"  # code yields
                        f"#EXOPLANET_NAME={self.p_dict['pl_name']}\n"  # code yields
                        f"#BINNING=1x1\n"  # uhhh i just put One. 
                        f"#EXPOSURE_TIME={self.i_dict.get('exposure', -1)}\n"  # UI 
                        f"#COMP_STAR-XC=null\n"
                        f"#NOTES=TESS Data\n"
                        "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"  # fixed
                        "#MEASUREMENT_TYPE=Rnflux\n"  # fixed
                        f"#FILTER=I\n" 
                        f"#FILTER-XC={dumps(filter_dict)}\n"
                        f"#PRIORS=Period={round_to_2(self.p_dict['pl_orbper'], self.p_dict['pl_orbpererr1'])} +/- {round_to_2(self.p_dict['pl_orbpererr1'])}"
                        f",a/R*={round_to_2(self.p_dict['pl_ratdor'], self.p_dict['pl_ratdorerr1'])} +/- {round_to_2(self.p_dict['pl_ratdorerr1'])}"
                        f",inc={round_to_2(self.p_dict['pl_orbincl'], self.p_dict['pl_orbinclerr1'])} +/- {round_to_2(self.p_dict['pl_orbinclerr1'])}"
                        f",ecc={round_to_2(self.p_dict['pl_orbeccen'])}"
                        f",u0={round_to_2(ld0)}"
                        f",u1={round_to_2(ld1)}"
                        f",u2={round_to_2(ld2)}"
                        f",u3={round_to_2(ld3)}\n"
                        f"#PRIORS-XC={dumps(priors_dict)}\n"  # code yields
                        f"#RESULTS=Tc={round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- {round_to_2(self.fit.errors['tmid'])}"
                        f",Rp/R*={round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- {round_to_2(self.fit.errors['rprs'])}"
                        f",Am1=0"
                        f",Am2=0\n"
                        f"#RESULTS-XC={dumps(results_dict)}\n")  # code yields

                f.write(
                    "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                    "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                    "This work is supported by NASA under award number NNX16AC65A to the "
                    "Space Telescope Science Institute.\n"
                    "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                    "aavso.org/data-usage-guidelines\n")

                f.write("#DATE,DIFF,ERR,DETREND_1,DETREND_2\n")
                for aavsoC in range(0, len(self.fit.time)):
                    # f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC] / self.fit.parameters['a1'], 7)},"
                    #         f"{round(self.fit.dataerr[aavsoC] / self.fit.parameters['a1'], 7)},{round(airmasses[aavsoC], 7)},"
                    #         f"{round(self.fit.airmass_model[aavsoC] / self.fit.parameters['a1'], 7)}\n")
                    f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC], 7)},"
                            f"{round(self.fit.dataerr[aavsoC], 7)},{round(airmasses[aavsoC], 7)},"
                            f"{round(self.fit.airmass_model[aavsoC], 7)}\n")


def aavso_dicts(planet_dict, fit, i_dict, durs, ld0, ld1, ld2, ld3):
    priors = {
        'Period': {
            'value': str(round_to_2(planet_dict['pl_orbper'], planet_dict['pl_orbpererr1'])),
            'uncertainty': str(round_to_2(planet_dict['pl_orbpererr1'])) if planet_dict['pl_orbpererr1'] else None,
            'units': "days"
        },
        'a/R*': {
            'value': str(round_to_2(planet_dict['pl_ratdor'], planet_dict['pl_ratdorerr1'])),
            'uncertainty': str(round_to_2(planet_dict['pl_ratdorerr1'])) if planet_dict['pl_ratdorerr1'] else planet_dict['pl_ratdorerr1'],
        },
        'inc': {
            'value': str(round_to_2(planet_dict['pl_orbincl'], planet_dict['pl_orbinclerr1'])),
            'uncertainty': str(round_to_2(planet_dict['pl_orbinclerr1'])) if planet_dict['pl_orbinclerr1'] else planet_dict['pl_orbinclerr1'],
            'units': "degrees"
        },
        'ecc': {
            'value': str(round_to_2(planet_dict['pl_orbeccen'])),
            'uncertainty': None,
        },
        'u0': {
            'value': str(round_to_2(ld0)),
            'uncertainty': None
        },
        'u1': {
            'value': str(round_to_2(ld1)),
            'uncertainty': None
        },
        'u2': {
            'value': str(round_to_2(ld2)),
            'uncertainty': None
        },
        'u3': {
            'value': str(round_to_2(ld3)),
            'uncertainty': None
        }
    }

    filter_type = {
        'name': "I",
        'fwhm': [{'value': 600, 'units': "nm"},
                {'value': 1000, 'units': "nm"}],
    }

    results = {
        'Tc': {
            'value': str(round_to_2(fit.parameters['tmid'], fit.errors['tmid'])),
            'uncertainty': str(round_to_2(fit.errors['tmid'])),
            'units': "BJD_TDB"
        },
        'Rp/R*': {
            'value': str(round_to_2(fit.parameters['rprs'], fit.errors['rprs'])),
            'uncertainty': str(round_to_2(fit.errors['rprs']))
        },
        'Am1': {
            'value': 0,
            'uncertainty': None
        },
        'Am2': {
            'value': 0,
            'uncertainty': None
        },
        'Duration': {
            'value': str(round_to_2(mean(durs))),
            'uncertainty': str(round_to_2(std(durs))),
            'units': "days"
        }
    }

    # try to add a/Rs if it exists
    if 'ars' in fit.errors:
        results['a/R*'] = {
            'value': str(round_to_2(fit.parameters['ars'], fit.errors['ars'])),
            'uncertainty': str(round_to_2(fit.errors['ars'])),
        }

    # check for inclination
    if 'inc' in fit.errors:
        results['inc'] = {
            'value': str(round_to_2(fit.parameters['inc'], fit.errors['inc'])),
            'uncertainty': str(round_to_2(fit.errors['inc'])),
            'units': "degrees"
        }

    return priors, filter_type, results
