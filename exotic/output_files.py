from json import dump, dumps
from numpy import mean, median, std
from pathlib import Path

try:
    from utils import round_to_2
except ImportError:
    from .utils import round_to_2
try:
    from version import __version__
except ImportError:
    from .version import __version__


class OutputFiles:
    def __init__(self, fit, p_dict, i_dict, durs):
        self.fit = fit
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.durs = durs
        self.dir = Path(self.i_dict['save'])

    def final_lightcurve(self, phase):
        params_file = self.dir / f"FinalLightCurve_{self.p_dict['pName']}_{self.i_dict['date']}.csv"

        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {self.p_dict['pName']}\n")
            f.write("# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n")

            for bjd, phase, flux, fluxerr, model, am in zip(self.fit.time, phase, self.fit.detrended,
                                                            self.fit.dataerr / self.fit.airmass_model,
                                                            self.fit.transit, self.fit.airmass_model):
                f.write(f"{bjd}, {phase}, {flux}, {fluxerr}, {model}, {am}\n")

    def final_planetary_params(self, phot_opt, comp_star=None, comp_coords=None, min_aper=None, min_annul=None):
        params_file = self.dir / f"FinalParams_{self.p_dict['pName']}_{self.i_dict['date']}.json"

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

        params_num["Transit Duration (day)"] = f"{round_to_2(mean(self.durs))} +/- {round_to_2(std(self.durs))}"
        final_params = {'FINAL PLANETARY PARAMETERS': params_num}

        with params_file.open('w') as f:
            dump(final_params, f, indent=4)

    def aavso(self, comp_star, airmasses, ld0, ld1, ld2, ld3):
        priors_dict, filter_dict, results_dict = aavso_dicts(self.p_dict, self.fit, self.i_dict, self.durs,
                                                             ld0, ld1, ld2, ld3)

        params_file = self.dir / f"AAVSO_{self.p_dict['pName']}_{self.i_dict['date']}.txt"

        with params_file.open('w') as f:
            f.write("#TYPE=EXOPLANET\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SECONDARY_OBSCODES={self.i_dict['second_obs']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSTYPE={self.i_dict['camera']}\n"
                    f"#STAR_NAME={self.p_dict['sName']}\n"  # code yields
                    f"#EXOPLANET_NAME={self.p_dict['pName']}\n"  # code yields
                    f"#BINNING={self.i_dict['pixel_bin']}\n"  # user input
                    f"#EXPOSURE_TIME={self.i_dict.get('exposure', -1)}\n"  # UI
                    f"#COMP_STAR-XC={dumps(comp_star)}\n"
                    f"#NOTES={self.i_dict['notes']}\n"
                    "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"  # fixed
                    "#MEASUREMENT_TYPE=Rnflux\n"  # fixed
                    f"#FILTER={self.i_dict['filter']}\n"
                    f"#FILTER-XC={dumps(filter_dict)}\n"
                    f"#PRIORS=Period={round_to_2(self.p_dict['pPer'], self.p_dict['pPerUnc'])} +/- {round_to_2(self.p_dict['pPerUnc'])}"
                    f",a/R*={round_to_2(self.p_dict['aRs'], self.p_dict['aRsUnc'])} +/- {round_to_2(self.p_dict['aRsUnc'])}"
                    f",inc={round_to_2(self.p_dict['inc'], self.p_dict['incUnc'])} +/- {round_to_2(self.p_dict['incUnc'])}"
                    f",ecc={round_to_2(self.p_dict['ecc'])}"
                    f",u0={round_to_2(ld0[0], ld0[1])} +/- {round_to_2(ld0[1])}"
                    f",u1={round_to_2(ld1[0], ld1[1])} +/- {round_to_2(ld1[1])}"
                    f",u2={round_to_2(ld2[0], ld2[1])} +/- {round_to_2(ld2[1])}"
                    f",u3={round_to_2(ld3[0], ld3[1])} +/- {round_to_2(ld3[1])}\n"
                    f"#PRIORS-XC={dumps(priors_dict)}\n"  # code yields
                    f"#RESULTS=Tc={round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- {round_to_2(self.fit.errors['tmid'])}"
                    f",Rp/R*={round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- {round_to_2(self.fit.errors['rprs'])}"
                    f",a/R*={round_to_2(self.fit.parameters['ars'], self.fit.errors['ars'])} +/- {round_to_2(self.fit.errors['ars'])}"
                    f",Am1={round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- {round_to_2(self.fit.errors['a1'])}"
                    f",Am2={round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- {round_to_2(self.fit.errors['a2'])}\n"
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


def aavso_dicts(planet_dict, fit, info_dict, durs, ld0, ld1, ld2, ld3):
    priors = {
        'Period': {
            'value': str(round_to_2(planet_dict['pPer'], planet_dict['pPerUnc'])),
            'uncertainty': str(round_to_2(planet_dict['pPerUnc'])) if planet_dict['pPerUnc'] else planet_dict['pPerUnc'],
            'units': "days"
        },
        'a/R*': {
            'value': str(round_to_2(planet_dict['aRs'], planet_dict['aRsUnc'])),
            'uncertainty': str(round_to_2(planet_dict['aRsUnc'])) if planet_dict['aRsUnc'] else planet_dict['aRsUnc'],
        },
        'inc': {
            'value': str(round_to_2(planet_dict['inc'], planet_dict['incUnc'])),
            'uncertainty': str(round_to_2(planet_dict['incUnc'])) if planet_dict['incUnc'] else planet_dict['incUnc'],
            'units': "degrees"
        },
        'ecc': {
            'value': str(round_to_2(planet_dict['ecc'])),
            'uncertainty': None,
        },
        'u0': {
            'value': str(round_to_2(ld0[0], ld0[1])),
            'uncertainty': str(round_to_2(ld0[1]))
        },
        'u1': {
            'value': str(round_to_2(ld1[0], ld1[1])),
            'uncertainty': str(round_to_2(ld1[1]))
        },
        'u2': {
            'value': str(round_to_2(ld2[0], ld2[1])),
            'uncertainty': str(round_to_2(ld2[1]))
        },
        'u3': {
            'value': str(round_to_2(ld3[0], ld3[1])),
            'uncertainty': str(round_to_2(ld3[1]))
        }
    }

    filter_type = {
        'name': info_dict['filter'],
        'fwhm': [{'value': str(info_dict['wl_min']) if info_dict['wl_min'] else info_dict['wl_min'], 'units': "nm"},
                 {'value': str(info_dict['wl_max']) if info_dict['wl_max'] else info_dict['wl_max'], 'units': "nm"}],
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
        'a/R*': {
            'value': str(round_to_2(fit.parameters['ars'], fit.errors['ars'])),
            'uncertainty': str(round_to_2(fit.errors['ars'])),
        },
        'Am1': {
            'value': str(round_to_2(fit.parameters['a1'], fit.errors['a1'])),
            'uncertainty': str(round_to_2(fit.errors['a1']))
        },
        'Am2': {
            'value': str(round_to_2(fit.parameters['a2'], fit.errors['a2'])),
            'uncertainty': str(round_to_2(fit.errors['a2']))
        },
        'Duration': {
            'value': str(round_to_2(mean(durs))),
            'uncertainty': str(round_to_2(std(durs))),
            'units': "days"
        }
    }

    return priors, filter_type, results
