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
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus


class OutputFiles:
    def __init__(self, fit, p_dict, i_dict, durs):
        self.fit = fit
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.durs = durs
        self.dir = Path(self.i_dict['save'])

    def final_lightcurve(self, phase):
        params_file = self.dir / "temp" / f"FinalLightCurve_{self.p_dict['pName']}_{self.i_dict['date']}.csv"

        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {self.p_dict['pName']}\n")
            f.write("# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n")

            for bjd, phase, flux, fluxerr, model, am in zip(self.fit.time, phase, self.fit.detrended,
                                                            self.fit.dataerr / self.fit.airmass_model,
                                                            self.fit.transit, self.fit.airmass_model):
                f.write(f"{bjd}, {phase}, {flux}, {fluxerr}, {model}, {am}\n")

    def final_planetary_params(self, phot_opt, vsp_params, comp_star=None, comp_coords=None, min_aper=None, min_annul=None):
        params_file = self.dir / "temp" / f"FinalParams_{self.p_dict['pName']}_{self.i_dict['date']}.json"

        params_num = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- "
                                       f"{round_to_2(self.fit.errors['tmid'])} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/R*)": f"{round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- "
                                                         f"{round_to_2(self.fit.errors['rprs'])}",
            "Transit depth (Rp/Rs)^2": f"{round_to_2(100. * (self.fit.parameters['rprs'] ** 2.))} +/- "
                                       f"{round_to_2(100. * 2. * self.fit.parameters['rprs'] * self.fit.errors['rprs'])} [%]",
            "Orbital Inclination (inc)": f"{round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- "
                                                  f"{round_to_2(self.fit.errors['inc'])} ",
            "Airmass coefficient 1 (a1)": f"{round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- "
                                          f"{round_to_2(self.fit.errors['a1'])}",
            "Airmass coefficient 2 (a2)": f"{round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- "
                                          f"{round_to_2(self.fit.errors['a2'])}",
            "Scatter in the residuals of the lightcurve fit is": f"{round_to_2(100. * std(self.fit.residuals / median(self.fit.data)))} %",
        }

        if vsp_params:
            params_num["Variable Reference Star"] = f"AAVSO Label: {vsp_params[0]['cname']}, " + \
                                                    f"Position: {vsp_params[0]['pos']}"

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

    def aavso(self, comp_star, airmasses, ld0, ld1, ld2, ld3, epw_md5):
        priors_dict, filter_dict, results_dict = aavso_dicts(self.p_dict, self.fit, self.i_dict, self.durs,
                                                             ld0, ld1, ld2, ld3)

        params_file = self.dir / f"AAVSO_{self.p_dict['pName']}_{self.i_dict['date']}.txt"

        with params_file.open('w', encoding="utf-8") as f:
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
                    f",Rp/R*={round_to_2(self.p_dict['rprs'], self.p_dict['rprsUnc'])} +/- {round_to_2(self.p_dict['rprsUnc'])}"
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
                    f",inc={round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- {round_to_2(self.fit.errors['inc'])}"
                    f",Am1={round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- {round_to_2(self.fit.errors['a1'])}"
                    f",Am2={round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- {round_to_2(self.fit.errors['a2'])}\n"
                    f"#RESULTS-XC={dumps(results_dict)}\n")  # code yields

            if epw_md5:
                f.write(f"#EPW_MD5-XC={dumps({'epw_checkout_md5': epw_md5})}\n")

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
    def plate_status(self, plate_status: PlateStatus):
        plate_status_file = self.dir / "temp" / f"PlateStatus_{self.p_dict['pName']}_{self.i_dict['date']}.csv"
        plate_status.writePlateStatus(plate_status_file)

class AIDOutputFiles:
    def __init__(self, fit, p_dict, i_dict, auid, chart_id, vsp_params):
        self.fit = fit
        self.auid = auid
        self.chart_id = chart_id
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.dir = Path(self.i_dict['save'])
        self.vsp_params = vsp_params

    def aavso(self):
        params_file = self.dir / f"AID_AAVSO_{self.p_dict['sName']}_{self.i_dict['date']}.txt"
        with params_file.open('w', encoding="utf-8") as f:
            f.write("#TYPE=EXTENDED\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE=JD\n"  # fixed
                    f"#OBSTYPE={self.i_dict['camera']}\n")
            f.write(
                "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                "This work is supported by NASA under award number NNX16AC65A to the "
                "Space Telescope Science Institute.\n"
                "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                "aavso.org/data-usage-guidelines\n")

            f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
            for vsp_p in self.vsp_params:
                f.write(f"{self.auid},{round(vsp_p['time'], 5)},{round(vsp_p['mag'], 5)},{round(vsp_p['mag_err'], 5)},"
                        f"{self.i_dict['filter']},NO,STD,{vsp_p['cname']},{round(vsp_p['cmag'], 5)},na,na," 
                        f"{round(vsp_p['airmass'], 7)},na,{self.chart_id},na\n")


def aavso_dicts(planet_dict, fit, info_dict, durs, ld0, ld1, ld2, ld3):
    priors = {
        'Period': {
            'value': str(round_to_2(planet_dict['pPer'], planet_dict['pPerUnc'])),
            'uncertainty': str(round_to_2(planet_dict['pPerUnc'])) if planet_dict['pPerUnc'] else planet_dict['pPerUnc'],
            'units': "days"
        },
        'Rp/R*': {
            'value': str(round_to_2(planet_dict['rprs'], planet_dict['rprsUnc'])),
            'uncertainty': str(round_to_2(planet_dict['rprsUnc'])) if planet_dict['rprsUnc'] else planet_dict['rprsUnc'],
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
        'desc': info_dict['filter_desc'],
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
        'inc': {
            'value': str(round_to_2(fit.parameters['inc'], fit.errors['inc'])),
            'uncertainty': str(round_to_2(fit.errors['inc'])),
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
