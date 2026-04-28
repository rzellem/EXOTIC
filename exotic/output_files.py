from json import dump, dumps
from numpy import mean, std
from pathlib import Path
import numpy as np

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


def aavso_airmass_results(fit):
    if getattr(fit, 'airmass_fit_skipped', False):
        return (
            ('Am1', '0', '0'),
            ('Am2', '0', '0'),
        )

    if 'a0' in fit.parameters:
        first_result = (
            'A0',
            str(round_to_2(fit.parameters['a0'], fit.errors['a0'])),
            str(round_to_2(fit.errors['a0'])),
        )
    else:
        first_result = (
            'Am1',
            str(round_to_2(fit.parameters['a1'], fit.errors['a1'])),
            str(round_to_2(fit.errors['a1'])),
        )

    return (
        first_result,
        (
            'Am2',
            str(round_to_2(fit.parameters.get('a2', 0), fit.errors.get('a2', 0))),
            str(round_to_2(fit.errors.get('a2', 0))),
        ),
    )


def aavso_detrend_model(fit):
    if getattr(fit, 'airmass_fit_skipped', False):
        return np.ones(len(fit.time), dtype=float)
    return np.asarray(fit.airmass_model, dtype=float)


def finite_float(value, default=np.nan):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if np.isfinite(value) else default


def format_parameter_with_error(value, error):
    value = finite_float(value)
    error = finite_float(error)
    if not np.isfinite(value):
        return None
    if np.isfinite(error) and error >= 0:
        return f"{round_to_2(value, error)} +/- {round_to_2(error)}"
    return f"{round_to_2(value)} +/- n/a"


def fit_impact_parameter_value_error(fit):
    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    sample_parameters = getattr(fit, 'sample_parameters', {}) or {}
    sample_errors = getattr(fit, 'sample_errors', {}) or {}

    if 'b' in sample_parameters:
        impact_parameter = finite_float(sample_parameters.get('b'))
        impact_error = finite_float(sample_errors.get('b'))
        if np.isfinite(impact_parameter):
            return impact_parameter, impact_error

    if 'b' in parameters:
        impact_parameter = finite_float(parameters.get('b'))
        impact_error = finite_float(errors.get('b'))
        if np.isfinite(impact_parameter):
            return impact_parameter, impact_error

    ars = finite_float(parameters.get('ars'))
    inc = finite_float(parameters.get('inc'))
    if not np.isfinite(ars) or not np.isfinite(inc):
        return np.nan, np.nan

    ecc = finite_float(parameters.get('ecc'), 0.0)
    omega = np.deg2rad(finite_float(parameters.get('omega'), 0.0))
    denominator = 1.0 + ecc * np.sin(omega)
    if not np.isfinite(denominator) or np.isclose(denominator, 0.0):
        return np.nan, np.nan

    scale_factor = (1.0 - ecc ** 2) / denominator
    inc_rad = np.deg2rad(inc)
    impact_parameter = scale_factor * ars * np.cos(inc_rad)

    ars_error = finite_float(errors.get('ars'))
    inc_error = finite_float(errors.get('inc'))
    if np.isfinite(ars_error) and np.isfinite(inc_error):
        impact_error = np.hypot(
            scale_factor * np.cos(inc_rad) * ars_error,
            scale_factor * ars * np.sin(inc_rad) * np.deg2rad(inc_error),
        )
    else:
        impact_error = np.nan

    return float(impact_parameter), float(impact_error) if np.isfinite(impact_error) else np.nan


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

    def final_planetary_params(self, phot_opt, vsp_params, comp_star=None, comp_coords=None, min_aper=None,
                               min_annul=None, adaptive_summary=None):
        params_file = self.dir / "temp" / f"FinalParams_{self.p_dict['pName']}_{self.i_dict['date']}.json"

        transit_qc = getattr(self.fit, 'transit_qc', None)
        qc_residual_scatter = np.nan
        if isinstance(transit_qc, dict):
            qc_residual_scatter = transit_qc.get('residual_scatter', np.nan)
        if not np.isfinite(qc_residual_scatter):
            residuals = np.asarray(getattr(self.fit, 'residuals', np.array([])), dtype=float)
            data = np.asarray(getattr(self.fit, 'data', np.array([])), dtype=float)
            if residuals.size and data.size:
                if residuals.shape == data.shape:
                    median_flux = np.nanmedian(data)
                    if np.isfinite(median_flux) and median_flux != 0:
                        qc_residual_scatter = float(np.std(residuals) / median_flux)
                elif residuals.size == 1:
                    median_flux = np.nanmedian(data)
                    if np.isfinite(median_flux) and median_flux != 0:
                        qc_residual_scatter = float(abs(residuals.reshape(-1)[0]) / median_flux)

        params_num = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- "
                                       f"{round_to_2(self.fit.errors['tmid'])} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/R*)": f"{round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- "
                                                         f"{round_to_2(self.fit.errors['rprs'])}",
            "Transit depth (Rp/Rs)^2": f"{round_to_2(100. * (self.fit.parameters['rprs'] ** 2.))} +/- "
                                       f"{round_to_2(100. * 2. * self.fit.parameters['rprs'] * self.fit.errors['rprs'])} [%]",
            "Orbital Inclination (inc)": f"{round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- "
                                                   f"{round_to_2(self.fit.errors['inc'])} ",
        }
        ars_text = format_parameter_with_error(
            self.fit.parameters.get('ars'),
            self.fit.errors.get('ars'),
        )
        if ars_text is not None:
            params_num["Ratio of Distance to Stellar Radius (a/Rs)"] = ars_text
        impact_parameter, impact_error = fit_impact_parameter_value_error(self.fit)
        impact_text = format_parameter_with_error(impact_parameter, impact_error)
        if impact_text is not None:
            params_num["Impact Parameter (b)"] = impact_text
        if np.isfinite(qc_residual_scatter):
            params_num["Residual scatter around full model fit"] = f"{qc_residual_scatter * 100.0:.4f} %"
        if getattr(self.fit, 'airmass_fit_skipped', False):
            params_num["Airmass correction"] = getattr(
                self.fit,
                'airmass_correction_note',
                "Skipped; no airmass correction applied.",
            )
        else:
            if 'a0' in self.fit.parameters:
                params_num["Baseline flux (a0)"] = (
                    f"{round_to_2(self.fit.parameters['a0'], self.fit.errors['a0'])} +/- "
                    f"{round_to_2(self.fit.errors['a0'])}"
                )
            else:
                params_num["Flux normalization (a1)"] = (
                    f"{round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- "
                    f"{round_to_2(self.fit.errors['a1'])}"
                )
            params_num["Airmass coefficient 2 (a2)"] = (
                f"{round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- "
                f"{round_to_2(self.fit.errors['a2'])}"
            )

        if isinstance(transit_qc, dict) and transit_qc:
            qc_status = transit_qc.get('status')
            qc_summary = transit_qc.get('summary')
            qc_notes = transit_qc.get('notes') or []
            qc_delta_bic = transit_qc.get('delta_bic', np.nan)
            qc_delta_chi2 = transit_qc.get('delta_chi2', np.nan)
            qc_rprs_sigma = transit_qc.get('rprs_sigma', np.nan)
            qc_duration_ratio = transit_qc.get('duration_ratio', np.nan)
            qc_eebls_depth_snr = transit_qc.get('eebls_depth_snr', np.nan)
            qc_deviation_metric = transit_qc.get('deviation_from_expected_value', np.nan)
            qc_tmid_deviation_sigma = transit_qc.get('tmid_deviation_sigma', np.nan)
            qc_tmid_deviation_minutes = transit_qc.get('tmid_deviation_minutes', np.nan)
            qc_tmid_threshold_minutes = transit_qc.get('tmid_deviation_threshold_minutes', np.nan)
            qc_expected_tmid_unc_minutes = transit_qc.get('expected_tmid_unc_minutes', np.nan)
            qc_rprs_deviation_sigma = transit_qc.get('rprs_deviation_sigma', np.nan)
            qc_sigma_threshold = transit_qc.get('deviation_sigma_threshold', np.nan)
            qc_ktmf = transit_qc.get('ktmf_metric', np.nan)
            qc_ktmf_contributions = transit_qc.get('ktmf_contributions') or []

            if qc_status:
                params_num["Transit detection QC"] = str(qc_status).upper()
            if qc_summary:
                params_num["Transit vs flat model"] = qc_summary
            if np.isfinite(qc_delta_bic):
                params_num["Transit vs flat Delta BIC"] = f"{qc_delta_bic:.2f}"
            if np.isfinite(qc_delta_chi2):
                params_num["Transit vs flat Delta chi2"] = f"{qc_delta_chi2:.2f}"
            if np.isfinite(qc_rprs_sigma):
                params_num["Transit depth significance"] = f"{qc_rprs_sigma:.2f} sigma"
            if np.isfinite(qc_duration_ratio):
                params_num["Transit duration consistency"] = f"{qc_duration_ratio:.2f}x modeled duration"
            if np.isfinite(qc_eebls_depth_snr):
                params_num["EEBLS depth SNR"] = f"{qc_eebls_depth_snr:.2f}"
            if np.isfinite(qc_deviation_metric):
                params_num["Deviation From Expected Value"] = f"{qc_deviation_metric:.2f} / 1.00"
            if np.isfinite(qc_sigma_threshold):
                params_num["Expected-value QC threshold"] = f"{qc_sigma_threshold:.2f} sigma"
            if np.isfinite(qc_tmid_deviation_sigma):
                params_num["Expected-value Tmid deviation"] = f"{qc_tmid_deviation_sigma:.2f} sigma"
            if np.isfinite(qc_tmid_deviation_minutes):
                params_num["Expected-value Tmid offset"] = f"{qc_tmid_deviation_minutes:.2f} minutes"
            if np.isfinite(qc_expected_tmid_unc_minutes):
                params_num["Expected-value Tmid uncertainty"] = f"{qc_expected_tmid_unc_minutes:.2f} minutes"
            if np.isfinite(qc_tmid_threshold_minutes):
                params_num["Expected-value Tmid QC window"] = f"{qc_tmid_threshold_minutes:.2f} minutes"
            if np.isfinite(qc_rprs_deviation_sigma):
                params_num["Expected-value Rp/R* deviation"] = f"{qc_rprs_deviation_sigma:.2f} sigma"
            if np.isfinite(qc_ktmf):
                params_num["KTMF"] = f"{qc_ktmf:.2f} / 5.00"
            for contribution_index, contribution in enumerate(qc_ktmf_contributions, start=1):
                label = contribution.get('label', f'Component {contribution_index}')
                detail = contribution.get('detail') or 'n/a'
                available = bool(contribution.get('available'))
                points = float(contribution.get('points', 0.0) or 0.0)
                max_points = float(contribution.get('max_points', 0.0) or 0.0)
                score = contribution.get('score', np.nan)
                if available and np.isfinite(score):
                    params_num[f"KTMF contribution {contribution_index}"] = (
                        f"{label}: +{points:.2f}/{max_points:.2f} (score={score:.2f}; {detail})"
                    )
                else:
                    params_num[f"KTMF contribution {contribution_index}"] = (
                        f"{label}: +0.00/0.00 (unavailable; {detail})"
                    )
            if qc_notes:
                params_num["Transit QC notes"] = " ".join(str(note) for note in qc_notes)

        if vsp_params:
            params_num["Variable Reference Star"] = f"AAVSO Label: {vsp_params[0]['cname']}, " + \
                                                    f"Position: {vsp_params[0]['pos']}"

        if phot_opt:
            phot_ext = {"Best Comparison Star": f"#{comp_star} - {comp_coords}" if min_aper >= 0 else str(comp_star)}
            if min_aper == 0:
                phot_ext["Optimal Method"] = "PSF photometry"
            else:
                if adaptive_summary:
                    phot_ext["Adaptive Aperture Scale"] = f"{adaptive_summary['aperture_sigma']:.2f} sigma"
                    phot_ext["Adaptive Annulus Scale"] = f"{adaptive_summary['annulus_sigma']:.2f} sigma"
                    phot_ext["Optimal Aperture"] = (
                        f"{adaptive_summary['aperture_median']:.2f} +/- {adaptive_summary['aperture_std']:.2f} px"
                    )
                    phot_ext["Aperture Range"] = (
                        f"{adaptive_summary['aperture_min']:.2f} to {adaptive_summary['aperture_max']:.2f} px"
                    )
                    phot_ext["Optimal Annulus"] = (
                        f"{adaptive_summary['annulus_median']:.2f} +/- {adaptive_summary['annulus_std']:.2f} px"
                    )
                    phot_ext["Annulus Range"] = (
                        f"{adaptive_summary['annulus_min']:.2f} to {adaptive_summary['annulus_max']:.2f} px"
                    )
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
        aavso_airmass_terms = aavso_airmass_results(self.fit)
        detrend_model = aavso_detrend_model(self.fit)
        obs_name = format_aavso_header_value(self.i_dict.get('obs_name'))
        obs_name_header = f"#OBSNAME={obs_name}\n" if obs_name else ""
        gaia_dist = format_aavso_header_value(self.p_dict.get('dist'))
        gaia_pmra = format_aavso_header_value(self.p_dict.get('pm_ra'))
        gaia_pmdec = format_aavso_header_value(self.p_dict.get('pm_dec'))
        gaia_dist_header = f"#GAIADIST={gaia_dist}\n" if gaia_dist else ""
        gaia_pmra_header = f"#GAIAPMRA={gaia_pmra}\n" if gaia_pmra else ""
        gaia_pmdec_header = f"#GAIAPMDEC={gaia_pmdec}\n" if gaia_pmdec else ""

        params_file = self.dir / f"AAVSO_{self.p_dict['pName']}_{self.i_dict['date']}.txt"

        with params_file.open('w', encoding="utf-8") as f:
            f.write("#TYPE=EXOPLANET\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SECONDARY_OBSCODES={self.i_dict['second_obs']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSDATE={format_aavso_header_value(self.i_dict.get('date'))}\n"
                    f"{obs_name_header}"
                    f"#OBSTYPE={self.i_dict['camera']}\n"
                    f"#STAR_NAME={self.p_dict['sName']}\n"  # code yields
                    f"#EXOPLANET_NAME={self.p_dict['pName']}\n"  # code yields
                    f"#BINNING={self.i_dict['pixel_bin']}\n"  # user input
                    f"#EXPOSURE_TIME={self.i_dict.get('exposure', -1)}\n"  # UI
                    f"#OBSLAT={format_aavso_header_value(self.i_dict.get('lat'))}\n"
                    f"#OBSLON={format_aavso_header_value(self.i_dict.get('long'))}\n"
                    f"#OBSELEV={format_aavso_header_value(self.i_dict.get('elev'))}\n"
                    f"{gaia_dist_header}"
                    f"{gaia_pmra_header}"
                    f"{gaia_pmdec_header}"
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
                    f",{aavso_airmass_terms[0][0]}={aavso_airmass_terms[0][1]} +/- {aavso_airmass_terms[0][2]}"
                    f",{aavso_airmass_terms[1][0]}={aavso_airmass_terms[1][1]} +/- {aavso_airmass_terms[1][2]}\n"
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
                        f"{round(detrend_model[aavsoC], 7)}\n")
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
                    f"#OBSDATE={format_aavso_header_value(self.i_dict.get('date'))}\n"
                    f"#OBSTYPE={self.i_dict['camera']}\n"
                    f"#OBSLAT={format_aavso_header_value(self.i_dict.get('lat'))}\n"
                    f"#OBSLON={format_aavso_header_value(self.i_dict.get('long'))}\n"
                    f"#OBSELEV={format_aavso_header_value(self.i_dict.get('elev'))}\n")
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
    aavso_airmass_terms = aavso_airmass_results(fit)
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
        'filter_width': {
            'left_side_wavelength': {
                'value': str(info_dict['wl_min']) if info_dict['wl_min'] else info_dict['wl_min'],
                'units': "nm"
            },
            'right_side_wavelength': {
                'value': str(info_dict['wl_max']) if info_dict['wl_max'] else info_dict['wl_max'],
                'units': "nm"
            }
        },
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
        'Am2': {
            'value': aavso_airmass_terms[1][1],
            'uncertainty': aavso_airmass_terms[1][2]
        },
        'Duration': {
            'value': str(round_to_2(mean(durs))),
            'uncertainty': str(round_to_2(std(durs))),
            'units': "days"
        }
    }

    results[aavso_airmass_terms[0][0]] = {
        'value': aavso_airmass_terms[0][1],
        'uncertainty': aavso_airmass_terms[0][2]
    }

    return priors, filter_type, results


def format_aavso_header_value(value):
    if value is None:
        return ""
    if isinstance(value, str):
        stripped = value.strip()
        return "" if stripped.lower() in ('', 'n/a', 'na', 'null', 'none') else stripped
    return str(value)


def save_comp_star_calibration_summary(save_dir, target_name, date, method_label, field_score,
                                       comp_summaries, best_comp_index):
    temp_dir = Path(save_dir) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    summary_file = temp_dir / f"CompStarCalibrationSummary_{target_name}_{date}.csv"

    with summary_file.open('w') as handle:
        handle.write(f"# Comparison-star calibration summary for {target_name}\n")
        handle.write(f"# Method,{method_label}\n")
        if field_score is not None and field_score == field_score:
            handle.write(f"# Field suitability score,{field_score}\n")
        else:
            handle.write("# Field suitability score,\n")
        handle.write(f"# Selected comparison star,{'' if best_comp_index is None else best_comp_index + 1}\n")
        handle.write("comp_star,x_pixel,y_pixel,selected,suitability_score,ensemble_score,pairwise_median_score,"
                     "pairwise_max_score,self_score,valid_pair_count,coverage_count,coverage_peer_median,"
                     "coverage_min_required,coverage_rejected,suitability_outlier_rejected\n")

        for summary in comp_summaries:
            position = summary.get('position') or [None, None]
            values = [
                summary.get('label', ''),
                position[0],
                position[1],
                str(bool(summary.get('selected'))).lower(),
                summary.get('aggregate_score'),
                summary.get('ensemble_score'),
                summary.get('pairwise_median_score'),
                summary.get('pairwise_max_score'),
                summary.get('self_score'),
                summary.get('valid_pair_count'),
                summary.get('coverage_count'),
                summary.get('coverage_reference_count'),
                summary.get('coverage_min_required_count'),
                summary.get('coverage_rejected'),
                summary.get('suitability_outlier_rejected'),
            ]
            handle.write(",".join("" if value is None else str(value) for value in values) + "\n")

    return summary_file
