from astropy.visualization import astropy_mpl_style, ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch, SquaredStretch, SqrtStretch, LogStretch
import inspect
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path

try:
    from utils import (
        filename_date_token,
        is_usable_apparent_magnitude,
        magnitude_text,
        normalized_magnitude_error,
        safe_output_filename,
    )
except ImportError:
    from .utils import (
        filename_date_token,
        is_usable_apparent_magnitude,
        magnitude_text,
        normalized_magnitude_error,
        safe_output_filename,
    )

try:
    from output_files import (
        empirical_red_noise_error_scale,
        fit_empirical_transit_uncertainty,
        fit_impact_parameter_value_error,
        fit_parameter_model_data_uncertainty,
    )
except ImportError:
    from .output_files import (
        empirical_red_noise_error_scale,
        fit_empirical_transit_uncertainty,
        fit_impact_parameter_value_error,
        fit_parameter_model_data_uncertainty,
    )

plt.style.use(astropy_mpl_style)


def _dated_plot_filename(prefix, *parts, date, extension):
    return safe_output_filename(prefix, *parts, filename_date_token(date), extension=extension)


def _working_artifacts_dir(save):
    output_dir = Path(save) / "working_artifacts"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# Plots of the centroid positions as a function of time
def plot_centroids(x_targ, y_targ, x_ref, y_ref, times, target_name, save, date):
    fig, axs = plt.subplots(3, 2, figsize=(12, 10))

    axs[0, 0].set_title(f"{target_name} X-Centroid Position", fontsize=14)
    axs[0, 0].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[0, 0].set_ylabel("X-Centroid [px]", fontsize=12)
    axs[0, 0].plot(times - np.nanmin(times), x_targ, 'k.')

    axs[0, 1].set_title(f"{target_name} Y-Centroid Position", fontsize=14)
    axs[0, 1].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[0, 1].set_ylabel("Y-Centroid [px]", fontsize=12)
    axs[0, 1].plot(times - np.nanmin(times), y_targ, 'k.')

    axs[1, 0].set_title(f"Comparison Star X-Centroid Position", fontsize=14)
    axs[1, 0].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[1, 0].set_ylabel("X-Centroid [px]", fontsize=12)
    axs[1, 0].plot(times - np.nanmin(times), x_ref, 'k.')

    axs[1, 1].set_title(f"Comparison Star Y-Centroid Position", fontsize=14)
    axs[1, 1].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[1, 1].set_ylabel("Y-Centroid [px]", fontsize=12)
    axs[1, 1].plot(times - np.nanmin(times), y_ref, 'k.')

    axs[2, 0].set_title("Distance between Target and Comparison X-Centroids", fontsize=14)
    axs[2, 0].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[2, 0].set_ylabel("X-Centroid Distance [px]", fontsize=12)
    for e in range(len(x_targ)):
        axs[2, 0].plot(times[e] - np.nanmin(times), abs(x_targ[e] - x_ref[e]), 'k.')

    axs[2, 1].set_title("Distance between Target and Comparison Y-Centroids", fontsize=14)
    axs[2, 1].set_xlabel(f"Time [BJD_TDB-{np.nanmin(times):.5f}]", fontsize=12)
    axs[2, 1].set_ylabel("Y-Centroid Distance [px]", fontsize=12)
    for e in range(len(y_targ)):
        axs[2, 1].plot(times[e] - np.nanmin(times), abs(y_targ[e] - y_ref[e]), 'k.')

    plt.tight_layout()
    plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename(
        "CentroidPositions&Distances",
        target_name,
        date=date,
        extension="pdf",
    ))
    plt.close()

def plot_fov(aper, annulus, sigma, x_targ, y_targ, x_ref, y_ref, image, image_scale, targ_name, save, date,
             opt_method, min_aper_fov, min_annulus_fov, sky_inner_radius=None, sky_outer_radius=None):

    ref_circle, ref_circle_sky = None, None
    picframe = 10. * (aper + 15. * sigma)

    pltx = [max([0, min([x_targ, x_ref]) - picframe]), min([np.shape(image)[1], max([x_targ, x_ref]) + picframe])]
    plty = [max([0, min([y_targ, y_ref]) - picframe]), min([np.shape(image)[0], max([y_targ, y_ref]) + picframe])]

    for stretch in [LinearStretch(), SquaredStretch(), SqrtStretch(), LogStretch()]:
        fig, ax = plt.subplots()

        # Set color for target and reference outer circles based on opt_method
        if opt_method == "Aperture":
            outer_circle_color = 'r'
        else:
            outer_circle_color = 'lime'

        # Create the target circles
        # We are using abs(aper) to account for a negative aperture in case EXOTIC is not using a comparison star
        if sky_inner_radius is None or sky_outer_radius is None:
            local_sky_inner_radius = abs(aper) + 2.0
            if np.isfinite(sigma) and sigma > 0:
                local_sky_inner_radius = max(local_sky_inner_radius, 3.0 * 2.355 * float(sigma))
            local_sky_outer_radius = max(
                local_sky_inner_radius + annulus,
                np.sqrt(local_sky_inner_radius ** 2 + 250.0 / np.pi),
            )
        else:
            local_sky_inner_radius = float(sky_inner_radius)
            local_sky_outer_radius = float(sky_outer_radius)

        target_circle = plt.Circle((x_targ, y_targ), abs(aper), color=outer_circle_color, fill=False, ls='-')
        target_circle_sky_inner = plt.Circle((x_targ, y_targ), local_sky_inner_radius, color=outer_circle_color, fill=False, ls='--')
        target_circle_sky_outer = plt.Circle((x_targ, y_targ), local_sky_outer_radius, color=outer_circle_color, fill=False, ls='-')

        # IF EXOTIC is using a comparison star, create its circles
        if aper >= 0:
            ref_circle = plt.Circle((x_ref, y_ref), aper, color=outer_circle_color, fill=False, ls='-')
            ref_circle_sky_inner = plt.Circle((x_ref, y_ref), local_sky_inner_radius, color=outer_circle_color, fill=False, ls='--')
            ref_circle_sky = plt.Circle((x_ref, y_ref), local_sky_outer_radius, color=outer_circle_color, fill=False, ls='-')

        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(image)

        norm = ImageNormalize(image, interval=interval, stretch=stretch, vmin=vmin, vmax=vmax)

        im = plt.imshow(image, norm=norm, origin='lower', cmap='Greys_r', interpolation=None)
        fig.colorbar(im)

        ax.add_artist(target_circle)
        ax.add_artist(target_circle_sky_inner)
        ax.add_artist(target_circle_sky_outer)
        ax.text(x_targ + local_sky_outer_radius + 5, y_targ, targ_name, color='w', fontsize=10,
                path_effects=[path_effects.withStroke(linewidth=2, foreground='black')])

        if aper >= 0: #EXOTIC is using a comparison star
            ax.add_artist(ref_circle)
            ax.add_artist(ref_circle_sky_inner)
            ax.add_artist(ref_circle_sky)
            ax.text(x_ref + local_sky_outer_radius + 5, y_ref, 'Comp Star', color='w', fontsize=10,
                    path_effects=[path_effects.withStroke(linewidth=2, foreground='black')])

        handles = []
        if opt_method == "PSF":
            label_aper = "PSF Photometry"
        else:
            label_aper = (
                f"{opt_method} Photometry\n"
                f"(Min Aper: {abs(min_aper_fov):.2f} px)\n"
                f"(Min Annulus: {min_annulus_fov:.2f} px)"
            )
        
        if opt_method == "Aperture":
            aperture_line = Line2D([], [], color=outer_circle_color, linestyle='-', label=label_aper)
            handles.append(aperture_line)
        elif opt_method == "PSF":
            psf_line = Line2D([], [], color=outer_circle_color, linestyle='-', label=label_aper)
            handles.append(psf_line)

        plt.title(f"FOV for {targ_name}\n({image_scale})")
        plt.xlabel("x-axis [pixel]")
        plt.ylabel("y-axis [pixel]")
        plt.xlim(pltx[0], pltx[1])
        plt.ylim(plty[0], plty[1])
        ax.grid(False)

        if handles:
            l = plt.legend(handles=handles, framealpha=0.75)
            for text in l.get_texts():
                text.set_color("k")
                text.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground='white')])

        apos = '\''
        Path(save).mkdir(parents=True, exist_ok=True)
        _working_artifacts_dir(save)

        stretch_name = str(stretch.__class__).split('.')[-1].split(apos)[0]
        plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename(
            "FOV",
            targ_name,
            stretch_name,
            date=date,
            extension="pdf",
        ), bbox_inches='tight')
        plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename(
            "FOV",
            targ_name,
            stretch_name,
            date=date,
            extension="png",
        ), bbox_inches='tight')
        plt.close()


def plot_flux(times, targ, targ_unc, ref, ref_unc, norm_flux, norm_unc, airmass, targ_name, save, date):
    plt.figure()
    plt.title(f"{targ_name} Raw Flux Values {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Flux [ADU]")
    plt.errorbar(times, targ, yerr=targ_unc, linestyle='None', fmt='-o')
    plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename("TargetRawFlux", targ_name, date=date, extension="pdf"))
    plt.close()

    plt.figure()
    plt.title(f"Comparison Star Raw Flux Values {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Flux [ADU]")
    plt.errorbar(times, ref, yerr=ref_unc, linestyle='None', fmt='-o')
    plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename("CompRawFlux", targ_name, date=date, extension="pdf"))
    plt.close()

    # Plots final reduced light curve (after the 3 sigma clip)
    plt.figure()
    plt.title(f"{targ_name} Normalized Flux vs. Time {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Normalized Flux")
    plt.errorbar(times, norm_flux, yerr=norm_unc, linestyle='None', fmt='-bo')
    plt.savefig(_working_artifacts_dir(save) / _dated_plot_filename("NormalizedFluxTime", targ_name, date=date, extension="pdf"))
    plt.close()

    # Save normalized flux to text file prior to NS
    params_file = _working_artifacts_dir(save) / _dated_plot_filename("NormalizedFlux", targ_name, date=date, extension="txt")
    with params_file.open('w') as f:
        f.write("BJD,Norm Flux,Norm Err,AM\n")

        for ti, fi, erri, ami in zip(times, norm_flux, norm_unc, airmass):
            f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")


def plot_comp_star_pairwise_matrix(pairwise_matrix, best_comp_index, targ_name, save, date, method_label):
    matrix = np.asarray(pairwise_matrix, dtype=float)
    if matrix.size == 0:
        return

    temp_dir = _working_artifacts_dir(save)

    fig, ax = plt.subplots(figsize=(max(6, matrix.shape[0] * 1.3), max(5, matrix.shape[0] * 1.1)))
    plot_matrix = np.ma.masked_invalid(matrix * 100.0)
    im = ax.imshow(plot_matrix, origin='upper', cmap='viridis')
    fig.colorbar(im, ax=ax, label="Residual Scatter [%]")

    labels = [f"Comp {index + 1}" for index in range(matrix.shape[0])]
    ax.set_xticks(np.arange(matrix.shape[0]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticklabels(labels)
    ax.set_title(f"{targ_name} Comparison-Star Pairwise Scatter\n{method_label}")

    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            if np.isfinite(value):
                ax.text(col, row, f"{value * 100.0:.3f}", ha='center', va='center', color='white', fontsize=8)

    if best_comp_index is not None and 0 <= best_comp_index < matrix.shape[0]:
        ax.add_patch(plt.Rectangle((best_comp_index - 0.5, best_comp_index - 0.5), 1, 1,
                                   fill=False, edgecolor='tomato', linewidth=2.5))

    ax.set_xlabel("Reference Comparison Star")
    ax.set_ylabel("Candidate Comparison Star")
    fig.tight_layout()
    fig.savefig(temp_dir / _dated_plot_filename("CompStarPairwiseScatter", targ_name, date=date, extension="png"), bbox_inches="tight")
    fig.savefig(temp_dir / _dated_plot_filename("CompStarPairwiseScatter", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_comp_star_calibration_series(times, comp_summaries, targ_name, save, date, method_label):
    if not comp_summaries:
        return

    times = np.asarray(times, dtype=float)
    temp_dir = _working_artifacts_dir(save)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, 10))

    fig_height = max(3.2, 2.4 * len(comp_summaries))
    fig, axes = plt.subplots(len(comp_summaries), 1, figsize=(12, fig_height), sharex=True)
    if len(comp_summaries) == 1:
        axes = [axes]

    for axis, summary in zip(axes, comp_summaries):
        _draw_comp_star_calibration_axis(axis, times, summary, colors)

    axes[-1].set_xlabel("Time [BJD_TDB]")
    fig.suptitle(f"{targ_name} Comparison-Star Calibration Curves\n{method_label}", y=1.01)
    fig.tight_layout()
    fig.savefig(temp_dir / _dated_plot_filename("CompStarCalibrationCurves", targ_name, date=date, extension="png"), bbox_inches="tight")
    fig.savefig(temp_dir / _dated_plot_filename("CompStarCalibrationCurves", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_individual_comp_star_calibration_series(times, comp_summaries, targ_name, save, date, method_label):
    if not comp_summaries:
        return

    times = np.asarray(times, dtype=float)
    temp_dir = _working_artifacts_dir(save)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, 10))

    for summary in comp_summaries:
        fig, axis = plt.subplots(figsize=(12, 4))
        _draw_comp_star_calibration_axis(axis, times, summary, colors)
        axis.set_xlabel("Time [BJD_TDB]")
        fig.suptitle(f"{targ_name} {summary['label']} Calibration Curves\n{method_label}")
        fig.tight_layout()
        label_slug = summary['label'].replace(" ", "")
        fig.savefig(temp_dir / _dated_plot_filename(
            "CompStarCalibrationCurve",
            label_slug,
            targ_name,
            date=date,
            extension="png",
        ), bbox_inches="tight")
        fig.savefig(temp_dir / _dated_plot_filename(
            "CompStarCalibrationCurve",
            label_slug,
            targ_name,
            date=date,
            extension="pdf",
        ), bbox_inches="tight")
        plt.close(fig)


def plot_comp_star_candidate_lightcurve_fits(candidate_fit_summaries, targ_name, save, date, method_label):
    if not candidate_fit_summaries:
        return

    temp_dir = _working_artifacts_dir(save)

    for summary in candidate_fit_summaries:
        fit = summary.get('fit')
        if fit is None:
            continue

        fig, (ax_lc, ax_res) = _plot_bestfit_for_lightcurve_png(
            fit,
            phase=False,
            show_flux_baseline_label=False,
        )
        selected_text = " selected" if summary.get('selected') else ""
        res_std = summary.get('res_std', np.nan)
        res_std_text = "n/a" if not np.isfinite(res_std) else f"{res_std * 100.0:.3f}%"
        ax_lc.set_title(f"{targ_name} vs {summary['label']}{selected_text}\n{method_label} | scatter={res_std_text}")
        ax_res.set_title("")

        label_slug = summary['label'].replace(" ", "")
        fig.savefig(temp_dir / _dated_plot_filename(
            "CompStarLightCurveFit",
            label_slug,
            targ_name,
            date=date,
            extension="png",
        ), bbox_inches="tight")
        fig.savefig(temp_dir / _dated_plot_filename(
            "CompStarLightCurveFit",
            label_slug,
            targ_name,
            date=date,
            extension="pdf",
        ), bbox_inches="tight")
        plt.close(fig)


def _callable_accepts_keyword(callable_object, keyword):
    try:
        signature = inspect.signature(callable_object)
    except (TypeError, ValueError):
        return False
    if keyword in signature.parameters:
        return True
    return any(
        parameter.kind == inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )


def _plot_bestfit_for_lightcurve_png(fit, **requested_kwargs):
    plotter = fit.plot_bestfit
    plot_kwargs = {
        key: value
        for key, value in requested_kwargs.items()
        if _callable_accepts_keyword(plotter, key)
    }
    return plotter(**plot_kwargs)


def _draw_comp_star_calibration_axis(axis, times, summary, colors):
    axis.axhline(1.0, color='lightgray', lw=1.0, zorder=1)
    ensemble_keep_mask = np.asarray(summary.get('ensemble_frame_keep_mask'), dtype=bool)
    has_ensemble_keep_mask = ensemble_keep_mask.shape == times.shape
    pairwise_series = summary.get('pairwise_ratio_series', {})
    for color_index, (other_label, ratio_series) in enumerate(pairwise_series.items()):
        ratio_series = np.asarray(ratio_series, dtype=float)
        line_ratio = ratio_series.copy()
        if has_ensemble_keep_mask and line_ratio.shape == times.shape:
            line_ratio[~ensemble_keep_mask] = np.nan
        valid_time = np.isfinite(times)
        valid_line = valid_time & np.isfinite(line_ratio)
        if np.any(valid_line):
            axis.plot(times[valid_time], line_ratio[valid_time], color=colors[color_index % len(colors)],
                      alpha=0.55, lw=1.0, label=other_label)

    ensemble_ratio = np.asarray(summary.get('ensemble_ratio_series'), dtype=float)
    ensemble_time_valid = np.isfinite(times)
    ensemble_valid = ensemble_time_valid & np.isfinite(ensemble_ratio)
    if np.any(ensemble_valid):
        line_ratio = ensemble_ratio.copy()
        if has_ensemble_keep_mask and line_ratio.shape == times.shape:
            line_ratio[~ensemble_keep_mask] = np.nan
        line_valid = ensemble_time_valid & np.isfinite(line_ratio)
        if np.any(line_valid):
            axis.plot(times[ensemble_time_valid], line_ratio[ensemble_time_valid], color='black', lw=1.8,
                      label='Ensemble')
        if has_ensemble_keep_mask:
            rejected = ensemble_valid & ~ensemble_keep_mask
            if np.any(rejected):
                axis.scatter(times[rejected], ensemble_ratio[rejected], marker='x', s=42,
                             color='red', linewidths=1.4, label='Ensemble clip')

    selected_text = " selected" if summary.get('selected') else ""
    aggregate = summary.get('aggregate_score', np.nan)
    aggregate_text = "n/a" if not np.isfinite(aggregate) else f"{aggregate * 100.0:.3f}%"
    axis.set_ylabel("Norm Ratio")
    axis.set_title(f"{summary['label']}{selected_text} | suitability={aggregate_text}", loc='left', fontsize=10)
    axis.grid(alpha=0.2)
    handles, labels = axis.get_legend_handles_labels()
    if handles and labels:
        axis.legend(ncol=4, fontsize=8, loc='upper right')


def plot_comp_star_suitability(comp_summaries, targ_name, save, date, method_label):
    if not comp_summaries:
        return

    temp_dir = _working_artifacts_dir(save)

    labels = [summary['label'] for summary in comp_summaries]
    positions = np.arange(len(labels))
    aggregate = np.array([summary.get('aggregate_score', np.nan) for summary in comp_summaries], dtype=float) * 100.0
    ensemble = np.array([summary.get('ensemble_score', np.nan) for summary in comp_summaries], dtype=float) * 100.0
    pairwise = np.array([summary.get('pairwise_median_score', np.nan) for summary in comp_summaries], dtype=float) * 100.0

    fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(labels)), 5))
    width = 0.25
    ax.bar(positions - width, aggregate, width=width, label='Suitability')
    ax.bar(positions, ensemble, width=width, label='Vs ensemble')
    ax.bar(positions + width, pairwise, width=width, label='Pairwise median')

    for position, summary in zip(positions, comp_summaries):
        if summary.get('selected'):
            ax.text(position - width, aggregate[position] if np.isfinite(aggregate[position]) else 0.0, 'selected',
                    rotation=90, va='bottom', ha='center', fontsize=8, color='tomato')

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Residual Scatter [%]")
    ax.set_title(f"{targ_name} Comparison-Star Suitability Summary\n{method_label}")
    ax.legend()
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(temp_dir / _dated_plot_filename("CompStarSuitability", targ_name, date=date, extension="png"), bbox_inches="tight")
    fig.savefig(temp_dir / _dated_plot_filename("CompStarSuitability", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_adaptive_aperture_diagnostics(times, aperture_series, annulus_series, fwhm_series, airmass,
                                       targ_name, save, date, aperture_sigma, annulus_sigma):
    times = np.asarray(times, dtype=float)
    aperture_series = np.asarray(aperture_series, dtype=float)
    annulus_series = np.asarray(annulus_series, dtype=float)
    fwhm_series = np.asarray(fwhm_series, dtype=float)
    airmass = np.asarray(airmass, dtype=float)

    plot_len = min(times.size, aperture_series.size, annulus_series.size, fwhm_series.size, airmass.size)
    if plot_len == 0:
        return

    times = times[:plot_len]
    aperture_series = aperture_series[:plot_len]
    annulus_series = annulus_series[:plot_len]
    fwhm_series = fwhm_series[:plot_len]
    airmass = airmass[:plot_len]

    valid_time = np.isfinite(times)
    valid_aperture = np.isfinite(aperture_series)
    valid_annulus = np.isfinite(annulus_series)
    valid_fwhm = np.isfinite(fwhm_series)
    valid_airmass = np.isfinite(airmass)

    temp_dir = _working_artifacts_dir(save)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    fig.suptitle(
        f"{targ_name} Adaptive Aperture Diagnostics\n"
        f"aper={aperture_sigma:.2f} sigma, annulus={annulus_sigma:.2f} sigma"
    )

    time_mask = valid_time & valid_aperture
    time_zero = np.nanmin(times[time_mask]) if np.any(time_mask) else 0.0
    axes[0, 0].set_title("Aperture Radius vs Time")
    axes[0, 0].set_xlabel(f"Time [BJD_TDB-{time_zero:.5f}]")
    axes[0, 0].set_ylabel("Aperture Radius [px]")
    if np.any(time_mask):
        axes[0, 0].plot(times[time_mask] - time_zero, aperture_series[time_mask], color='tab:blue',
                        marker='o', ms=3, lw=1.1)
    axes[0, 0].grid(alpha=0.25)

    annulus_mask = valid_time & valid_annulus
    axes[0, 1].set_title("Annulus Width vs Time")
    axes[0, 1].set_xlabel(f"Time [BJD_TDB-{time_zero:.5f}]")
    axes[0, 1].set_ylabel("Annulus Width [px]")
    if np.any(annulus_mask):
        axes[0, 1].plot(times[annulus_mask] - time_zero, annulus_series[annulus_mask], color='tab:orange',
                        marker='o', ms=3, lw=1.1)
    axes[0, 1].grid(alpha=0.25)

    fwhm_mask = valid_aperture & valid_fwhm
    axes[1, 0].set_title("Aperture Radius vs Target FWHM")
    axes[1, 0].set_xlabel("Target PSF FWHM [px]")
    axes[1, 0].set_ylabel("Aperture Radius [px]")
    if np.any(fwhm_mask):
        axes[1, 0].scatter(fwhm_series[fwhm_mask], aperture_series[fwhm_mask], color='tab:green', s=18, alpha=0.8)
        order = np.argsort(fwhm_series[fwhm_mask])
        axes[1, 0].plot(fwhm_series[fwhm_mask][order], aperture_series[fwhm_mask][order], color='tab:green',
                        alpha=0.35, lw=1.0)
    axes[1, 0].grid(alpha=0.25)

    airmass_mask = valid_aperture & valid_airmass
    axes[1, 1].set_title("Aperture Radius vs Airmass")
    axes[1, 1].set_xlabel("Airmass")
    axes[1, 1].set_ylabel("Aperture Radius [px]")
    if np.any(airmass_mask):
        axes[1, 1].scatter(airmass[airmass_mask], aperture_series[airmass_mask], color='tab:red', s=18, alpha=0.8)
        order = np.argsort(airmass[airmass_mask])
        axes[1, 1].plot(airmass[airmass_mask][order], aperture_series[airmass_mask][order], color='tab:red',
                        alpha=0.35, lw=1.0)
    axes[1, 1].grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(temp_dir / _dated_plot_filename("AdaptiveApertureDiagnostics", targ_name, date=date, extension="png"), bbox_inches="tight")
    fig.savefig(temp_dir / _dated_plot_filename("AdaptiveApertureDiagnostics", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_variable_residuals(save):
    plt.title("Stellar Variability Residuals")
    plt.ylabel("Residuals (flux)")
    plt.xlabel("Time [JD]")
    plt.legend()
    plt.savefig(_working_artifacts_dir(save) / "Variable_Residuals.png")
    plt.close()


def _finite_plot_float(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) else None


def _stellar_variability_reference_label(vsp_param, comparison_label):
    comp_ra = _finite_plot_float(vsp_param.get('comp_ra'))
    comp_dec = _finite_plot_float(vsp_param.get('comp_dec'))
    details = []
    comparison_label = str(comparison_label).strip() if comparison_label else ""

    if comparison_label and not comparison_label.lower().startswith("ra="):
        details.append(f"Label: {comparison_label}")

    if comp_ra is not None and comp_dec is not None:
        details.extend((f"RA={comp_ra:.6f}", f"Dec={comp_dec:.6f}"))

    if not details and comparison_label:
        details.append(comparison_label)

    return "\n".join(details)


def _stellar_variability_comparison_metadata_label(vsp_param):
    details = []
    observed_filter = vsp_param.get('observed_filter')
    if observed_filter:
        details.append(f"Original filter: {observed_filter}")

    comparison_mag = magnitude_text(
        vsp_param.get('catalog_mag_band') or vsp_param.get('mag_band') or 'V',
        vsp_param.get('cmag'),
        vsp_param.get('cmag_err'),
    )
    if comparison_mag is not None:
        details.append(f"Comparison mag: {comparison_mag}")

    return " | ".join(details)


def plot_stellar_variability(vsp_params, save, s_name, vsp_auid_comp):
    if not vsp_params:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted_points = 0
    for vsp_p in vsp_params:
        if not is_usable_apparent_magnitude(vsp_p.get('mag')):
            continue
        mag_err = normalized_magnitude_error(vsp_p.get('mag_err'))
        ax.errorbar(vsp_p['time'], vsp_p['mag'], yerr=mag_err, color="tomato", fmt='.')
        plotted_points += 1

    if plotted_points == 0:
        plt.close(fig)
        return

    first_param = vsp_params[0]
    band = first_param.get('mag_band') or 'V'
    reference_label = _stellar_variability_reference_label(first_param, vsp_auid_comp)
    title_lines = [s_name]
    if reference_label:
        title_lines.append(reference_label)
    metadata_label = _stellar_variability_comparison_metadata_label(first_param)
    if metadata_label:
        title_lines.append(metadata_label)
    ax.set_title("\n".join(title_lines), fontsize=11)
    ax.set_ylabel(f"Magnitude ({band})")
    ax.set_xlabel("Time [JD]")
    fig.tight_layout()
    output_dir = _working_artifacts_dir(save)
    fig.savefig(output_dir / f"Stellar_Variability.png", bbox_inches="tight")
    plt.close(fig)


def _stellar_variability_magnitude_series(vsp_params):
    rows = []
    for vsp_p in vsp_params or []:
        time_value = _finite_plot_float(vsp_p.get('time'))
        mag_value = _finite_plot_float(vsp_p.get('mag'))
        mag_err = normalized_magnitude_error(vsp_p.get('mag_err'))
        if (
            time_value is None
            or mag_value is None
            or mag_err is None
            or not is_usable_apparent_magnitude(mag_value)
        ):
            continue
        rows.append((time_value, mag_value, mag_err, vsp_p))

    if not rows:
        return None

    rows.sort(key=lambda row: row[0])
    times = np.array([row[0] for row in rows], dtype=float)
    magnitudes = np.array([row[1] for row in rows], dtype=float)
    magnitude_errors = np.array([row[2] for row in rows], dtype=float)
    return times, magnitudes, magnitude_errors, rows[0][3]


def _stellar_variability_apparent_magnitude_calibration(fit):
    series = _stellar_variability_magnitude_series(
        getattr(fit, 'stellar_variability_params', None)
    )
    if series is None:
        return None
    _, magnitudes, _, first_param = series
    finite = np.isfinite(magnitudes)
    if not np.any(finite):
        return None
    return {
        'baseline_magnitude': float(np.nanmedian(magnitudes[finite])),
        'band': first_param.get('mag_band') or 'V',
    }


def _add_apparent_magnitude_axis(ax_lc, fit):
    calibration = _stellar_variability_apparent_magnitude_calibration(fit)
    if calibration is None:
        return False
    baseline_magnitude = calibration['baseline_magnitude']

    def flux_to_magnitude(flux):
        flux = np.asarray(flux, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            return baseline_magnitude - (2.5 * np.log10(flux))

    def magnitude_to_flux(magnitude):
        magnitude = np.asarray(magnitude, dtype=float)
        with np.errstate(over='ignore', invalid='ignore'):
            return 10 ** ((baseline_magnitude - magnitude) / 2.5)

    secondary_axis = ax_lc.secondary_yaxis(
        'right',
        functions=(flux_to_magnitude, magnitude_to_flux),
    )
    secondary_axis.set_ylabel(f"Apparent Magnitude ({calibration['band']})")
    return True


# Observation statistics series selection
def _select_plot_rows(rows, sort_index=None, sigma_mask=None, relative_flux_mask=None):
    rows = np.asarray(rows)

    if sort_index is not None:
        rows = rows[np.asarray(sort_index)]

    if sigma_mask is not None:
        sigma_mask = np.asarray(sigma_mask)
        if sigma_mask.dtype == bool and rows.shape[0] == sigma_mask.shape[0]:
            rows = rows[sigma_mask]

    if relative_flux_mask is not None:
        relative_flux_mask = np.asarray(relative_flux_mask)
        if relative_flux_mask.dtype == bool and rows.shape[0] == relative_flux_mask.shape[0]:
            rows = rows[relative_flux_mask]

    return rows


def plot_obs_stats(fit, comp_stars, psf, si, gi, target_name, save, date, relative_flux_mask=None,
                   background_series=None):
    fit_time = np.asarray(fit.time)
    fit_airmass = np.asarray(fit.airmass)
    temp_dir = _working_artifacts_dir(save)

    for i in range(len(comp_stars) + 1):
        if i == 0:
            title, key = target_name, "target"
        else:
            title, key = f"Comp Star {i}", f"comp{i}"

        fig, axs = plt.subplots(3, 2, figsize=(12, 10))
        fig.suptitle(f"Observing Statistics - {title} - {date}")

        star_stats = _select_plot_rows(psf[key], sort_index=si, sigma_mask=gi,
                                       relative_flux_mask=relative_flux_mask)
        background_data = None
        if background_series is not None and key in background_series:
            background_data = _select_plot_rows(
                background_series[key],
                sort_index=si,
                sigma_mask=gi,
                relative_flux_mask=relative_flux_mask,
            )

        plot_len_inputs = [fit_time.shape[0], fit_airmass.shape[0], star_stats.shape[0]]
        if background_data is not None:
            plot_len_inputs.append(background_data.shape[0])
        plot_len = min(plot_len_inputs)
        if plot_len == 0:
            plt.close(fig)
            continue

        time_data = fit_time[:plot_len]
        airmass_data = fit_airmass[:plot_len]
        star_stats = star_stats[:plot_len]
        if background_data is None:
            background_data = star_stats[:, 6]
        else:
            background_data = np.asarray(background_data)[:plot_len]

        axs[0, 0].set(xlabel="Time [BJD_TDB]", ylabel="X-Centroid [px]")
        axs[0, 0].plot(time_data, star_stats[:, 0], 'k.')

        axs[0, 1].set(xlabel="Time [BJD_TDB]", ylabel="Y-Centroid [px]")
        axs[0, 1].plot(time_data, star_stats[:, 1], 'k.')

        axs[1, 0].set(xlabel="Time [BJD_TDB]", ylabel="Seeing [px]")
        axs[1, 0].plot(time_data, 2.355 * 0.5 * (star_stats[:, 3] + star_stats[:, 4]), 'k.')

        axs[1, 1].set(xlabel="Time [BJD_TDB]", ylabel="Airmass")
        axs[1, 1].plot(time_data, airmass_data, 'k.')

        axs[2, 0].set(xlabel="Time [BJD_TDB]", ylabel="Amplitude [ADU]")
        axs[2, 0].plot(time_data, star_stats[:, 2], 'k.')

        axs[2, 1].set(xlabel="Time [BJD_TDB]", ylabel="Background [ADU]")
        axs[2, 1].plot(time_data, background_data, 'k.')

        plt.tight_layout()

        try:
            fig.savefig(temp_dir / _dated_plot_filename("Observing_Statistics", key, date=date, extension="png"), bbox_inches="tight")
            fig.savefig(temp_dir / _dated_plot_filename("Observing_Statistics", key, date=date, extension="pdf"), bbox_inches="tight")
        except Exception:
            pass
        plt.close()


# Plotting Final Lightcurve
def _final_lightcurve_model_grid(fit, high_res):
    if hasattr(fit, 'phase_upsample') and hasattr(fit, 'transit_upsample'):
        x_values = np.asarray(fit.phase_upsample, dtype=float)
        model = np.asarray(fit.transit_upsample, dtype=float)
        times = getattr(fit, 'time_upsample', None)
        if times is not None:
            times = np.asarray(times, dtype=float)
            if times.shape != model.shape:
                times = None
        return x_values, model, times

    phase = np.asarray(getattr(fit, 'phase', np.array([])), dtype=float)
    model = np.asarray(high_res, dtype=float)
    if phase.size == 0 or model.size == 0:
        return None, None, None
    x_values = np.linspace(np.nanmin(phase), np.nanmax(phase), model.size)
    return x_values, model, None


def _transit_model_uncertainty_envelope_for_grid(fit, times, model_shape):
    if times is None or times.shape != model_shape:
        return None

    uncertainty_func = getattr(fit, 'transit_model_uncertainty', None)
    if not callable(uncertainty_func):
        return None

    try:
        envelope = uncertainty_func(times)
    except Exception:
        return None
    if envelope is None or len(envelope) != 2:
        return None

    lower = np.asarray(envelope[0], dtype=float)
    upper = np.asarray(envelope[1], dtype=float)
    if lower.shape != model_shape or upper.shape != model_shape:
        return None
    return lower, upper


def _plot_final_data_scatter_uncertainty_band(ax_lc, fit, high_res):
    empirical_uncertainty = getattr(fit, 'empirical_transit_uncertainty', None)
    if not isinstance(empirical_uncertainty, dict) or not empirical_uncertainty.get('available'):
        empirical_uncertainty = fit_empirical_transit_uncertainty(fit)
    if not isinstance(empirical_uncertainty, dict) or not empirical_uncertainty.get('available'):
        return False

    depth_uncertainty = empirical_uncertainty.get('depth_uncertainty_fraction')
    try:
        depth_uncertainty = float(depth_uncertainty)
    except (TypeError, ValueError):
        return False
    if not np.isfinite(depth_uncertainty) or depth_uncertainty <= 0:
        return False

    x_values, model, times = _final_lightcurve_model_grid(fit, high_res)
    if x_values is None or model is None or x_values.shape != model.shape:
        return False

    finite = np.isfinite(x_values) & np.isfinite(model)
    if not np.any(finite):
        return False

    empirical_lower = model - depth_uncertainty
    empirical_upper = model + depth_uncertainty
    sort_index = np.argsort(x_values)
    x_sorted = x_values[sort_index]
    finite_sorted = finite[sort_index]
    empirical_lower_sorted = empirical_lower[sort_index]
    empirical_upper_sorted = empirical_upper[sort_index]

    model_envelope = _transit_model_uncertainty_envelope_for_grid(fit, times, model.shape)
    drew_band = False
    def next_label():
        nonlocal drew_band
        drew_band = True
        return '_nolegend_'

    if model_envelope is not None:
        model_lower, model_upper = model_envelope
        model_lower_sorted = np.asarray(model_lower, dtype=float)[sort_index]
        model_upper_sorted = np.asarray(model_upper, dtype=float)[sort_index]

        upper_region = (
            finite_sorted
            & np.isfinite(empirical_upper_sorted)
            & np.isfinite(model_upper_sorted)
            & (empirical_upper_sorted > model_upper_sorted)
        )
        lower_region = (
            finite_sorted
            & np.isfinite(empirical_lower_sorted)
            & np.isfinite(model_lower_sorted)
            & (empirical_lower_sorted < model_lower_sorted)
        )
        if np.any(upper_region):
            ax_lc.fill_between(
                x_sorted,
                model_upper_sorted,
                empirical_upper_sorted,
                where=upper_region,
                interpolate=True,
                color='#6a1b9a',
                alpha=0.16,
                linewidth=0,
                zorder=2.35,
                label=next_label(),
            )
        if np.any(lower_region):
            ax_lc.fill_between(
                x_sorted,
                empirical_lower_sorted,
                model_lower_sorted,
                where=lower_region,
                interpolate=True,
                color='#6a1b9a',
                alpha=0.16,
                linewidth=0,
                zorder=2.35,
                label=next_label(),
            )
        return drew_band

    ax_lc.fill_between(
        x_sorted,
        empirical_lower_sorted,
        empirical_upper_sorted,
        where=finite_sorted,
        interpolate=True,
        color='#6a1b9a',
        alpha=0.14,
        linewidth=0,
        zorder=2.0,
        label=next_label(),
    )
    return drew_band


def _plot_final_residual_rejected_points(ax_lc, ax_res, fit):
    rejection = getattr(fit, 'final_residual_rejection', None)
    if not isinstance(rejection, dict) or not rejection.get('applied'):
        return

    phase = np.asarray(rejection.get('rejected_phase', []), dtype=float)
    flux = np.asarray(rejection.get('rejected_flux', []), dtype=float)
    residual_percent = np.asarray(rejection.get('rejected_residual_percent', []), dtype=float)
    plot_count = min(phase.size, flux.size, residual_percent.size)
    if plot_count == 0:
        return

    phase = phase[:plot_count]
    flux = flux[:plot_count]
    residual_percent = residual_percent[:plot_count]
    finite = np.isfinite(phase) & np.isfinite(flux) & np.isfinite(residual_percent)
    if not np.any(finite):
        return

    ax_lc.scatter(
        phase[finite],
        flux[finite],
        marker='x',
        s=58,
        linewidths=1.6,
        color='red',
        zorder=1200,
        label='_nolegend_',
    )
    ax_res.scatter(
        phase[finite],
        residual_percent[finite],
        marker='x',
        s=58,
        linewidths=1.6,
        color='red',
        zorder=1200,
        label='_nolegend_',
    )


def plot_final_lightcurve(fit, high_res, targ_name, save, date):
    if getattr(fit, 'stellar_variability_only', False):
        series = _stellar_variability_magnitude_series(
            getattr(fit, 'stellar_variability_params', None)
        )
        if series is None:
            return

        obs_time, magnitudes, magnitude_errors, first_param = series
        f, ax_lc = plt.subplots(figsize=(8, 5))
        title_name = getattr(fit, 'stellar_variability_target_name', targ_name)
        title_lines = [title_name]
        reference_label = _stellar_variability_reference_label(
            first_param,
            getattr(fit, 'stellar_variability_reference_label', first_param.get('cname')),
        )
        if reference_label:
            title_lines.append(reference_label)
        metadata_label = _stellar_variability_comparison_metadata_label(first_param)
        if metadata_label:
            title_lines.append(metadata_label)

        ax_lc.set_title("\n".join(title_lines), fontsize=11)
        ax_lc.errorbar(
            obs_time,
            magnitudes,
            yerr=magnitude_errors,
            color="tomato",
            fmt='.',
        )
        band = first_param.get('mag_band') or 'V'
        ax_lc.set_ylabel(f"Magnitude ({band})")
        ax_lc.set_xlabel("Time [BJD_TDB]")
        f.tight_layout()

        Path(save).mkdir(parents=True, exist_ok=True)
        try:
            f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="png"), bbox_inches="tight")
            f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="pdf"), bbox_inches="tight")
        except Exception:
            pass
        plt.close(f)
        return

    empirical_uncertainty = getattr(fit, 'empirical_transit_uncertainty', None)
    if not isinstance(empirical_uncertainty, dict) or not empirical_uncertainty.get('available'):
        empirical_uncertainty = fit_empirical_transit_uncertainty(fit)
        if isinstance(empirical_uncertainty, dict) and empirical_uncertainty.get('available'):
            try:
                fit.empirical_transit_uncertainty = empirical_uncertainty
            except Exception:
                pass

    f, (ax_lc, ax_res) = _plot_bestfit_for_lightcurve_png(
        fit,
        show_flux_baseline_label=False,
        show_model_uncertainty=True,
        show_baseline_uncertainty=True,
    )

    ax_lc.set_title(targ_name)
    drew_data_scatter_band = _plot_final_data_scatter_uncertainty_band(ax_lc, fit, high_res)
    if hasattr(fit, 'phase_upsample') and hasattr(fit, 'transit_upsample'):
        ax_lc.plot(fit.phase_upsample, fit.transit_upsample, 'r', zorder=1000, lw=2)
    else:
        ax_lc.plot(np.linspace(np.nanmin(fit.phase), np.nanmax(fit.phase), 1000), high_res, 'r', zorder=1000, lw=2)
    _plot_final_residual_rejected_points(ax_lc, ax_res, fit)
    if drew_data_scatter_band:
        ax_lc.legend(loc='best')
    _add_apparent_magnitude_axis(ax_lc, fit)

    Path(save).mkdir(parents=True, exist_ok=True)
    try:
        f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="png"), bbox_inches="tight")
        f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    except Exception:
        pass
    plt.close()


def _plot_scalar(value, default=np.nan):
    try:
        result = np.asarray(value, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return default
    if result.size == 0:
        return default
    result = float(result[0])
    return result if np.isfinite(result) else default


def _plot_positive_error(value):
    value = _plot_scalar(value)
    if not np.isfinite(value) or value < 0:
        return np.nan
    return value


def _decimal_places_for_two_sigfig_error(error):
    error = _plot_positive_error(error)
    if not np.isfinite(error) or error == 0:
        return None

    exponent = int(np.floor(np.log10(abs(error))))
    return max(0, 1 - exponent)


def _format_parameter_value(value, error=None, unit="", split_error=False):
    value = _plot_scalar(value)
    if not np.isfinite(value):
        return "n/a"

    suffix = f" {unit}" if unit else ""
    if error is None:
        return f"{value:.6f}".rstrip('0').rstrip('.') + suffix

    error = _plot_positive_error(error)
    if np.isfinite(error):
        decimal_places = _decimal_places_for_two_sigfig_error(error)
        if decimal_places is None:
            decimal_places = 0
        value_text = f"{value:.{decimal_places}f}"
        error_text = f"{error:.{decimal_places}f}"
        if split_error:
            return f"{value_text}\n+/- {error_text}{suffix}"
        return f"{value_text} +/- {error_text}{suffix}"
    return f"{value:.6f}".rstrip('0').rstrip('.') + suffix


def _prior_impact_parameter_value_error(planet_dict):
    ars = _plot_scalar(planet_dict.get('aRs'))
    inc = _plot_scalar(planet_dict.get('inc'))
    if not np.isfinite(ars) or not np.isfinite(inc):
        return np.nan, np.nan

    ecc = _plot_scalar(planet_dict.get('ecc'), 0.0)
    omega = np.deg2rad(_plot_scalar(planet_dict.get('omega'), 0.0))
    denominator = 1.0 + ecc * np.sin(omega)
    if not np.isfinite(denominator) or np.isclose(denominator, 0.0):
        return np.nan, np.nan

    scale_factor = (1.0 - ecc ** 2) / denominator
    inc_rad = np.deg2rad(inc)
    impact_parameter = scale_factor * ars * np.cos(inc_rad)

    ars_error = _plot_positive_error(planet_dict.get('aRsUnc'))
    inc_error = _plot_positive_error(planet_dict.get('incUnc'))
    if np.isfinite(ars_error) and np.isfinite(inc_error):
        impact_error = np.hypot(
            scale_factor * np.cos(inc_rad) * ars_error,
            scale_factor * ars * np.sin(inc_rad) * np.deg2rad(inc_error),
        )
    else:
        impact_error = np.nan

    return float(impact_parameter), float(impact_error) if np.isfinite(impact_error) else np.nan


def _ephemeris_prior_at_posterior_epoch(planet_dict, posterior_tmid):
    mid_t = _plot_scalar(planet_dict.get('midT'))
    period = _plot_scalar(planet_dict.get('pPer'))
    posterior_tmid = _plot_scalar(posterior_tmid)
    if not np.isfinite(mid_t) or not np.isfinite(period) or period <= 0 or not np.isfinite(posterior_tmid):
        return mid_t, _plot_positive_error(planet_dict.get('midTUnc')), None

    epoch = int(np.round((posterior_tmid - mid_t) / period))
    expected_tmid = mid_t + epoch * period

    error_terms = []
    mid_t_error = _plot_positive_error(planet_dict.get('midTUnc'))
    period_error = _plot_positive_error(planet_dict.get('pPerUnc'))
    if np.isfinite(mid_t_error):
        error_terms.append(mid_t_error)
    if np.isfinite(period_error):
        error_terms.append(abs(epoch) * period_error)

    if error_terms:
        expected_error = float(np.sqrt(np.sum(np.square(error_terms))))
    else:
        expected_error = np.nan
    return float(expected_tmid), expected_error, epoch


def _posterior_parameter_value_error(fit, parameter_key, empirical_uncertainty):
    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}

    if parameter_key == 'b':
        errors_override = {}
        errors = getattr(fit, 'errors', {}) or {}
        sample_errors = getattr(fit, 'sample_errors', {}) or {}
        b_error = _plot_positive_error(errors.get('b'))
        if not np.isfinite(b_error):
            b_error = _plot_positive_error(sample_errors.get('b'))
        if np.isfinite(b_error):
            errors_override['b'] = float(b_error * empirical_red_noise_error_scale(empirical_uncertainty))
        ars_error = fit_parameter_model_data_uncertainty(
            fit,
            'ars',
            empirical_uncertainty=empirical_uncertainty,
        )
        inc_error = fit_parameter_model_data_uncertainty(
            fit,
            'inc',
            empirical_uncertainty=empirical_uncertainty,
        )
        if np.isfinite(ars_error):
            errors_override['ars'] = ars_error
        if np.isfinite(inc_error):
            errors_override['inc'] = inc_error
        return fit_impact_parameter_value_error(fit, errors_override=errors_override)

    value = _plot_scalar(parameters.get(parameter_key))
    if parameter_key == 'rprs':
        error = _plot_positive_error(
            (empirical_uncertainty or {}).get('combined_rprs_uncertainty')
        )
        if not np.isfinite(error):
            error = _plot_positive_error(errors.get(parameter_key))
        return value, error

    error = fit_parameter_model_data_uncertainty(
        fit,
        parameter_key,
        empirical_uncertainty=empirical_uncertainty,
    )
    if not np.isfinite(error):
        error = _plot_positive_error(errors.get(parameter_key))
    return value, error


def _prior_posterior_comparison_rows(fit, planet_dict):
    empirical_uncertainty = getattr(fit, 'empirical_transit_uncertainty', None)
    if not isinstance(empirical_uncertainty, dict) or not empirical_uncertainty.get('available'):
        empirical_uncertainty = fit_empirical_transit_uncertainty(fit)
        if isinstance(empirical_uncertainty, dict) and empirical_uncertainty.get('available'):
            try:
                fit.empirical_transit_uncertainty = empirical_uncertainty
            except Exception:
                pass

    definitions = [
        ("Tmid", "tmid", "midT", "midTUnc", "", True),
        ("Rp/R*", "rprs", "rprs", "rprsUnc", "", False),
        ("a/Rs", "ars", "aRs", "aRsUnc", "", False),
        ("Inc.", "inc", "inc", "incUnc", "deg", False),
        ("b", "b", None, None, "", False),
    ]

    rows = []
    rprs_prior_fallback = bool(
        getattr(fit, 'rprs_prior_fallback_applied', False)
        or (isinstance(empirical_uncertainty, dict)
            and empirical_uncertainty.get('rprs_prior_fallback_applied'))
        or (isinstance(empirical_uncertainty, dict)
            and empirical_uncertainty.get('rprs_uncertainty_basis') == 'prior_assumed_data_only')
    )
    omitted_notes = []
    for label, parameter_key, prior_key, prior_error_key, unit, split_error in definitions:
        if parameter_key == 'rprs' and rprs_prior_fallback:
            prior_value = _plot_scalar(planet_dict.get(prior_key))
            prior_error = _plot_positive_error(planet_dict.get(prior_error_key))
            posterior_value, posterior_error = _posterior_parameter_value_error(
                fit,
                parameter_key,
                empirical_uncertainty,
            )
            omitted_notes.append(
                "Rp/R* omitted: prior value assumed, not measured "
                f"({_format_parameter_value(prior_value, prior_error)}; "
                f"data-only uncertainty {_format_parameter_value(posterior_value, posterior_error)})."
            )
            continue

        posterior_value, posterior_error = _posterior_parameter_value_error(
            fit,
            parameter_key,
            empirical_uncertainty,
        )

        if parameter_key == 'b':
            prior_value, prior_error = _prior_impact_parameter_value_error(planet_dict)
            prior_label = "Prior"
        elif parameter_key == 'tmid':
            prior_value, prior_error, _ = _ephemeris_prior_at_posterior_epoch(
                planet_dict,
                posterior_value,
            )
            prior_label = "Prior"
        else:
            prior_value = _plot_scalar(planet_dict.get(prior_key))
            prior_error = _plot_positive_error(planet_dict.get(prior_error_key))
            prior_label = "Prior"

        if not np.isfinite(prior_value) or not np.isfinite(posterior_value):
            continue

        error_terms = [
            term for term in (prior_error, posterior_error)
            if np.isfinite(term) and term > 0
        ]
        if error_terms:
            combined_sigma = float(np.sqrt(np.sum(np.square(error_terms))))
        else:
            separation = abs(posterior_value - prior_value)
            combined_sigma = float(separation) if separation > 0 else np.nan
        if not np.isfinite(combined_sigma) or combined_sigma <= 0:
            continue

        posterior_offset = (posterior_value - prior_value) / combined_sigma
        prior_error_sigma = prior_error / combined_sigma if np.isfinite(prior_error) else 0.0
        posterior_error_sigma = (
            posterior_error / combined_sigma if np.isfinite(posterior_error) else 0.0
        )
        prior_assumed = parameter_key == 'rprs' and rprs_prior_fallback

        rows.append({
            "label": label,
            "parameter_key": parameter_key,
            "posterior_offset": float(posterior_offset),
            "prior_error_sigma": float(prior_error_sigma),
            "posterior_error_sigma": float(posterior_error_sigma),
            "prior_text": _format_parameter_value(
                prior_value,
                prior_error,
                unit=unit,
                split_error=split_error,
            ),
            "posterior_text": _format_parameter_value(
                posterior_value,
                posterior_error,
                unit=unit,
                split_error=split_error,
            ),
            "prior_label": prior_label,
            "prior_assumed": prior_assumed,
        })

    return rows, omitted_notes


def plot_prior_posterior_comparison(fit, planet_dict, targ_name, save, date):
    rows, omitted_notes = _prior_posterior_comparison_rows(fit, planet_dict)
    if not rows:
        return None

    note_height = 0.34 * len(omitted_notes)
    row_spacing = 1.35
    fig_height = max(5.0, 1.02 * len(rows) + 2.0 + note_height)
    fig, ax = plt.subplots(figsize=(11.8, fig_height))

    y_positions = np.arange(len(rows), dtype=float) * row_spacing
    posterior_offsets = np.array([row["posterior_offset"] for row in rows], dtype=float)
    prior_errors = np.array([row["prior_error_sigma"] for row in rows], dtype=float)
    posterior_errors = np.array([row["posterior_error_sigma"] for row in rows], dtype=float)

    xmin = min(-3.5, np.nanmin(np.r_[posterior_offsets - posterior_errors, -prior_errors]) - 0.45)
    xmax = max(3.5, np.nanmax(np.r_[posterior_offsets + posterior_errors, prior_errors]) + 0.45)

    ax.axvspan(-1.0, 1.0, color='#2e7d32', alpha=0.08, linewidth=0)
    ax.axvspan(-3.0, 3.0, color='#f9a825', alpha=0.06, linewidth=0)
    ax.axvline(0.0, color='0.25', lw=1.2, ls='--', zorder=1)

    ax.errorbar(
        np.zeros_like(y_positions),
        y_positions + 0.13,
        xerr=prior_errors,
        fmt='o',
        ms=6,
        color='#1565c0',
        ecolor='#1565c0',
        elinewidth=1.4,
        capsize=3,
        label='Prior',
        zorder=5,
    )
    ax.errorbar(
        posterior_offsets,
        y_positions - 0.13,
        xerr=posterior_errors,
        fmt='s',
        ms=6,
        color='#c62828',
        ecolor='#c62828',
        elinewidth=1.4,
        capsize=3,
        label='Posterior',
        zorder=6,
    )

    for y_position, row in zip(y_positions, rows):
        annotation = (
            f"{row['prior_label']}\n"
            f"{row['prior_text']}\n"
            "Posterior\n"
            f"{row['posterior_text']}"
        )
        ax.text(
            1.015,
            y_position,
            annotation,
            transform=ax.get_yaxis_transform(),
            ha='left',
            va='center',
            fontsize=8.5,
            linespacing=1.12,
            color='0.18',
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels([row["label"] for row in rows])
    ax.invert_yaxis()
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Posterior offset from prior [combined sigma]")
    ax.set_title(f"{targ_name} Prior vs Posterior Transit Parameters")
    ax.grid(axis='x', alpha=0.28)
    if omitted_notes:
        ax.text(
            0.0,
            -0.16,
            "\n".join(omitted_notes),
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=9,
            color='0.22',
        )
    ax.legend(
        handles=[
            Line2D([0], [0], marker='o', color='none', markerfacecolor='#1565c0',
                   markeredgecolor='#1565c0', markersize=7, label='Prior'),
            Line2D([0], [0], marker='s', color='none', markerfacecolor='#c62828',
                   markeredgecolor='#c62828', markersize=7, label='Posterior'),
        ],
        loc='lower right',
    )
    fig.subplots_adjust(right=0.64)

    Path(save).mkdir(parents=True, exist_ok=True)
    png_path = Path(save) / _dated_plot_filename(
        "PriorPosteriorComparison",
        targ_name,
        date=date,
        extension="png",
    )
    pdf_path = Path(save) / _dated_plot_filename(
        "PriorPosteriorComparison",
        targ_name,
        date=date,
        extension="pdf",
    )
    try:
        fig.savefig(png_path, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
    except Exception:
        png_path = None
    plt.close(fig)
    return png_path


def _fit_ktmf_metric_contributions_status(fit):
    transit_qc = getattr(fit, 'transit_qc', None)
    if not isinstance(transit_qc, dict):
        transit_qc = {}

    metric = _plot_scalar(
        getattr(fit, 'transit_qc_ktmf_metric', transit_qc.get('ktmf_metric', np.nan))
    )
    contributions = getattr(fit, 'transit_qc_ktmf_contributions', None)
    if not contributions:
        contributions = transit_qc.get('ktmf_contributions', [])
    if not isinstance(contributions, (list, tuple)):
        contributions = []

    status = _ktmf_status_from_metric(metric)
    if not status:
        status = getattr(fit, 'transit_qc_status', transit_qc.get('status', None))
    return metric, list(contributions), status


def _ktmf_status_from_metric(metric):
    metric = _plot_scalar(metric)
    if not np.isfinite(metric):
        return None
    if metric >= 4.0:
        return "pass"
    if metric >= 3.0:
        return "marginal"
    return "fail"


def _short_ktmf_label(label):
    replacements = {
        "Deviation From Expected Value": "Expected Rp/R*",
        "Residual Scatter Around Full Model Fit": "Residual scatter",
        "Residual Flatness": "Residual flatness",
        "Tmid Posterior Gaussianity": "Tmid Gaussianity",
        "Duration Consistency": "Duration",
        "EEBLS Depth SNR": "EEBLS SNR",
        "Sampling / Cadence": "Sampling",
    }
    return replacements.get(str(label), str(label))


def _ktmf_marker_color(score):
    score = _plot_scalar(score)
    if not np.isfinite(score):
        return '0.45'
    if score >= 0.8:
        return '#2e7d32'
    if score >= 0.6:
        return '#f9a825'
    return '#c62828'


def _format_ktmf_metric(value, maximum=5.0):
    value = _plot_scalar(value)
    maximum = _plot_scalar(maximum)
    if not np.isfinite(value):
        return "n/a"
    if np.isfinite(maximum) and maximum > 0:
        return f"{value:.2f} / {maximum:.2f}"
    return f"{value:.2f}"


def _format_ktmf_component_annotation(row):
    if row.get('kind') == 'total':
        status = row.get('status')
        status_text = f"\n{status.upper()}" if status else ""
        uncertainty = _plot_positive_error(row.get('score_uncertainty'))
        uncertainty_text = f"\nscore spread +/- {uncertainty:.2f}" if np.isfinite(uncertainty) else ""
        return f"KTMF\n{_format_ktmf_metric(row.get('points'), row.get('max_points'))}{status_text}{uncertainty_text}"

    if not row.get('available', True):
        return "Not scored"

    score_uncertainty = _plot_positive_error(row.get('score_uncertainty'))
    if np.isfinite(score_uncertainty):
        score_text = _format_parameter_value(row.get('score'), score_uncertainty)
    else:
        score_text = _format_parameter_value(row.get('score'))
    return (
        f"Score\n{score_text}\n"
        f"Points\n{_format_ktmf_metric(row.get('points'), row.get('max_points'))}"
    )


def _compact_ktmf_detail(row):
    detail = row.get('detail')
    if not detail:
        return None
    detail = str(detail)
    if row.get('label') == "Expected Rp/R*":
        if "fixed to the input prior" in detail or "prior" in detail.lower():
            return "Rp/R* prior assumed; not scored."
        keep = []
        for part in detail.split(','):
            part = part.strip()
            if part.startswith("Rp/R* sigma="):
                keep.append(part)
        return ", ".join(keep) if keep else detail
    return detail


def _ktmf_plot_rows(fit):
    metric, contributions, status = _fit_ktmf_metric_contributions_status(fit)
    component_rows = []
    for contribution in contributions:
        if not isinstance(contribution, dict):
            continue
        available = bool(contribution.get('available', True))
        score = _plot_scalar(contribution.get('score'))
        if not available or not np.isfinite(score):
            score = np.nan
        component_rows.append({
            "kind": "component",
            "label": _short_ktmf_label(contribution.get('label', 'KTMF component')),
            "score": float(np.clip(score, 0.0, 1.0)) if np.isfinite(score) else np.nan,
            "score_uncertainty": _plot_positive_error(contribution.get('score_uncertainty')),
            "points": _plot_scalar(contribution.get('points'), 0.0),
            "max_points": _plot_scalar(contribution.get('max_points'), 0.0),
            "available": available and np.isfinite(score),
            "detail": contribution.get('detail'),
        })

    rows = []
    if np.isfinite(metric):
        total_score = float(np.clip(metric / 5.0, 0.0, 1.0))
        available_component_rows = [
            row for row in component_rows
            if row.get('available')
            and np.isfinite(row.get('score', np.nan))
            and np.isfinite(row.get('max_points', np.nan))
            and row.get('max_points', 0.0) > 0
        ]
        score_uncertainty = np.nan
        if len(available_component_rows) > 1:
            scores = np.asarray([row['score'] for row in available_component_rows], dtype=float)
            weights = np.asarray([row['max_points'] for row in available_component_rows], dtype=float)
            if np.isfinite(weights).all() and np.sum(weights) > 0:
                score_uncertainty = float(
                    np.sqrt(np.average((scores - total_score) ** 2, weights=weights))
                )
        rows.append({
            "kind": "total",
            "label": "KTMF total",
            "score": total_score,
            "score_uncertainty": score_uncertainty,
            "points": float(metric),
            "max_points": 5.0,
            "available": True,
            "status": status,
        })
    rows.extend(component_rows)
    return rows


def _score_errorbar_limits(score, uncertainty):
    score = _plot_scalar(score)
    uncertainty = _plot_positive_error(uncertainty)
    if not np.isfinite(score) or not np.isfinite(uncertainty) or uncertainty <= 0:
        return None
    lower = min(uncertainty, max(score, 0.0))
    upper = min(uncertainty, max(1.0 - score, 0.0))
    if lower <= 0 and upper <= 0:
        return None
    return np.asarray([[lower], [upper]], dtype=float)


def _draw_ktmf_score_background(ax):
    ax.axvspan(0.0, 0.6, color='#c62828', alpha=0.055, linewidth=0)
    ax.axvspan(0.6, 0.8, color='#f9a825', alpha=0.09, linewidth=0)
    ax.axvspan(0.8, 1.0, color='#2e7d32', alpha=0.08, linewidth=0)
    ax.axvline(0.6, color='0.55', lw=1.0, ls=':', zorder=1)
    ax.axvline(0.8, color='0.45', lw=1.1, ls='--', zorder=1)
    ax.grid(axis='x', alpha=0.28)


def _plot_ktmf_score_marker(ax, row, y_position):
    score = _plot_scalar(row.get('score'))
    color = _ktmf_marker_color(score)
    marker = 'D' if row.get('kind') == 'total' else 's'
    marker_size = 62 if row.get('kind') == 'total' else 48
    marker_scale = 3.0
    if np.isfinite(score):
        ax.errorbar(
            [score],
            [y_position],
            xerr=_score_errorbar_limits(score, row.get('score_uncertainty')),
            fmt=marker,
            ms=np.sqrt(marker_size) * marker_scale,
            color=color,
            ecolor=color,
            elinewidth=1.8,
            capsize=4,
            markeredgecolor='white',
            markeredgewidth=1.2,
            zorder=5,
        )
    else:
        ax.scatter(
            [0.0],
            [y_position],
            marker='x',
            s=52 * marker_scale ** 2,
            color='0.45',
            linewidths=2.0,
            zorder=5,
        )


def plot_ktmf_qc_metrics(fit, targ_name, save, date):
    rows = _ktmf_plot_rows(fit)
    if not rows:
        return None

    total_rows = [row for row in rows if row.get('kind') == 'total']
    component_rows = [row for row in rows if row.get('kind') != 'total']
    row_spacing = 1.35
    component_height = max(3.6, 0.98 * max(len(component_rows), 1) + 1.3)
    fig_height = component_height + (1.55 if total_rows else 0.0)
    if total_rows:
        fig, (ax_total, ax_components) = plt.subplots(
            2,
            1,
            figsize=(11.8, fig_height),
            sharex=True,
            gridspec_kw={'height_ratios': [1.0, component_height]},
        )
        axes = [ax_total, ax_components]
    else:
        fig, ax_components = plt.subplots(figsize=(11.8, fig_height))
        ax_total = None
        axes = [ax_components]

    for axis in axes:
        _draw_ktmf_score_background(axis)
        axis.set_xlim(-0.05, 1.05)

    if total_rows:
        total_row = total_rows[0]
        _plot_ktmf_score_marker(ax_total, total_row, 0.0)
        ax_total.text(
            1.025,
            0.0,
            _format_ktmf_component_annotation(total_row),
            transform=ax_total.get_yaxis_transform(),
            ha='left',
            va='center',
            fontsize=8.5,
            linespacing=1.12,
            color='0.18',
        )
        ax_total.set_yticks([0.0])
        ax_total.set_yticklabels([total_row["label"]])
        ax_total.set_ylim(0.65, -0.65)
        ax_total.tick_params(axis='x', labelbottom=False)
        ax_total.set_title(f"{targ_name} KTMF QC Metrics")

    component_positions = np.arange(len(component_rows), dtype=float) * row_spacing
    for y_position, row in zip(component_positions, component_rows):
        _plot_ktmf_score_marker(ax_components, row, y_position)
        ax_components.text(
            1.025,
            y_position,
            _format_ktmf_component_annotation(row),
            transform=ax_components.get_yaxis_transform(),
            ha='left',
            va='center',
            fontsize=8.5,
            linespacing=1.12,
            color='0.18',
        )

    ax_components.set_yticks(component_positions)
    ax_components.set_yticklabels([row["label"] for row in component_rows])
    ax_components.invert_yaxis()
    ax_components.set_xlabel("KTMF component score fraction")
    if not total_rows:
        ax_components.set_title(f"{targ_name} KTMF QC Metrics")
    fig.subplots_adjust(right=0.62, hspace=0.12)

    Path(save).mkdir(parents=True, exist_ok=True)
    png_path = Path(save) / _dated_plot_filename(
        "KTMF_QC",
        targ_name,
        date=date,
        extension="png",
    )
    pdf_path = Path(save) / _dated_plot_filename(
        "KTMF_QC",
        targ_name,
        date=date,
        extension="pdf",
    )
    try:
        fig.savefig(png_path, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
    except Exception:
        png_path = None
    plt.close(fig)
    return png_path
