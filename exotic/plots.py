from astropy.visualization import astropy_mpl_style, ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch, SquaredStretch, SqrtStretch, LogStretch
import inspect
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path

try:
    from utils import filename_date_token, safe_output_filename
except ImportError:
    from .utils import filename_date_token, safe_output_filename

plt.style.use(astropy_mpl_style)


def _dated_plot_filename(prefix, *parts, date, extension):
    return safe_output_filename(prefix, *parts, filename_date_token(date), extension=extension)


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
    plt.savefig(Path(save) / "temp" / _dated_plot_filename(
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
                local_sky_inner_radius = max(local_sky_inner_radius, 2.0 * 2.355 * float(sigma))
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
        Path(save, "temp").mkdir(parents=True, exist_ok=True)

        stretch_name = str(stretch.__class__).split('.')[-1].split(apos)[0]
        plt.savefig(Path(save) / "temp" / _dated_plot_filename(
            "FOV",
            targ_name,
            stretch_name,
            date=date,
            extension="pdf",
        ), bbox_inches='tight')
        plt.savefig(Path(save) / "temp" / _dated_plot_filename(
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
    plt.savefig(Path(save) / "temp" / _dated_plot_filename("TargetRawFlux", targ_name, date=date, extension="pdf"))
    plt.close()

    plt.figure()
    plt.title(f"Comparison Star Raw Flux Values {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Flux [ADU]")
    plt.errorbar(times, ref, yerr=ref_unc, linestyle='None', fmt='-o')
    plt.savefig(Path(save) / "temp" / _dated_plot_filename("CompRawFlux", targ_name, date=date, extension="pdf"))
    plt.close()

    # Plots final reduced light curve (after the 3 sigma clip)
    plt.figure()
    plt.title(f"{targ_name} Normalized Flux vs. Time {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Normalized Flux")
    plt.errorbar(times, norm_flux, yerr=norm_unc, linestyle='None', fmt='-bo')
    plt.savefig(Path(save) / "temp" / _dated_plot_filename("NormalizedFluxTime", targ_name, date=date, extension="pdf"))
    plt.close()

    # Save normalized flux to text file prior to NS
    params_file = Path(save) / "temp" / _dated_plot_filename("NormalizedFlux", targ_name, date=date, extension="txt")
    with params_file.open('w') as f:
        f.write("BJD,Norm Flux,Norm Err,AM\n")

        for ti, fi, erri, ami in zip(times, norm_flux, norm_unc, airmass):
            f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")


def plot_comp_star_pairwise_matrix(pairwise_matrix, best_comp_index, targ_name, save, date, method_label):
    matrix = np.asarray(pairwise_matrix, dtype=float)
    if matrix.size == 0:
        return

    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

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
    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
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
    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
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

    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

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
    pairwise_series = summary.get('pairwise_ratio_series', {})
    for color_index, (other_label, ratio_series) in enumerate(pairwise_series.items()):
        ratio_series = np.asarray(ratio_series, dtype=float)
        valid = np.isfinite(times) & np.isfinite(ratio_series)
        if np.any(valid):
            axis.plot(times[valid], ratio_series[valid], color=colors[color_index % len(colors)],
                      alpha=0.55, lw=1.0, label=other_label)

    ensemble_ratio = np.asarray(summary.get('ensemble_ratio_series'), dtype=float)
    ensemble_valid = np.isfinite(times) & np.isfinite(ensemble_ratio)
    if np.any(ensemble_valid):
        axis.plot(times[ensemble_valid], ensemble_ratio[ensemble_valid], color='black', lw=1.8,
                  label='Ensemble')

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

    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

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

    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

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
    plt.savefig(Path(save) / "temp" / f"Variable_Residuals.png")
    plt.close()


def _finite_plot_float(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) else None


def _stellar_variability_reference_label(vsp_param, comparison_label):
    band = vsp_param.get('mag_band') or 'V'
    observed_filter = vsp_param.get('observed_filter')
    cmag = _finite_plot_float(vsp_param.get('cmag'))
    cmag_err = _finite_plot_float(vsp_param.get('cmag_err'))
    comp_ra = _finite_plot_float(vsp_param.get('comp_ra'))
    comp_dec = _finite_plot_float(vsp_param.get('comp_dec'))

    comparison_parts = []
    if vsp_param.get('is_aavso_vsp', True) and comparison_label:
        comparison_parts.append(f"Label={comparison_label}")
    if comp_ra is not None and comp_dec is not None:
        comparison_parts.append(f"RA={comp_ra:.7f}")
        comparison_parts.append(f"Dec={comp_dec:.7f}")
    elif comparison_label:
        comparison_parts.append(str(comparison_label))

    detail_parts = []
    if observed_filter not in (None, ''):
        detail_parts.append(f"Observed filter={observed_filter}")

    if cmag is not None and cmag_err is not None:
        detail_parts.append(f"{band}={cmag:.5f} +/- {cmag_err:.5f}")
    elif cmag is not None:
        detail_parts.append(f"{band}={cmag:.5f}")
    else:
        detail_parts.append(f"{band}=na")

    label_lines = []
    if comparison_parts:
        label_lines.append(", ".join(comparison_parts))
    if detail_parts:
        label_lines.append(", ".join(detail_parts))

    return "\n".join(label_lines)


def plot_stellar_variability(vsp_params, save, s_name, vsp_auid_comp):
    if not vsp_params:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    for vsp_p in vsp_params:
        ax.errorbar(vsp_p['time'], vsp_p['mag'], yerr=vsp_p['mag_err'], color="tomato", fmt='.')

    first_param = vsp_params[0]
    band = first_param.get('mag_band') or 'V'
    reference_label = _stellar_variability_reference_label(first_param, vsp_auid_comp)
    ax.set_title(f"{s_name}\nComparison: {reference_label}")
    ax.set_ylabel(f"Magnitude ({band})")
    ax.set_xlabel("Time [JD]")
    fig.tight_layout()
    output_dir = Path(save) / "temp"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"Stellar_Variability.png")
    plt.close(fig)


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
    temp_dir = Path(save) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

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
def plot_final_lightcurve(fit, high_res, targ_name, save, date):
    f, (ax_lc, ax_res) = _plot_bestfit_for_lightcurve_png(
        fit,
        show_flux_baseline_label=False,
        show_model_uncertainty=True,
        show_baseline_uncertainty=True,
    )

    ax_lc.set_title(targ_name)
    if hasattr(fit, 'phase_upsample') and hasattr(fit, 'transit_upsample'):
        ax_lc.plot(fit.phase_upsample, fit.transit_upsample, 'r', zorder=1000, lw=2)
    else:
        ax_lc.plot(np.linspace(np.nanmin(fit.phase), np.nanmax(fit.phase), 1000), high_res, 'r', zorder=1000, lw=2)

    Path(save).mkdir(parents=True, exist_ok=True)
    try:
        f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="png"), bbox_inches="tight")
        f.savefig(Path(save) / _dated_plot_filename("FinalLightCurve", targ_name, date=date, extension="pdf"), bbox_inches="tight")
    except Exception:
        pass
    plt.close()
