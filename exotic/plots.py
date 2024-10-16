from astropy.visualization import astropy_mpl_style, ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch, SquaredStretch, SqrtStretch, LogStretch
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path

plt.style.use(astropy_mpl_style)


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
    plt.savefig(Path(save) / "temp" / f"CentroidPositions&Distances_{target_name}_{date}.pdf")
    plt.close()

def plot_fov(aper, annulus, sigma, x_targ, y_targ, x_ref, y_ref, image, image_scale, targ_name, save, date, opt_method, min_aper_fov, min_annulus_fov):
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

        # Create the target and reference circles
        target_circle = plt.Circle((x_targ, y_targ), aper, color='r', fill=False, ls='-')
        target_circle_sky = plt.Circle((x_targ, y_targ), aper + annulus, color=outer_circle_color, fill=False, ls='-')

        if aper >= 0:
            ref_circle = plt.Circle((x_ref, y_ref), aper, color='r', fill=False, ls='-')
            ref_circle_sky = plt.Circle((x_ref, y_ref), aper + annulus, color=outer_circle_color, fill=False, ls='-')

        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(image)

        norm = ImageNormalize(image, interval=interval, stretch=stretch, vmin=vmin, vmax=vmax)

        im = plt.imshow(image, norm=norm, origin='lower', cmap='Greys_r', interpolation=None)
        fig.colorbar(im)

        ax.add_artist(target_circle)
        ax.add_artist(target_circle_sky)
        ax.text(x_targ + aper + annulus + 5, y_targ, targ_name, color='w', fontsize=10,
                path_effects=[path_effects.withStroke(linewidth=2, foreground='black')])

        if aper >= 0:
            ax.add_artist(ref_circle)
            ax.add_artist(ref_circle_sky)
            ax.text(x_ref + aper + annulus + 5, y_ref, 'Comp Star', color='w', fontsize=10,
                    path_effects=[path_effects.withStroke(linewidth=2, foreground='black')])

        handles = []
        label_aper = f"{opt_method} Photometry\n(Min Aper: {min_aper_fov:.2f} px)\n(Min Annulus: {min_annulus_fov:.2f} px)"
        
        if opt_method == "Aperture":
            aperture_line = Line2D([], [], color='r', linestyle='-', label=label_aper)
            handles.append(aperture_line)
        elif opt_method == "PSF":
            psf_line = Line2D([], [], color='lime', linestyle='-', label=label_aper)
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

        plt.savefig(Path(save) / "temp" / f"FOV_{targ_name}_{date}_"
                    f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.pdf", bbox_inches='tight')
        plt.savefig(Path(save) / "temp" / f"FOV_{targ_name}_{date}_"
                    f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.png", bbox_inches='tight')
        plt.close()


def plot_flux(times, targ, targ_unc, ref, ref_unc, norm_flux, norm_unc, airmass, targ_name, save, date):
    plt.figure()
    plt.title(f"{targ_name} Raw Flux Values {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Flux [ADU]")
    plt.errorbar(times, targ, yerr=targ_unc, linestyle='None', fmt='-o')
    plt.savefig(Path(save) / "temp" / f"TargetRawFlux_{targ_name}_{date}.pdf")
    plt.close()

    plt.figure()
    plt.title(f"Comparison Star Raw Flux Values {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Flux [ADU]")
    plt.errorbar(times, ref, yerr=ref_unc, linestyle='None', fmt='-o')
    plt.savefig(Path(save) / "temp" / f"CompRawFlux_{targ_name}_{date}.pdf")
    plt.close()

    # Plots final reduced light curve (after the 3 sigma clip)
    plt.figure()
    plt.title(f"{targ_name} Normalized Flux vs. Time {date}")
    plt.xlabel("Time [BJD_TDB]")
    plt.ylabel("Normalized Flux")
    plt.errorbar(times, norm_flux, yerr=norm_unc, linestyle='None', fmt='-bo')
    plt.savefig(Path(save) / "temp" / f"NormalizedFluxTime_{targ_name}_{date}.pdf")
    plt.close()

    # Save normalized flux to text file prior to NS
    params_file = Path(save) / "temp" / f"NormalizedFlux_{targ_name}_{date}.txt"
    with params_file.open('w') as f:
        f.write("BJD,Norm Flux,Norm Err,AM\n")

        for ti, fi, erri, ami in zip(times, norm_flux, norm_unc, airmass):
            f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")


def plot_variable_residuals(save):
    plt.title("Stellar Variability Residuals")
    plt.ylabel("Residuals (flux)")
    plt.xlabel("Time [JD]")
    plt.legend()
    plt.savefig(Path(save) / "temp" / f"Variable_Residuals.png")
    plt.close()


def plot_stellar_variability(vsp_params, save, s_name, vsp_auid_comp):
    for vsp_p in vsp_params:
        plt.errorbar(vsp_p['time'], vsp_p['mag'], yerr=vsp_p['mag_err'], color="tomato", fmt='.')

    plt.title(f"{s_name} (Label: {vsp_auid_comp})")
    plt.ylabel("Vmag")
    plt.xlabel("Time [JD]")
    plt.savefig(Path(save) / "temp" / f"Stellar_Variability.png")
    plt.close()


# Observation statistics from PSF data
def plot_obs_stats(fit, comp_stars, psf, si, gi, target_name, save, date):
    for i in range(len(comp_stars) + 1):
        if i == 0:
            title, key = target_name, "target"
        else:
            title, key = f"Comp Star {i}", f"comp{i}"

        fig, axs = plt.subplots(3, 2, figsize=(12, 10))
        fig.suptitle(f"Observing Statistics - {title} - {date}")

        axs[0, 0].set(xlabel="Time [BJD_TDB]", ylabel="X-Centroid [px]")
        axs[0, 0].plot(fit.time, psf[key][si, 0][gi], 'k.')

        axs[0, 1].set(xlabel="Time [BJD_TDB]", ylabel="Y-Centroid [px]")
        axs[0, 1].plot(fit.time, psf[key][si, 1][gi], 'k.')

        axs[1, 0].set(xlabel="Time [BJD_TDB]", ylabel="Seeing [px]")
        axs[1, 0].plot(fit.time, 2.355 * 0.5 * (psf[key][si, 3][gi] + psf[key][si, 4][gi]), 'k.')

        axs[1, 1].set(xlabel="Time [BJD_TDB]", ylabel="Airmass")
        axs[1, 1].plot(fit.time, fit.airmass, 'k.')

        axs[2, 0].set(xlabel="Time [BJD_TDB]", ylabel="Amplitude [ADU]")
        axs[2, 0].plot(fit.time, psf[key][si, 2][gi], 'k.')

        axs[2, 1].set(xlabel="Time [BJD_TDB]", ylabel="Background [ADU]")
        axs[2, 1].plot(fit.time, psf[key][si, 6][gi], 'k.')

        plt.tight_layout()

        try:
            fig.savefig(Path(save) / "temp" / f"Observing_Statistics_{key}_{date}.png", bbox_inches="tight")
            fig.savefig(Path(save) / "temp" / f"Observing_Statistics_{key}_{date}.pdf", bbox_inches="tight")
        except Exception:
            pass
        plt.close()


# Plotting Final Lightcurve
def plot_final_lightcurve(fit, high_res, targ_name, save, date):
    f, (ax_lc, ax_res) = fit.plot_bestfit()

    ax_lc.set_title(targ_name)
    ax_lc.plot(np.linspace(np.nanmin(fit.phase), np.nanmax(fit.phase), 1000), high_res, 'r', zorder=1000, lw=2)

    Path(save).mkdir(parents=True, exist_ok=True)
    try:
        f.savefig(Path(save) / f"FinalLightCurve_{targ_name}_{date}.png", bbox_inches="tight")
        f.savefig(Path(save) / f"FinalLightCurve_{targ_name}_{date}.pdf", bbox_inches="tight")
    except Exception:
        pass
    plt.close()
