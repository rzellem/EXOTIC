from astropy.visualization import astropy_mpl_style, ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch, SquaredStretch, SqrtStretch, LogStretch
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.ndimage import median_filter

plt.style.use(astropy_mpl_style)


# Plots of the centroid positions as a function of time
def plot_centroids(x_targ, y_targ, x_ref, y_ref, times, target_name, save, date):
    fig, axs = plt.subplots(3, 2, figsize=(12, 10))

    axs[0, 0].set_title(f"{target_name} X Centroid Position {date}", fontsize=14)
    axs[0, 0].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[0, 0].set_ylabel("X Pixel Position", fontsize=12)
    axs[0, 0].plot(times - np.nanmin(times), x_targ, '-bo')

    axs[0, 1].set_title(f"{target_name} Y Centroid Position {date}", fontsize=14)
    axs[0, 1].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[0, 1].set_ylabel("Y Pixel Position", fontsize=12)
    axs[0, 1].plot(times - np.nanmin(times), y_targ, '-bo')

    axs[1, 0].set_title(f"Comp Star X Centroid Position {date}", fontsize=14)
    axs[1, 0].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[1, 0].set_ylabel("X Pixel Position", fontsize=12)
    axs[1, 0].plot(times - np.nanmin(times), x_ref, '-ro')

    axs[1, 1].set_title(f"Comp Star Y Centroid Position {date}", fontsize=14)
    axs[1, 1].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[1, 1].set_ylabel("Y Pixel Position", fontsize=12)
    axs[1, 1].plot(times - np.nanmin(times), y_ref, '-ro')

    axs[2, 0].set_title("Distance between Target and Comparison X position", fontsize=14)
    axs[2, 0].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[2, 0].set_ylabel("X Pixel Distance", fontsize=12)
    for e in range(len(x_targ)):
        axs[2, 0].plot(times[e] - np.nanmin(times), abs(x_targ[e] - x_ref[e]), 'bo')

    axs[2, 1].set_title("Distance between Target and Comparison Y position", fontsize=14)
    axs[2, 1].set_xlabel(f"Time (JD-{np.nanmin(times)})", fontsize=12)
    axs[2, 1].set_ylabel("Y Pixel Distance", fontsize=12)
    for e in range(len(y_targ)):
        axs[2, 1].plot(times[e] - np.nanmin(times), abs(y_targ[e] - y_ref[e]), 'bo')

    plt.tight_layout()
    plt.savefig(Path(save) / "temp" / f"CentroidPositions&Distances_{target_name}_{date}.png")
    plt.close()


def plot_fov(aper, annulus, sigma, x_targ, y_targ, x_ref, y_ref, image, image_scale, targ_name, save, date):
    ref_circle, ref_circle_sky = None, None
    picframe = 10 * (aper + 15 * sigma)

    pltx = [max([0, min([x_targ, x_ref]) - picframe]), min([np.shape(image)[1], max([x_targ, x_ref]) + picframe])]
    plty = [max([0, min([y_targ, y_ref]) - picframe]), min([np.shape(image)[0], max([y_targ, y_ref]) + picframe])]

    for stretch in [LinearStretch(), SquaredStretch(), SqrtStretch(), LogStretch()]:
        fig, ax = plt.subplots()

        # Draw apertures and sky annuli
        target_circle = plt.Circle((x_targ, y_targ), aper, color='lime', fill=False, ls='-', label='Target')
        target_circle_sky = plt.Circle((x_targ, y_targ), aper + annulus, color='lime', fill=False, ls='--', lw=.5)
        if aper >= 0:
            ref_circle = plt.Circle((x_ref, y_ref), aper, color='r', fill=False, ls='-.', label='Comp')
            ref_circle_sky = plt.Circle((x_ref, y_ref), aper + annulus, color='r', fill=False, ls='--', lw=.5)
        med_img = median_filter(image, (4, 4))[int(pltx[0]):round(int(pltx[1])), int(plty[0]):round(int(plty[1]))]
        norm = ImageNormalize(image, interval=ZScaleInterval(), stretch=stretch, vmin=np.nanpercentile(med_img, 5),
                              vmax=np.nanpercentile(med_img, 99))
        plt.imshow(image, norm=norm, origin='lower', cmap='Greys_r', interpolation=None)
        plt.plot(x_targ, y_targ, marker='+', color='lime')
        ax.add_artist(target_circle)
        ax.add_artist(target_circle_sky)
        if aper >= 0:
            ax.add_artist(ref_circle)
            ax.add_artist(ref_circle_sky)
            plt.plot(x_ref, y_ref, '+r')
        plt.xlabel("x-axis [pixel]")
        plt.ylabel("y-axis [pixel]")
        plt.title(f"FOV for {targ_name}\n({image_scale})")
        plt.xlim(pltx[0], pltx[1])
        plt.ylim(plty[0], plty[1])
        ax.grid(False)
        plt.plot(0, 0, color='lime', ls='-', label='Target')
        if aper >= 0:
            plt.plot(0, 0, color='r', ls='-.', label='Comp')
        l = plt.legend(framealpha=0.75)
        for text in l.get_texts():
            text.set_color("k")
            text.set_path_effects([path_effects.withStroke(linewidth=1, foreground='white')])

        apos = '\''
        Path(save).mkdir(parents=True, exist_ok=True)
        Path(save,"temp").mkdir(parents=True, exist_ok=True)

        plt.savefig(Path(save) / "temp" / f"FOV_{targ_name}_{date}_"
                    f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.pdf", bbox_inches='tight')
        plt.savefig(Path(save) / "temp" / f"FOV_{targ_name}_{date}_"
                    f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.png", bbox_inches='tight')
        plt.close()


def plot_flux(times, targ, targ_unc, ref, ref_unc, norm_flux, norm_unc, airmass, targ_name, save, date):
    plt.figure()
    plt.title(f"{targ_name} Raw Flux Values {date}")
    plt.xlabel("Time (BJD)")
    plt.ylabel("Total Flux")
    plt.errorbar(times, targ, yerr=targ_unc, linestyle='None', fmt='-o')
    plt.savefig(Path(save) / "temp" / f"TargetRawFlux_{targ_name}_{date}.png")
    plt.close()

    plt.figure()
    plt.title(f"Comparison Star Raw Flux Values {date}")
    plt.xlabel("Time (BJD)")
    plt.ylabel("Total Flux")
    plt.errorbar(times, ref, yerr=ref_unc, linestyle='None', fmt='-o')
    plt.savefig(Path(save) / "temp" / f"CompRawFlux_{targ_name}_{date}.png")
    plt.close()

    # Plots final reduced light curve (after the 3 sigma clip)
    plt.figure()
    plt.title(f"{targ_name} Normalized Flux vs. Time {date}")
    plt.xlabel("Time (BJD)")
    plt.ylabel("Normalized Flux")
    plt.errorbar(times, norm_flux, yerr=norm_unc, linestyle='None', fmt='-bo')
    plt.savefig(Path(save) / "temp" / f"NormalizedFluxTime_{targ_name}_{date}.png")
    plt.close()

    # Save normalized flux to text file prior to NS
    params_file = Path(save) / f"NormalizedFlux_{targ_name}_{date}.txt"
    with params_file.open('w') as f:
        f.write("BJD,Norm Flux,Norm Err,AM\n")

        for ti, fi, erri, ami in zip(times, norm_flux, norm_unc, airmass):
            f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")


# Observation statistics from PSF data
def plot_obs_stats(fit, comp_stars, psf, si, gi, target_name, save, date):
    for i in range(len(comp_stars) + 1):
        if i == 0:
            title, key = target_name, "target"
        else:
            title, key = f"Comp Star {i}", f"comp{i}"

        fig, axs = plt.subplots(3, 2, figsize=(12, 10))
        fig.suptitle(f"Observing Statistics - {title} - {date}")

        axs[0, 0].set(xlabel="Time [BJD_TBD]", ylabel="X-Centroid [px]")
        axs[0, 0].plot(fit.time, psf[key][si, 0][gi], 'k.')

        axs[0, 1].set(xlabel="Time [BJD_TBD]", ylabel="Y-Centroid [px]")
        axs[0, 1].plot(fit.time, psf[key][si, 1][gi], 'k.')

        axs[1, 0].set(xlabel="Time [BJD_TBD]", ylabel="Seeing [px]")
        axs[1, 0].plot(fit.time, 2.355 * 0.5 * (psf[key][si, 3][gi] + psf[key][si, 4][gi]), 'k.')

        axs[1, 1].set(xlabel="Time [BJD_TBD]", ylabel="Airmass")
        axs[1, 1].plot(fit.time, fit.airmass, 'k.')

        axs[2, 0].set(xlabel="Time [BJD_TBD]", ylabel="Amplitude [ADU]")
        axs[2, 0].plot(fit.time, psf[key][si, 2][gi], 'k.')

        axs[2, 1].set(xlabel="Time [BJD_TBD]", ylabel="Background [ADU]")
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
