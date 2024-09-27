import numpy as np
import matplotlib.pyplot as plt

# gaussian function
def gaussian(a, x, mu, sig):
    return a*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

if __name__ == "__main__":
    # compute the flux ratio between objects at different mags
    m1 = np.random.normal(11, 1)
    m2 = np.random.normal(13, 1)

    # compute the flux ratio
    f2f1 = 10**((m1 - m2) / 2.5)
    print(f"Flux ratio: {np.mean(f2f1):.2f}")
    print(f"m1 = {np.mean(m1):.2f}")
    print(f"m2 = {np.mean(m2):.2f}")

    distance = 55 # arc seconds
    seeing = 20 # arc seconds
    aperture = 75 # arc seconds

    # compute the probability of contamination
    x = np.linspace(-150, 150, 1000)
    psf1 = gaussian(1, x, 0, seeing)
    psf2 = gaussian(f2f1, x, distance, seeing)
    psf_total = psf1 + psf2

    # integrate psf1 between aperture
    amask = np.abs(x) < aperture
    flux1 = np.trapz(psf1[amask], x[amask])
    #print(f"Flux1: {flux1:.2f} within aperture")

    # integrate psf2 between aperture
    flux2 = np.trapz(psf2[amask], x[amask])
    #print(f"Flux2: {flux2:.2f} within aperture")

    # total flux:
    total_flux = flux1 + flux2
    #print(f"Total flux: {total_flux:.2f}")

    # fraction of flux from object 2
    fraction = flux2 / total_flux
    print(f"Fraction of flux from object 2: {fraction:.3f}")

    # contamination factor for modifying transit depth
    contamination = 1./(1 - fraction)
    print(f"Contamination factor: {contamination:.2f} (transit appears this much shallower)")

    # plot the distribution
    fig,ax = plt.subplots(1, figsize=(8,6))
    ax.plot(x, psf1, label='PSF1')
    ax.plot(x, psf2, label='PSF2')
    ax.plot(x, psf_total, label='Total PSF', ls='--')

    # draw vertical lines for aperture
    ax.axvline(-aperture, color='black', ls='--', lw=1, label='Aperture')
    ax.axvline(aperture, color='black', ls='--', lw=1)

    # overlap
    ax.fill_between(x[amask], psf2[amask], color='gray', alpha=0.5, label='Contamination')

    ax.set_xlabel('Distance (arc seconds)')
    ax.set_ylabel('Rel. Flux')
    ax.legend()
    plt.show()
