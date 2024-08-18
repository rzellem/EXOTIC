import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import (AsinhStretch, LogStretch, LinearStretch, ImageNormalize, ZScaleInterval, MinMaxInterval)
from matplotlib.widgets import Slider, RadioButtons
import matplotlib.image as mpimg

def load_and_autostretch_fits(filename):
    # open FITS file
    hdu = fits.open(filename)
    image_data = hdu[0].data
    return image_data

def update(val):
    # get selected stretch/interval
    selected_stretch = stretch_map[stretch_radio.value_selected]
    selected_interval = interval_map[interval_radio.value_selected]

    vmin, vmax = vmin_slider.val, vmax_slider.val

    # normalize image w/ selections
    norm = ImageNormalize(image_data, interval=selected_interval, stretch=selected_stretch, vmin=vmin, vmax=vmax)
    img.set_norm(norm)
    fig.canvas.draw_idle()

def on_hover(event, ax, image_data):
    if event.inaxes == ax:
        x, y = int(event.xdata), int(event.ydata)
        if 0 <= x < image_data.shape[1] and 0 <= y < image_data.shape[0]:
            ax.set_title(f"x={x}, y={y}")
            fig.canvas.draw_idle()

if __name__ == "__main__":
    # Load the FITS file
    # filename = "lum-bin-1-min-5-0001-20.fit"  # Replace with your FITS file
    filename = "lum-bin-1-min-5-0001-20.fit"  # Replace with your FITS file
    image_data = load_and_autostretch_fits(filename)

    # Define the available stretches and intervals
    stretch_map = {
        'Asinh': AsinhStretch(),
        'Log': LogStretch(),
        'Linear': LinearStretch()
    }

    interval_map = {
        'ZScale': ZScaleInterval(),
        'MinMax': MinMaxInterval()
    }

    # implemented this to be able to change a default setting for the image stretch on load
    # not sure how to normalize this for everyone's varying images. reminder to ask kyle/john
    # how they landed on a default with Bokeh. also thinking about a randomize button people can
    # just click through to populate a seemingly or actually random set of displays until they get
    # one that is close enough. sliders would need to update on button click to settings as well
    # so a user could hone in I guess. will get feedback. - Ira
    initial_stretch = 'Log'
    initial_interval = 'ZScale'
    initial_vmin = 219
    initial_vmax = 976

    norm = ImageNormalize(image_data, interval=interval_map[initial_interval], stretch=stretch_map[initial_stretch], vmin=initial_vmin, vmax=initial_vmax)
    
    # display image
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.subplots_adjust(left=0.25, bottom=0.35)
    img = ax.imshow(image_data, norm=norm, cmap='gray', origin='lower')

    # exoplanet watch logo
    logo = mpimg.imread('exoplanet-watch-logo.png')
    fig.figimage(logo, xo=10, yo=fig.bbox.ymax - logo.shape[0] - 10, zorder=10)

    # Asinh, Log, Linear radio buttons (need to decide if these are even needed, maybe just default to log or linear)
    ax_stretch = plt.axes([0.05, 0.5, 0.15, 0.15])
    stretch_radio = RadioButtons(ax_stretch, ('Asinh', 'Log', 'Linear'))
    stretch_radio.set_active(list(stretch_map.keys()).index(initial_stretch))

    # ZScale, MinMax radio buttons
    ax_interval = plt.axes([0.05, 0.25, 0.15, 0.15])
    interval_radio = RadioButtons(ax_interval, ('ZScale', 'MinMax'))
    interval_radio.set_active(list(interval_map.keys()).index(initial_interval))

    # sliders for vmin and vmax
    ax_vmin = plt.axes([0.25, 0.15, 0.65, 0.03])
    vmin_slider = Slider(ax_vmin, 'vmin', np.min(image_data), np.max(image_data), valinit=initial_vmin)
    
    ax_vmax = plt.axes([0.25, 0.1, 0.65, 0.03])
    vmax_slider = Slider(ax_vmax, 'vmax', np.min(image_data), np.max(image_data), valinit=initial_vmax)

    # this ties in the update function to the RadioButtons and Sliders
    stretch_radio.on_clicked(update)
    interval_radio.on_clicked(update)
    vmin_slider.on_changed(update)
    vmax_slider.on_changed(update)

    # connect hover event to the function
    fig.canvas.mpl_connect('motion_notify_event', lambda event: on_hover(event, ax, image_data))

    plt.show()
