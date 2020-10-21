import numpy as np
from io import BytesIO
from astropy.io import fits
from astroscrappy import detect_cosmics
from scipy.ndimage import label
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Viridis256
from bokeh.models import ColorBar, LinearColorMapper
from bokeh.io import output_notebook

def plot_image(filename):

    hdu = fits.open(filename)
    print(hdu.info())
    print(hdu[0].header)
    data = hdu[0].data

    # quick hot pixel/ cosmic ray mask
    mask, cdata = detect_cosmics(
        data, psfmodel='gauss',
        psffwhm=4, psfsize=2*round(4)+1, # just a guess
        sepmed=False, sigclip = 4.25,
        niter=3, objlim=10, cleantype='idw', verbose=False
        )


    # show how many pixels are saturated
    SATURATION = 65536
    mmask = cdata >= SATURATION*0.9
    labels, ngroups = label(mmask)
    print('Saturated Areas:',ngroups)
    labeli, counts = np.unique(labels, return_counts=True)
    bad_pix = {'x':[], 'y':[], 'value':[]}
    # loop through each group to find position
    for i in range(1,labeli[-1]+1):
        imask = labels == i
        yc,xc = np.argwhere(imask).mean(0)
        bad_pix['x'].append(xc)
        bad_pix['y'].append(yc)
        bad_pix['value'].append(cdata[imask].mean())

    print(bad_pix)


    # create a figure with text on mouse hover
    fig = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")], plot_width=800, plot_height=800)
    fig.x_range.range_padding = fig.y_range.range_padding = 0

    # set up a colobar + data range
    mapper = LinearColorMapper( palette=Viridis256, low=np.percentile(cdata,55), high=np.percentile(cdata,99.5))

    # must give a vector of image data for image parameter
    fig.image(
        image=[cdata],
          x=0, y=0, dw=hdu[0].data.shape[1], dh=hdu[0].data.shape[0],
          level="image", color_mapper=mapper
    )

    # plot saturated stars
    fig.x(bad_pix['x'], bad_pix['y'], size=25, color='red')
    # TODO figure out hover value

    fig.grid.grid_line_width = 0.5
    # TODO show a color bar

    #output_file("image.html", title="fts example")
    show(fig)



##Testing
plot_image("sample-data/HatP32Dec202017/HATP-32171220013343.fits")
