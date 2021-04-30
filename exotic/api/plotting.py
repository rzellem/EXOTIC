import json
import numpy as np
from io import BytesIO
from astropy.io import fits
from astroscrappy import detect_cosmics
from scipy.ndimage import label
from skimage.transform import downscale_local_mean
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Viridis256
from bokeh.models import ColorBar, LinearColorMapper, LogColorMapper, LogTicker
from bokeh.models import BoxZoomTool,WheelZoomTool,ResetTool,HoverTool,PanTool,FreehandDrawTool
from bokeh.io import output_notebook
from pprint import pprint

def plot_image(filename, save=False, bg_min=60, bg_max=99):

    hdu = fits.open(filename)

    extension = 0
    image_header = hdu[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdu[extension].header

    dheader = dict(hdu[extension].header)
    djson = {'filename':filename}
    for k in dheader:
        if len(k) >= 2:
            print(f"{k}: {dheader[k]}")
        djson[k] = str(dheader[k])

    data = hdu[extension].data

    with open('header.json', 'w') as json_file: 
        json.dump(djson, json_file, indent=4)
        print("Image header written to header.json")

    if data.shape[0] > 6000:
        image_downscaled = downscale_local_mean(data, (4, 4)).astype(int)
    elif data.shape[0] > 2000:
        image_downscaled = downscale_local_mean(data, (2, 2)).astype(int)

    # quick hot pixel/ cosmic ray mask
    mask, cdata = detect_cosmics(
        data, psfmodel='gauss',
        psffwhm=4, psfsize=2*round(4)+1, # just a guess
        sepmed=False, sigclip = 4.25,
        niter=3, objlim=10, cleantype='idw', verbose=False
    )

    # show how many pixels are saturated
    SATURATION = 2**(hdu[extension].header['bitpix'])
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

    pprint(bad_pix)

    # create a figure with text on mouse hover\
    print("Saturated pixels are marked with red. These are pixels which have exceeded the maximum value for brightness, and are thus not suitable for use as comparison stars.")
    fig = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")], plot_width=800, plot_height=800,
        tools=[PanTool(),BoxZoomTool(),WheelZoomTool(),ResetTool(),HoverTool()])
    fig.x_range.range_padding = fig.y_range.range_padding = 0

    r = fig.multi_line('x', 'y', source={'x':[],'y':[]},color='white',line_width=3)
    fig.add_tools(FreehandDrawTool(renderers=[r]))

    # set up a colobar + data range
    color_mapper = LogColorMapper(palette="Cividis256", low=np.percentile(data, bg_min), high=np.percentile(data, bg_max))

    # must give a vector of image data for image parameter
    fig.image(
        image=[image_downscaled],
          x=0, y=0, dw=hdu[extension].data.shape[1], dh=hdu[extension].data.shape[0],
          level="image", color_mapper=color_mapper
    )

    # plot saturated stars
    fig.x(bad_pix['x'], bad_pix['y'], size=25, color='red', line_width=3)
    fig.x(bad_pix['x'], bad_pix['y'], size=25, color='white', line_width=1)
    # TODO figure out hover value

    fig.grid.grid_line_width = 0.5

    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                         label_standoff=12, border_line_color=None, location=(0,0))

    fig.add_layout(color_bar, 'right')

    if save:
        output_file("interactivefits.html")
    else:
        show(fig)