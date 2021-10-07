import json
import numpy as np
from io import BytesIO
from astropy.io import fits
# from astroscrappy import detect_cosmics
from scipy.ndimage import label
from skimage.transform import downscale_local_mean
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Viridis256
from bokeh.models import ColorBar, LinearColorMapper, LogColorMapper, LogTicker
from bokeh.models import BoxZoomTool,WheelZoomTool,ResetTool,HoverTool,PanTool,FreehandDrawTool
from bokeh.io import output_notebook
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

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
    else:
        image_downscaled = downscale_local_mean(data, (1, 1)).astype(int)

    # quick hot pixel/ cosmic ray mask
    # mask, cdata = detect_cosmics(
    #     data, psfmodel='gauss',
    #     psffwhm=4, psfsize=2*round(4)+1, # just a guess
    #     sepmed=False, sigclip = 4.25,
    #     niter=3, objlim=10, cleantype='idw', verbose=False
    # )
    mask, cdata = None, None  # temp til astroscrappy works again for python v3.9

    # show how many pixels are saturated

    SATURATION = 2**(hdu[extension].header['bitpix'])
    if SATURATION == None:
        SATURATION = np.max(data)
    mmask = image_downscaled >= SATURATION*0.9
    labels, ngroups = label(mmask)
    print('Saturated Areas:',ngroups)
    print('Saturation level:',SATURATION)
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


# rip off of corner.py so we can cmap scatter to chi2
def corner(xs, bins=20, range=None, weights=None, color="k", hist_bin_factor=1,
           smooth=None, smooth1d=None, levels=[1],
           labels=None, label_kwargs=None,
           titles=[], title_fmt=".2f", title_kwargs=None,
           truths=None, truth_color="#4682b4",
           scale_hist=False, quantiles=None, verbose=False, fig=None,
           max_n_ticks=5, top_ticks=False, use_math_text=False, reverse=False,
           hist_kwargs=None, **hist2d_kwargs):

    if quantiles is None:
        quantiles = []
    if title_kwargs is None:
        title_kwargs = dict()
    if label_kwargs is None:
        label_kwargs = dict()

    # Try filling in labels from pandas.DataFrame columns.
    if labels is None:
        try:
            labels = xs.columns
        except AttributeError:
            pass

    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
                                       "dimensions than samples!"

    # Parse the weight array.
    if weights is not None:
        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("Weights must be 1-D")
        if xs.shape[1] != weights.shape[0]:
            raise ValueError("Lengths of weights must match number of samples")

    # Parse the parameter ranges.
    if range is None:
        if "extents" in hist2d_kwargs:
            logging.warn("Deprecated keyword argument 'extents'. "
                         "Use 'range' instead.")
            range = hist2d_kwargs.pop("extents")
        else:
            range = [[x.min(), x.max()] for x in xs]
            # Check for parameters that never change.
            m = np.array([e[0] == e[1] for e in range], dtype=bool)
            if np.any(m):
                raise ValueError(("It looks like the parameter(s) in "
                                  "column(s) {0} have no dynamic range. "
                                  "Please provide a `range` argument.")
                                 .format(", ".join(map(
                                     "{0}".format, np.arange(len(m))[m]))))

    else:
        # If any of the extents are percentiles, convert them to ranges.
        # Also make sure it's a normal list.
        range = list(range)
        for i, _ in enumerate(range):
            try:
                emin, emax = range[i]
            except TypeError:
                q = [0.5 - 0.5*range[i], 0.5 + 0.5*range[i]]
                range[i] = quantile(xs[i], q, weights=weights)

    if len(range) != xs.shape[0]:
        raise ValueError("Dimension mismatch between samples and range")

    # Parse the bin specifications.
    try:
        bins = [int(bins) for _ in range]
    except TypeError:
        if len(bins) != len(range):
            raise ValueError("Dimension mismatch between bins and range")
    try:
        hist_bin_factor = [float(hist_bin_factor) for _ in range]
    except TypeError:
        if len(hist_bin_factor) != len(range):
            raise ValueError("Dimension mismatch between hist_bin_factor and "
                             "range")

    # Some magic numbers for pretty axis layout.
    K = len(xs)
    factor = 2.0           # size of one side of one panel
    if reverse:
        lbdim = 0.2 * factor   # size of left/bottom margin
        trdim = 0.5 * factor   # size of top/right margin
    else:
        lbdim = 0.5 * factor   # size of left/bottom margin
        trdim = 0.2 * factor   # size of top/right margin
    whspace = 0.05         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    # Create a new figure if one wasn't provided.
    if fig is None:
        fig, axes = plt.subplots(K, K, figsize=(dim, dim))
    else:
        try:
            axes = np.array(fig.axes).reshape((K, K))
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                             "dimensions K={1}".format(len(fig.axes), K))

    # Format the figure.
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    # Set up the default histogram keywords.
    if hist_kwargs is None:
        hist_kwargs = dict()
    hist_kwargs["color"] = hist_kwargs.get("color", color)
    if smooth1d is None:
        hist_kwargs["histtype"] = hist_kwargs.get("histtype", "step")

    for i, x in enumerate(xs):
        # Deal with masked arrays.
        if hasattr(x, "compressed"):
            x = x.compressed()

        if np.shape(xs)[0] == 1:
            ax = axes
        else:
            if reverse:
                ax = axes[K-i-1, K-i-1]
            else:
                ax = axes[i, i]
        # Plot the histograms.
        if smooth1d is None:
            bins_1d = int(max(1, np.round(hist_bin_factor[i] * bins[i])))
            n, _, _ = ax.hist(x, bins=bins_1d, weights=weights,
                              range=np.sort(range[i]), **hist_kwargs)
        else:
            if gaussian_filter is None:
                raise ImportError("Please install scipy for smoothing")
            n, b = np.histogram(x, bins=bins[i], weights=weights,
                                range=np.sort(range[i]))
            n = gaussian_filter(n, smooth1d)
            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            ax.plot(x0, y0, **hist_kwargs)

        if truths is not None and truths[i] is not None:
            ax.axvline(truths[i], color=truth_color)

        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            qvalues = quantile(x, quantiles, weights=weights)
            for q in qvalues:
                ax.axvline(q, ls="dashed", color=color)

            if verbose:
                print("Quantiles:")
                print([item for item in zip(quantiles, qvalues)])

        if len(titles):
            title = None
            ax.set_title(titles[i], **title_kwargs)

        # Set up the axes.
        ax.set_xlim(range[i])
        if scale_hist:
            maxn = np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        if max_n_ticks == 0:
            ax.xaxis.set_major_locator(NullLocator())
            ax.yaxis.set_major_locator(NullLocator())
        else:
            ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
            ax.yaxis.set_major_locator(NullLocator())

        if i < K - 1:
            if top_ticks:
                ax.xaxis.set_ticks_position("top")
                [l.set_rotation(45) for l in ax.get_xticklabels()]
            else:
                ax.set_xticklabels([])
        else:
            if reverse:
                ax.xaxis.tick_top()
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                if reverse:
                    ax.set_title(labels[i], y=1.25, **label_kwargs)
                else:
                    ax.set_xlabel(labels[i], **label_kwargs)

            # use MathText for axes ticks
            ax.xaxis.set_major_formatter(
                ScalarFormatter(useMathText=use_math_text))

        for j, y in enumerate(xs):
            if np.shape(xs)[0] == 1:
                ax = axes
            else:
                if reverse:
                    ax = axes[K-i-1, K-j-1]
                else:
                    ax = axes[i, j]
            if j > i:
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            elif j == i:
                continue

            # Deal with masked arrays.
            if hasattr(y, "compressed"):
                y = y.compressed()

            hist2d(y, x, ax=ax, range=[range[j], range[i]], weights=weights,
                    smooth=smooth, bins=[bins[j], bins[i]], levels=levels,
                    **hist2d_kwargs)

            if truths is not None:
                if truths[i] is not None and truths[j] is not None:
                    ax.plot(truths[j], truths[i], "s", color=truth_color)
                if truths[j] is not None:
                    ax.axvline(truths[j], color=truth_color)
                if truths[i] is not None:
                    ax.axhline(truths[i], color=truth_color)

            if max_n_ticks == 0:
                ax.xaxis.set_major_locator(NullLocator())
                ax.yaxis.set_major_locator(NullLocator())
            else:
                ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                       prune="lower"))
                ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                       prune="lower"))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                if reverse:
                    ax.xaxis.tick_top()
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j], **label_kwargs)
                    if reverse:
                        ax.xaxis.set_label_coords(0.5, 1.4)
                    else:
                        ax.xaxis.set_label_coords(0.5, -0.3)

                # use MathText for axes ticks
                ax.xaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text))

            if j > 0:
                ax.set_yticklabels([])
            else:
                if reverse:
                    ax.yaxis.tick_right()
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    if reverse:
                        ax.set_ylabel(labels[i], rotation=-90, **label_kwargs)
                        ax.yaxis.set_label_coords(1.3, 0.5)
                    else:
                        ax.set_ylabel(labels[i], **label_kwargs)
                        ax.yaxis.set_label_coords(-0.3, 0.5)

                # use MathText for axes ticks
                ax.yaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text))

    return fig

def quantile(x, q, weights=None):
    """
    Compute sample quantiles with support for weighted samples.

    Note
    ----
    When ``weights`` is ``None``, this method simply calls numpy's percentile
    function with the values of ``q`` multiplied by 100.

    Parameters
    ----------
    x : array_like[nsamples,]
       The samples.

    q : array_like[nquantiles,]
       The list of quantiles to compute. These should all be in the range
       ``[0, 1]``.

    weights : Optional[array_like[nsamples,]]
        An optional weight corresponding to each sample. These

    Returns
    -------
    quantiles : array_like[nquantiles,]
        The sample quantiles computed at ``q``.

    Raises
    ------
    ValueError
        For invalid quantiles; ``q`` not in ``[0, 1]`` or dimension mismatch
        between ``x`` and ``weights``.

    """
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)

    if np.any(q < 0.0) or np.any(q > 1.0):
        raise ValueError("Quantiles must be between 0 and 1")

    if weights is None:
        return np.percentile(x, list(100.0 * q))
    else:
        weights = np.atleast_1d(weights)
        if len(x) != len(weights):
            raise ValueError("Dimension mismatch: len(weights) != len(x)")
        idx = np.argsort(x)
        sw = weights[idx]
        cdf = np.cumsum(sw)[:-1]
        cdf /= cdf[-1]
        cdf = np.append(0, cdf)
        return np.interp(q, cdf, x[idx]).tolist()

def hist2d(x, y, bins=20, range=None, levels=[2],
           ax=None, plot_datapoints=True, plot_contours=True, 
           contour_kwargs=None, contourf_kwargs=None, data_kwargs=None,
            **kwargs):
    if ax is None:
        ax = plt.gca()

    if plot_datapoints:
        if data_kwargs is None:
            data_kwargs = dict()
        data_kwargs["s"] = data_kwargs.get("s", 2.0)
        data_kwargs["alpha"] = data_kwargs.get("alpha", 0.2)
        ax.scatter(x, y, marker="o", zorder=-1, rasterized=True, **data_kwargs)

    # Plot the contour edge colors.
    if plot_contours:
        if contour_kwargs is None:
            contour_kwargs = dict()

        # mask data in range + chi2
        maskx = (x > range[0][0]) & (x < range[0][1])
        masky = (y > range[1][0]) & (y < range[1][1])
        mask = maskx & masky & (data_kwargs['c'] < data_kwargs['vmax']*1.2)
        
        try: # contour
            # approx posterior + smooth
            xg, yg = np.meshgrid( np.linspace(x[mask].min(),x[mask].max(),256), np.linspace(y[mask].min(),y[mask].max(),256) )
            cg = griddata(np.vstack([x[mask],y[mask]]).T, data_kwargs['c'][mask], (xg,yg), method='nearest', rescale=True)
            scg = gaussian_filter(cg,sigma=15)

            ax.contour(xg, yg, scg*np.nanmin(cg)/np.nanmin(scg), np.sort(levels), **contour_kwargs, vmin=data_kwargs['vmin'], vmax=data_kwargs['vmax'])        
        except Exception as err:
            print(err)
            print("contour plotting failed")
    
    ax.set_xlim(range[0])
    ax.set_ylim(range[1])
