''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve sliders.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders
in your browser.
'''
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure

from exotic.api.elca import transit

prior = {
    'rprs':0.1,        # Rp/Rs
    'ars':14.25,        # a/Rs
    'per':3.336817,     # Period [day]
    'inc':87.5,        # Inclination [deg]
    'u0': 0, 'u1': 0, 'u2': 0, 'u3': 0,  # limb darkening (nonlinear)
    'ecc':0,            # Eccentricity
    'omega':0,          # Arg of periastron
    'tmid':0.5,         # time of mid transit [day]
}

time = np.linspace(0.4,0.6,200) # [day]
data = transit(time, prior)
source = ColumnDataSource(data=dict(x=time, y=data))

data += np.random.normal(0, 500e-6, len(time))
source_noisy = ColumnDataSource(data=dict(x=time, y=data))

# Set up plot
plot = figure(plot_height=400, plot_width=400, title="Interactive Transit",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0.4,0.6], y_range=[0.98, 1.01])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
plot.circle('x', 'y', source=source_noisy, color='black')
plot.xaxis.axis_label = "Time [day]"
plot.yaxis.axis_label = "Relative Flux"

# Set up widgets
#text = TextInput(title="title", value='my transit')

rprs = Slider(title="Transit Depth (Rp/Rs)", value=prior['rprs'], start=0.01, end=0.15, step=0.001)
inc = Slider(title="Inclination (i)", value=90, start=85, end=90, step=0.1)
tmid = Slider(title="Mid-Transit (Tmid)", value=0.5, start=0.45, end=0.55, step=0.001)
per = Slider(title="Period", value=prior['per'], start=prior['per']*0.8, end=prior['per']*1.2, step=prior['per']*0.01)

# Set up callbacks
#def update_title(attrname, old, new):
#    plot.title.text = text.value

#text.on_change('value', update_title)

def update_data(attrname, old, new):
    prior['rprs'] = rprs.value
    prior['inc'] = inc.value
    prior['tmid'] = tmid.value
    prior['per'] = per.value
    source.data = dict(x=time, y=transit(time, prior))

for w in [rprs, inc, tmid, per]:
    w.on_change('value', update_data)

# Set up layouts and add to document
inputs = column(rprs, tmid, per, inc)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"