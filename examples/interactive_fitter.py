''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve interactive_fitter.py
at your command prompt. Then navigate to the URL
http://localhost:5006/interactive_fitter in your browser.'''

import json
import numpy as np
import urllib.request
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure

from exotic.api.elca import transit
from exotic.exotic import NASAExoplanetArchive, get_wcs, find_target

# physical constants
rsun = 6.955e8 #m
msun = 1.989e30 #kg
au = 1.496e11 #m
G = 0.00029591220828559104 # day, AU, Msun

# load data from file
data = json.load(open("sampledata.json", "r"))
data['time'] = np.array(data['time'])

# download priors
target = NASAExoplanetArchive(data["name"])

if target.resolve_name():
    planet_name, _, parameters = target.planet_info()
else:
    raise(Exception(f"can not resolve planet name: {data['name']} in NEA"))

prior = {
    'rprs':parameters['rprs'],     # Rp/Rs
    'ars':parameters['aRs'],       # a/Rs
    'per':parameters['pPer'],      # Period [day]
    'inc':parameters['inc'],       # Inclination [deg]
    'u0': 0, 'u1': 0, 'u2': 0, 'u3': 0,  # limb darkening (nonlinear)
    'ecc': parameters['ecc'],            # Eccentricity
    'omega': 0,                     # Arg of periastron
    'tmid':parameters['midT'],      # time of mid transit [day]
}

print(json.dumps(parameters,indent=4))


# generate transit model
tdata = transit(data['time'], prior)
source_transit = ColumnDataSource(data=dict(x=data['time'], y=tdata))

# format data for bokeh
source_noisy = ColumnDataSource(data=dict(x=data['time'], y=data['flux']))

# Light curve plot
plot = figure(plot_height=400, plot_width=400, title="Transit Light Curve",
              tools="crosshair,pan,reset,save,wheel_zoom",
              y_range=[0.95, 1.01])

plot.line('x', 'y', source=source_transit, line_width=3, line_alpha=0.6)
plot.circle('x', 'y', source=source_noisy, color='black')
plot.xaxis.axis_label = "Time [day]"
plot.yaxis.axis_label = "Relative Flux"


# Sun + Planet plot
plot2 = figure(plot_height=400, plot_width=400, x_range=[-2,2], y_range=[-2, 2], title="Exoplanet Transit Simulator")

#sun
sundata = ColumnDataSource(data={'x':[0], 'y':[0]})
plot2.circle('x', 'y', source=sundata, fill_color='orange', radius=1, line_color='yellow')

#planet
planetdata = ColumnDataSource(data={'x':[0], 'y':[0], 's':[prior['rprs']]})
plot2.circle('x', 'y', source=planetdata, fill_color='black', radius='s', line_color='green')

# Set up widgets
rprs2 = Slider(title="Transit Depth (Rp/Rs)^2", value=prior['rprs']**2, start=0.00, end=0.04, step=0.001)
inc = Slider(title="Inclination (i)", value=prior['inc'], start=80, end=90, step=0.1)
tmid = Slider(title="Mid-Transit (Tmid)", value=prior['tmid'], start=min(data['time']), end=max(data['time']), step=0.001)

# Example API call from DIY Planet Search: https://waps.cfa.harvard.edu/microobservatory/diy/api/v1/datasets/dataset.php?target=tres-3-b
# TODO add text input for a url             (https://docs.bokeh.org/en/latest/docs/user_guide/interaction/widgets.html#textinput)
# TODO add button called "Load from URL"    (https://docs.bokeh.org/en/latest/docs/user_guide/interaction/widgets.html#button)
# TODO create a function that loads a json from the text box url and updates the data plot
#  when updating the plot data, x-axis use key: "timeJD", y-axis use key: "brightness"
# webtext = urllib.request.urlopen(url).read()
# newdata = json.loads(webtext)

def update_data(attrname, old, new):
    prior['rprs'] = rprs2.value**0.5
    prior['inc'] = inc.value
    prior['tmid'] = tmid.value
    offset = prior['ars']*np.cos(np.deg2rad(prior['inc'])) #hypotenuse #b

    # update plot data
    source_transit.data = dict(x=data['time'], y=transit(data['time'], prior))
    planetdata.data = dict(x=[0], y=[offset], s=[prior['rprs']])

for w in [rprs2, inc, tmid]:
    w.on_change('value', update_data)

# TODO set button.on_change to function above

# Set up layouts and add to document
inputs = column(rprs2, tmid, inc) # TODO add text box and button to column
plots = row(plot, plot2)
curdoc().add_root(column(plots, inputs, width=800))
curdoc().title = "Exoplanet Transit"
