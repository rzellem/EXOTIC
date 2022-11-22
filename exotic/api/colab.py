
# If the user presses enter to run the sample data, download sample data if needed and
# put it into a sample-data directory at the top level of the user's Gdrive.  Count
# the .fits files (images) and .json files (inits files) in the directory entered 
# by the user (or in the sample-data directory if the user pressed enter).  If 
# there are at least 20 .fits files, assume this is a directory of images and display
# the first one in the series.  If there is exactly one inits file in the directory, 
# show the specified target and comp coords so that the user can check these against
# the displayed image.  Otherwise, prompt for target / comp coords and make an inits 
# file based  on those (save this new inits file in the folder with the output files 
# so that the student can consult it later).  Finally, run EXOTIC with the newly-made 
# or pre-existing inits file, plus any other inits files in the directory.

#########################################################
from IPython.display import display, HTML
from astropy.time import Time
from barycorrpy import utc_tdb
import numpy as np
from io import BytesIO
from astropy.io import fits
from scipy.ndimage import label
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Viridis256
from bokeh.models import ColorBar, LinearColorMapper, LogColorMapper, LogTicker
from bokeh.models import BoxZoomTool,WheelZoomTool,ResetTool,HoverTool,PanTool,FreehandDrawTool
#import bokeh.io
#from bokeh.io import output_notebook
from pprint import pprint
#from IPython.display import Image
#from ipywidgets import widgets, HBox
from skimage.transform import rescale, resize, downscale_local_mean
#import copy
import os
import re
import json
#import subprocess
import time


def display_image(filename):
    #print(f"{filename}")
    hdu = fits.open(filename)

    extension = 0
    image_header = hdu[extension].header
    while image_header["NAXIS"] == 0:
      extension += 1
      image_header = hdu[extension].header

    dheader = dict(hdu[extension].header)
  
    data = hdu[extension].data
    megapixel_factor = (data.shape[0])*(data.shape[1])/1000000.0
    if megapixel_factor > 5:
      print(f"Downsampling image because it has {megapixel_factor} megapixels.")
      image_downscaled = downscale_local_mean(data, (2, 2)).astype(int)
      data = image_downscaled
    
    max_y = len(data)
    max_x = len(data[0])
    p_height = 500
    p_width = int((p_height/max_y) * max_x)

    # quick hot pixel/ cosmic ray mask
    # mask, cdata = detect_cosmics(
    #     data, psfmodel='gauss',
    #     psffwhm=4, psfsize=2*round(4)+1, # just a guess
    #     sepmed=False, sigclip = 4.25,
    #     niter=3, objlim=10, cleantype='idw', verbose=False
    # )

    # create a figure with text on mouse hover
    fig = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")], plot_width=p_width, plot_height=p_height,
        tools=[PanTool(),BoxZoomTool(),WheelZoomTool(),ResetTool(),HoverTool()])
    fig.x_range.range_padding = fig.y_range.range_padding = 0

    r = fig.multi_line('x', 'y', source={'x':[],'y':[]},color='white',line_width=3)
    fig.add_tools(FreehandDrawTool(renderers=[r]))

    # set up a colobar + data range
    color_mapper = LogColorMapper(palette="Cividis256", low=np.percentile(data, 55), high=np.percentile(data, 99))

    # must give a vector of image data for image parameter
    fig.image(
        image=[data],
          x=0, y=0, dw=hdu[extension].data.shape[1], dh=hdu[extension].data.shape[0],
          level="image", color_mapper=color_mapper
    )
    fig.grid.grid_line_width = 0.5

    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                         label_standoff=12, border_line_color=None, location=(0,0))

    fig.add_layout(color_bar, 'right')

    show(fig)

#########################################################

def floats_to_ints(l):
  while (True):
#    print (l)
    m = re.search(r"^(.*?)(\d+\.\d+)(.*?)$", l)
    if m:
      start, fl, end = m.group(1), float(m.group(2)), m.group(3)
      l = start+str("%.0f" % fl)+end
    else:
      return(l)
  
#########################################################

# Find a field in the image fits header or prompt the user to enter the corresponding
# value.

def check_dir(p):
  p = p.replace("\\", "/")

  if not(os.path.isdir(p)):
    print(f"Problem: the directory {p} doesn't seem to exist")
    print("on your Gdrive filesystem.")
    return("")
  return(p)

#########################################################

def add_sign(var):
  str_var = str(var)
  m=re.search(r"^[\+\-]", str_var)
  if m:
    return(str_var)
  if float(var) >= 0:
    return(str("+%.6f" % float(var)))
  else:
    return(str("-%.6f" % float(var)))

#########################################################

def get_val(hdr, ks):
  for key in ks:
    if key in hdr.keys():
      return hdr[key]
    if key.lower() in hdr.keys():
      return hdr[key.lower()]
    new_key = key[0]+key[1:len(key)].lower()  # first letter capitalized
    if new_key in hdr.keys():
      return hdr[new_key]
  return("")

#########################################################

def process_lat_long(val, key):
  m = re.search(r"\'?([+-]?\d+)[\s\:](\d+)[\s\:](\d+\.?\d*)", val)
  if m:
    deg, min, sec = float(m.group(1)), float(m.group(2)), float(m.group(3))
    if deg < 0:
      v = deg - (((60*min) + sec)/3600)
    else:
      v = deg + (((60*min) + sec)/3600)
    return(add_sign(v))
  m = re.search("^\'?([+-]?\d+\.\d+)", val)
  if m:
    v = float(m.group(1))
    return(add_sign(v))
  else:
    print(f"Cannot match value {val}, which is meant to be {key}.")

#########################################################

# Convert a MicroObservatory timestamp (which is in local time) to UTC.

def convert_Mobs_to_utc(datestamp, latitude, longitude, height):

#  print(datestamp)
  t = Time(datestamp[0:21], format='isot', scale='utc')
  t -= 0.33

  return(str(t)[0:10])

#########################################################

def find (hdr, ks, obs):
  # Special stuff for MObs and Boyce-Astro Observatories
  boyce = {"FILTER": "ip", "LATITUDE": "+32.6135", "LONGITUD": "-116.3334", "HEIGHT": 1405 }
  mobs = {"FILTER": "V", "LATITUDE": "+37.04", "LONGITUD": "-110.73", "HEIGHT": 2606 }

  if "OBSERVAT" in hdr.keys() and hdr["OBSERVAT"] == 'Whipple Observatory':
    obs = "MObs"

#  if "USERID" in hdr.keys() and hdr["USERID"] == 'PatBoyce':
#    obs = "Boyce"

  if obs == "Boyce":
    boyce_val = get_val(boyce, ks)
    if (boyce_val != ""):
      return(boyce_val)
  if obs == "MObs":
    mobs_val = get_val(mobs, ks)
    if (mobs_val != ""):
      return(mobs_val)

  val = get_val(hdr, ks)

  if ks[0] == "LATITUDE" and val != "":         # Because EXOTIC needs these with signs shown.
    return(process_lat_long(str(val), "latitude"))
  if ks[0] == "LONGITUD" and val != "":
    return(process_lat_long(str(val), "longitude"))

  if (val != ""):
    return(val)

  print(f"\nI cannot find a field with any of these names in your image header: \n{ks}.")
  print("Please enter the value (not the name of the header field, the actual value) that should")
  print("be used for the value associated with this field.\n")
  if ks[0] == "HEIGHT":
    print("The units of elevation are meters.")
  
  value = input("")

  return(value)

###############################################

def look_for_calibration(image_dir):
  darks_dir, flats_dir, biases_dir = "null", "null", "null"

  m = re.search(r"(.*?)(\d\d\d\d\-\d\d\-\d\d)\/images", image_dir)  # This handles the way I set up the MObs image paths for my seminar teams.
  if m:
    prefix, date = m.group(1), m.group(2)
    darks = prefix+date+"/darks"
    if os.path.isdir(darks):
      darks_dir = str("\""+darks+"\"")
      
  d_names = ["dark", "darks", "DARK", "DARKS", "Dark", "Darks"]  # Possible names for calibration image directories.
  f_names = ["flat", "flats", "FLAT", "FLATS", "Flat", "Flats"]
  b_names = ["bias", "biases", "BIAS", "BIASES", "Bias", "Biases"]

  for d in d_names:
    if os.path.isdir(os.path.join(image_dir, d)):
      darks_dir = str("\""+os.path.join(image_dir, d)+"\"")
      break

  for f in f_names:
    if os.path.isdir(os.path.join(image_dir, f)):
      flats_dir = str("\""+os.path.join(image_dir, f)+"\"")
      break

  for b in b_names:
    if os.path.isdir(os.path.join(image_dir, b)):
      biases_dir = str("\""+os.path.join(image_dir, b)+"\"")
      break

  return(darks_dir, flats_dir, biases_dir)

###############################################

# Writes a new inits file into the directory with the output plots.  This prompts
# for needed information that it cannot find in the fits header of the first image.

def make_inits_file(planetary_params, image_dir, output_dir, first_image, targ_coords, comp_coords, obs, aavso_obs_code, sec_obs_code, sample_data):
  inits_file_path = output_dir+"inits.json"
  hdul = fits.open(first_image)
  hdr = dict(hdul[0].header)

  min, max = "null", "null"
  filter = find(hdr, ['FILTER', 'FILT'], obs)
  if filter == "w":
    filter = "PanSTARRS-w"
    min = "404"
    max = "846"
  if filter == "Clear":
    filter = "V"
  if filter == "ip":
    min = "690"
    max = "819"
  if filter == "EXO":
    filter = "CBB"
  if re.search(r"Green", filter, re.IGNORECASE):
    filter = "SG"
    
  date_obs = find(hdr,["DATE", "DATE_OBS", "DATE-OBS"], obs)
  date_obs = date_obs.replace("/", "_")
  longitude = find(hdr,['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG'],obs)
  latitude = find(hdr,['LATITUDE', 'LAT', 'SITELAT'],obs)
  height = int(find(hdr, ['HEIGHT', 'ELEVATION', 'ELE', 'EL', 'OBSGEO-H', 'ALT-OBS', 'SITEELEV'], obs))
  obs_notes = "N/A"

  mobs_data = False
  # For MObs, the date is local rather than UTC, so convert.
  if "OBSERVAT" in hdr.keys() and hdr["OBSERVAT"] == 'Whipple Observatory':
    date_obs = convert_Mobs_to_utc(date_obs, latitude, longitude, height)
    weather = hdr["WEATHER"] 
    temps = float(hdr["TELTEMP"]) - float(hdr["CAMTEMP"])
    obs_notes = str("First image seeing %s (0: poor, 100: excellent), Teltemp - Camtemp %.1f.  These observations were conducted with MicroObservatory, a robotic telescope network managed by the Harvard-Smithsonian Center for Astrophysics on behalf of NASA's Universe of Learning. This work is supported by NASA under award number NNX16AC65A to the Space Telescope Science Institute." % (weather, temps))
    sec_obs_code = "MOBS"  
    mobs_data = True
  
  if aavso_obs_code == "":
      aavso_obs_code = "N/A"
  if sec_obs_code == "":
      sec_obs_code = "N/A"

  obs_date = date_obs[0:10]
  (darks_dir, flats_dir, biases_dir) = look_for_calibration(image_dir)

  with open(inits_file_path, 'w') as inits_file:
    inits_file.write("""
{
  %s,
    "user_info": {
            "Directory with FITS files": "%s",
            "Directory to Save Plots": "%s",
            "Directory of Flats": %s,
            "Directory of Darks": %s,
            "Directory of Biases": %s,

            "AAVSO Observer Code (N/A if none)": "%s",
            "Secondary Observer Codes (N/A if none)": "%s",

            "Observation date": "%s",
            "Obs. Latitude": "%s",
            "Obs. Longitude": "%s",
            "Obs. Elevation (meters)": %d,
            "Camera Type (CCD or DSLR)": "CCD",
            "Pixel Binning": "1x1",
            "Filter Name (aavso.org/filters)": "%s",
            "Observing Notes": "%s",

            "Plate Solution? (y/n)": "y",
            "Align Images? (y/n)": "y",

            "Target Star X & Y Pixel": %s,
            "Comparison Star(s) X & Y Pixel": %s
    },    
    "optional_info": {
            "Pixel Scale (Ex: 5.21 arcsecs/pixel)": null,
            "Filter Minimum Wavelength (nm)": %s,
            "Filter Maximum Wavelength (nm)": %s
    }
}
""" % (planetary_params, image_dir, output_dir, flats_dir, darks_dir, biases_dir, 
       aavso_obs_code, sec_obs_code, obs_date, latitude, longitude, height, filter, 
       obs_notes, targ_coords, comp_coords, min, max))

  display(HTML('<p class="output"><b>Initialization File Created.</b></p>'))
  print(f'Created: {inits_file_path}')
  print('This folder will also contain the output files when EXOTIC finishes running.')

  if not mobs_data:  
    print(f"\nThe inits.json file currently says that your observatory latitude was {latitude} deg,")
    print(f"longitude was {longitude} deg, and elevation was {height}m.  \n")
    print("*** If any of these are incorrect, please change them in the inits.json file. ***")
    print("*** (Please make sure that Western longitudes have a negative sign! ***")
    print("*** TheSkyX sometimes stamps Western longitudes as positive; this needs to be switched! ***\n")

  display(HTML('<p class="output"><br /><b>If you want to change anything in the inits file, such as planetary parameters or user info, please do that now.</b></p>'))
  print('You can edit the file by clicking the folder icon in the left nav,')
  print(f'navigating to the inits file at {inits_file_path}, and double-clicking the file.')
  print('\nWhen you are done, save your changes, and proceed to the next step.')
  

  return(inits_file_path)
  
##############################################################

def fix_planetary_params (p_param_dict):
  for param in p_param_dict.keys():
    if param == "Target Star RA" or param == "Target Star Dec" or param == "Planet Name" or param == "Host Star Name" or param == "Argument of Periastron (deg)":
      continue
    val = p_param_dict[param]
    if val == 0.0 or np.isnan(float(val)):
      if param == "Orbital Eccentricity (0 if null)":
        continue
      if param == "Ratio of Planet to Stellar Radius (Rp/Rs)":
        p_param_dict[param] = 0.151
      if param == "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty":
        p_param_dict[param] = 0.151
        if p_param_dict["Host Star Name"] == "Qatar-6":
          p_param_dict[param] = 0.01
      print(f"\nIn the planetary parameters from the NASA Exoplanet Archive, \n\"{param}\" is listed as {val}.\n\n**** This might make EXOTIC crash. ****\n\nIf the parameter is *not* changed below, please edit it\nin the inits file before running EXOTIC.\n")
  p_param_string = json.dumps(p_param_dict)

  planetary_params = "\"planetary_parameters\": {\n"
  num_done, num_total = 0, len(p_param_dict.keys())
  for key, value in p_param_dict.items():
    num_done += 1
    if key == "Target Star RA" or key == "Target Star Dec" or key == "Planet Name" or key == "Host Star Name":
      planetary_params = planetary_params + str(f"    \"{key}\": \"{value}\",\n")
    else:
      if num_done < num_total:
        planetary_params = planetary_params + str(f"    \"{key}\": {value},\n")
      else:
        planetary_params = planetary_params + str(f"    \"{key}\": {value}\n")
  planetary_params = planetary_params + "}"

  return(planetary_params)
