import importlib
import json
import sys
import types

import numpy as np
from astropy.io import fits


def import_colab(monkeypatch):
    fake_barycorrpy = types.ModuleType("barycorrpy")
    fake_barycorrpy.utc_tdb = object()
    fake_ipython = types.ModuleType("IPython")
    fake_display = types.ModuleType("IPython.display")
    fake_display.display = lambda *args, **kwargs: None
    fake_display.HTML = lambda value: value

    monkeypatch.setitem(sys.modules, "barycorrpy", fake_barycorrpy)
    monkeypatch.setitem(sys.modules, "IPython", fake_ipython)
    monkeypatch.setitem(sys.modules, "IPython.display", fake_display)

    return importlib.import_module("exotic.api.colab")


def test_mobs_defaults_match_exotic_metadata(monkeypatch, tmp_path):
    colab = import_colab(monkeypatch)
    fits_path = tmp_path / "mobs.fits"
    fits.PrimaryHDU(data=np.zeros((2, 2))).writeto(fits_path)

    with fits.open(fits_path, mode="update") as hdul:
        hdr = hdul[0].header
        hdr["OBSERVAT"] = "Whipple Observatory"
        hdr["DATE"] = "2017-12-20T01:33:43.000"
        hdr["WEATHER"] = 80
        hdr["TELTEMP"] = 20.0
        hdr["CAMTEMP"] = 10.0

    hdr = fits.getheader(fits_path)
    assert colab.find(hdr, ["FILTER", "FILT"], "MObs") == "CV"
    assert colab.find(hdr, ["LATITUDE", "LAT", "SITELAT"], "MObs") == "+31.675467"
    assert colab.find(hdr, ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"], "MObs") == "-110.951376"
    assert colab.find(hdr, ["HEIGHT", "ELEVATION", "ELE", "EL", "OBSGEO-H", "ALT-OBS", "SITEELEV"], "MObs") == 1268

    image_dir = tmp_path.as_posix()
    output_dir = f"{tmp_path.as_posix()}/"
    inits_path = colab.make_inits_file(
        '"planetary_parameters": {}',
        image_dir,
        output_dir,
        fits_path.as_posix(),
        "[424, 286]",
        "[[465, 183], [512, 263]]",
        "MObs",
        "RTZ",
        "",
        False,
    )
    inits = json_loads(inits_path)
    user_info = inits["user_info"]

    assert user_info["Filter Name (aavso.org/filters)"] == "CV"
    assert user_info["Obs. Latitude"] == "+31.675467"
    assert user_info["Obs. Longitude"] == "-110.951376"
    assert user_info["Obs. Elevation (meters)"] == 1268
    assert user_info["Pixel Binning"] == "2x2"
    assert user_info["Secondary Observer Codes (N/A if none)"] == "MOBS"


def json_loads(path):
    with open(path) as inits_file:
        return json.load(inits_file)
