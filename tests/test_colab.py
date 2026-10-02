"""Colab helper: MicroObservatory site metadata and header-driven binning.

Ported from PR #1382 (S-Bhattacharya240611) onto the magic branch, with the
pixel binning read from the frame header instead of a fixed template value.
"""
import importlib
import json
import sys
import types

import numpy as np
import pytest
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


def _mobs_frame(tmp_path, binning=(2, 2)):
    fits_path = tmp_path / "mobs.fits"
    fits.PrimaryHDU(data=np.zeros((2, 2))).writeto(fits_path)
    with fits.open(fits_path, mode="update") as hdul:
        hdr = hdul[0].header
        hdr["OBSERVAT"] = "Whipple Observatory"
        hdr["DATE"] = "2017-12-20T01:33:43.000"
        hdr["WEATHER"] = 80
        hdr["TELTEMP"] = 20.0
        hdr["CAMTEMP"] = 10.0
        if binning is not None:
            hdr["XBINNING"], hdr["YBINNING"] = binning
    return fits_path


def test_mobs_defaults_match_exotic_metadata(monkeypatch, tmp_path):
    colab = import_colab(monkeypatch)
    fits_path = _mobs_frame(tmp_path)

    hdr = fits.getheader(fits_path)
    assert colab.find(hdr, ["FILTER", "FILT"], "MObs") == "CV"
    assert colab.find(hdr, ["LATITUDE", "LAT", "SITELAT"], "MObs") == "+31.675467"
    assert colab.find(hdr, ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"], "MObs") == "-110.951376"
    assert colab.find(hdr, ["HEIGHT", "ELEVATION", "ELE", "EL", "OBSGEO-H", "ALT-OBS", "SITEELEV"], "MObs") == 1268

    inits_path = colab.make_inits_file(
        '"planetary_parameters": {}',
        tmp_path.as_posix(),
        f"{tmp_path.as_posix()}/",
        fits_path.as_posix(),
        "[424, 286]",
        "[[465, 183], [512, 263]]",
        "MObs",
        "RTZ",
        "",
        False,
    )
    with open(inits_path) as inits_file:
        user_info = json.load(inits_file)["user_info"]

    assert user_info["Filter Name (aavso.org/filters)"] == "CV"
    assert user_info["Obs. Latitude"] == "+31.675467"
    assert user_info["Obs. Longitude"] == "-110.951376"
    assert user_info["Obs. Elevation (meters)"] == 1268
    assert user_info["Pixel Binning"] == "2x2"
    assert user_info["Secondary Observer Codes (N/A if none)"] == "MOBS"


@pytest.mark.parametrize(
    "cards, expected",
    [
        ({"XBINNING": 2, "YBINNING": 2}, "2x2"),
        ({"XBINNING": 1, "YBINNING": 1}, "1x1"),
        ({"XBINNING": "3"}, "3x3"),
        ({"CCDXBIN": 2, "CCDYBIN": 4}, "2x4"),
        ({"BINNING": "2X2"}, "2x2"),
        ({"XBINNING": "junk"}, "1x1"),
        ({}, "1x1"),
    ],
)
def test_pixel_binning_from_header(monkeypatch, cards, expected):
    colab = import_colab(monkeypatch)
    hdr = fits.Header()
    for key, value in cards.items():
        hdr[key] = value
    assert colab.pixel_binning_from_header(hdr) == expected
