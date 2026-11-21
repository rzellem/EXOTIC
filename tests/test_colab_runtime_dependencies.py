"""Exercise the LDTK JIT and Bokeh image APIs used with retained Colab packages."""
import json

from astropy.io import fits
from bokeh.embed import json_item
from ldtk.ldmodel import NonlinearModel
import numpy as np
import pytest


def test_ldtk_nonlinear_numba_kernel_matches_profile_formula():
    mu = np.linspace(.01, 1., 101)
    coefficients = np.array([1.77248, -2.96806, 3.44175, -1.24617])
    expected = 1. - sum(value * (1. - mu ** ((i + 1) / 2.))
                        for i, value in enumerate(coefficients))
    actual = NonlinearModel.evaluate(mu, coefficients)
    np.testing.assert_allclose(actual[0, 0], expected, rtol=0, atol=3.e-15)


@pytest.mark.parametrize('module_name', ['exotic.api.colab', 'exotic.api.plotting'])
def test_bokeh_image_and_drawing_tools_serialize(module_name, monkeypatch, tmp_path):
    from importlib import import_module
    if module_name == 'exotic.api.colab':
        pytest.importorskip('IPython.display', reason='This helper runs in the Colab/IPython runtime')
    module = import_module(module_name)
    image = (100 + np.arange(32 * 48).reshape(32, 48)).astype(np.uint16)
    path = tmp_path / 'synthetic.fits'
    fits.PrimaryHDU(image).writeto(path)
    figures = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(module, 'show', figures.append)
    if module_name == 'exotic.api.colab':
        module.display_image(str(path))
    else:
        module.plot_image(str(path))
    assert len(figures) == 1
    payload = json.dumps(json_item(figures[0]))
    assert 'FreehandDrawTool' in payload
    assert 'LogColorMapper' in payload
    assert 'image' in payload
