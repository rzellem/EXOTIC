"""Exercise image, plotting and profile APIs with retained Colab packages."""
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


def test_alignment_apis_recover_known_image_translation():
    from exotic import exotic as reduction

    reference = np.random.default_rng(19403).normal(size=(64, 80))
    shifted = np.roll(reference, (3, -5), axis=(0, 1))
    shift, error, _ = reduction.phase_cross_correlation(
        shifted, reference, upsample_factor=4)
    np.testing.assert_array_equal(shift, [3., -5.])
    assert np.isfinite(error)
    transform = reduction.SimilarityTransform(translation=[-shift[1], -shift[0]])
    np.testing.assert_allclose(transform([[10., 20.], [40., 50.]]),
                               [[15., 17.], [45., 47.]], rtol=0, atol=1.e-13)


def test_bundled_imreg_translation_recovers_known_shift():
    from exotic import exotic as reduction
    image = np.random.default_rng(19403).normal(size=(64, 80))
    result = reduction.ird.translation(image, np.roll(image, (3, -5), axis=(0, 1)))
    np.testing.assert_allclose(result['tvec'], [-3., 5.], rtol=0, atol=1.e-12)


@pytest.mark.parametrize('algorithm', ['bilinear', 'malvar2004', 'menon2007'])
def test_slim_demosaicing_preserves_constant_colour_field(algorithm):
    from exotic.api.demosaicing import reconstruct_bayer
    image = np.full((24, 28), 50.)
    image[::2, ::2] = 100.
    image[1::2, 1::2] = 10.
    actual = reconstruct_bayer(image, 'RGGB', algorithm)
    assert actual.shape == (24, 28, 3)
    expected = np.broadcast_to([100., 50., 10.], actual[4:-4, 4:-4].shape)
    np.testing.assert_allclose(actual[4:-4, 4:-4], expected, rtol=0, atol=1.e-12)


@pytest.mark.parametrize('module_name', ['exotic.api.colab', 'exotic.api.plotting'])
def test_display_downsampling_preserves_block_means(module_name):
    from importlib import import_module
    if module_name == 'exotic.api.colab':
        pytest.importorskip('IPython.display', reason='This helper runs in the Colab/IPython runtime')
    module = import_module(module_name)
    image = np.arange(32 * 48, dtype=np.uint16).reshape(32, 48)
    expected = image.reshape(8, 4, 12, 4).mean(axis=(1, 3))
    np.testing.assert_array_equal(module.downscale_local_mean(image, (4, 4)), expected)
