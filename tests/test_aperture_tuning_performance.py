import numpy as np
from astropy.io import fits

import exotic.exotic as exotic_module


def test_evenly_spaced_aperture_tuning_indices_span_full_sequence():
    indices = exotic_module.evenly_spaced_aperture_tuning_indices(524, max_frames=24)

    assert len(indices) == 24
    assert indices[0] == 0
    assert indices[-1] == 523
    assert np.all(np.diff(indices) >= 22)
    assert np.all(np.diff(indices) <= 23)


def test_centered_numpy_cutout_uses_local_slice_coordinates():
    image = np.arange(100, dtype=float).reshape(10, 10)

    cutout, local_x, local_y = exotic_module.centered_numpy_cutout(
        image,
        xc=5.25,
        yc=4.75,
        radius=2.0,
    )

    np.testing.assert_array_equal(cutout, image[2:8, 3:9])
    assert local_x == 2.25
    assert local_y == 2.75


def test_memmap_cutouts_require_identity_frame_processing():
    empty = np.empty((0, 0))

    assert exotic_module.can_memmap_aperture_tuning_cutouts(
        generalDark=empty,
        generalBias=empty,
        generalFlat=empty,
        demosaic_fmt=None,
        bad_pixel_reference=None,
    )
    assert not exotic_module.can_memmap_aperture_tuning_cutouts(
        generalDark=np.ones((2, 2)),
    )
    assert not exotic_module.can_memmap_aperture_tuning_cutouts(
        demosaic_fmt="RGGB",
    )
    assert not exotic_module.can_memmap_aperture_tuning_cutouts(
        bad_pixel_reference={"coord_x": np.array([1]), "coord_y": np.array([1])},
    )


def test_fits_header_memmap_guard_rejects_scaled_images():
    assert exotic_module.fits_header_supports_memmap({"BITPIX": -32})
    assert not exotic_module.fits_header_supports_memmap({"BSCALE": 2.0})
    assert not exotic_module.fits_header_supports_memmap({"BZERO": 32768})


def test_alignment_worker_memmap_path_skips_full_frame_calibration(tmp_path, monkeypatch):
    frame_path = tmp_path / "frame.fits"
    expected = np.arange(100, dtype=np.float32).reshape(10, 10)
    fits.PrimaryHDU(expected).writeto(frame_path)
    empty = np.empty((0, 0))
    exotic_module._ALIGNMENT_POOL_CONTEXT = {
        "generalDark": empty,
        "generalBias": empty,
        "generalFlat": empty,
        "demosaic_fmt": None,
        "demosaic_out": None,
        "demosaic_mult": None,
        "bad_pixel_reference": None,
    }
    monkeypatch.setattr(
        exotic_module,
        "apply_cals",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("identity memmap path must not calibrate the full image")
        ),
    )

    header, image = exotic_module._load_alignment_worker_frame(str(frame_path))

    assert header["NAXIS1"] == 10
    np.testing.assert_array_equal(image[2:5, 3:7], expected[2:5, 3:7])


def test_aperture_tuning_cutout_grid_matches_direct_local_measurement():
    y, x = np.indices((41, 41), dtype=float)
    image = 100.0 + 5000.0 * np.exp(-((x - 20.2) ** 2 + (y - 19.8) ** 2) / (2.0 * 2.0 ** 2))
    frame = {
        "frame_sigma": 2.0,
        "stars": {
            "comp1": {
                "data": image,
                "xc": 20.2,
                "yc": 19.8,
                "sigma": 2.0,
            }
        },
    }
    apertures = np.array([2.5, 3.0])
    annuli = np.array([6.0, 8.0])

    result = exotic_module.populate_aperture_tuning_data_from_cutouts(
        [frame],
        apertures,
        annuli,
        comparison_indices=(0,),
        adaptive_apertures=True,
        reference_sigma=2.0,
    )
    expected_flux, expected_bg = exotic_module.compute_star_aperture_grid(
        image,
        1,
        20.2,
        19.8,
        apertures * 2.0,
        annuli * 2.0,
        sigma_hint=2.0,
    )

    np.testing.assert_allclose(result["comp1"][0], expected_flux)
    np.testing.assert_allclose(result["comp1_bg"][0], expected_bg)


def test_image_process_pool_uses_spawn_context_on_windows(monkeypatch):
    captured = {}

    class FakeProcessPool:
        def __init__(self, *args, **kwargs):
            captured["args"] = args
            captured["kwargs"] = kwargs

    fake_context = object()
    monkeypatch.setattr(exotic_module.sys, "platform", "win32")
    monkeypatch.setattr(exotic_module.multiprocessing, "get_context", lambda mode: fake_context)
    monkeypatch.setattr(exotic_module, "_ProcessPoolExecutor", FakeProcessPool)

    executor = exotic_module.ImageProcessPoolExecutor(max_workers=3)

    assert isinstance(executor, FakeProcessPool)
    assert captured["kwargs"]["max_workers"] == 3
    assert captured["kwargs"]["mp_context"] is fake_context
