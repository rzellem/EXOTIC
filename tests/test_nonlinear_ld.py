import pytest

import exotic.exotic as exotic_module


def test_nonlinear_ld_non_interactive_treats_g_as_photographic_g_without_prompt(monkeypatch):
    ld = exotic_module.LimbDarkening({})
    monkeypatch.setattr(ld, "calculate_ld", lambda: None)
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda *_args, **_kwargs: pytest.fail("recognized G filter must not prompt"),
    )
    info_dict = {
        "filter": "G",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": "y",
    }

    exotic_module.nonlinear_ld(ld, info_dict, non_interactive_run=True)

    assert info_dict["filter"] == "PG"
    assert info_dict["filter_desc"] == "Photographic G"
    assert info_dict["wl_min"] == 502.8
    assert info_dict["wl_max"] == 586.8


class UnrecognizedFilterLimbDarkening:
    fwhm_names_nonspecific = {}

    @staticmethod
    def check_fwhm(_observed_filter):
        return False

    @staticmethod
    def check_standard(_observed_filter):
        return False


def test_nonlinear_ld_non_interactive_rejects_unrecognized_filter_without_prompt(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda *_args, **_kwargs: pytest.fail("non-interactive limb-darkening selection must not prompt"),
    )
    info_dict = {
        "filter": "mystery-band",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": "y",
    }

    with pytest.raises(ValueError, match="did not recognize the filter 'mystery-band'"):
        exotic_module.nonlinear_ld(
            UnrecognizedFilterLimbDarkening(),
            info_dict,
            non_interactive_run=True,
        )
