import pytest

import exotic.exotic as exotic_module


@pytest.mark.parametrize('telescope', ['Ben', 'Donald', ' donald '])
@pytest.mark.parametrize('style', ['original', 'lower', 'padded'])
@pytest.mark.parametrize('observed_filter, code, band, wl_min, wl_max', [
    ('Red', 'R', 'Cousins R', 561.7, 719.7),
    ('Green', 'V', 'Johnson V', 502.8, 586.8),
    ('Blue', 'B', 'Johnson B', 391.6, 480.6),
    ('R', 'R', 'Cousins R', 561.7, 719.7),
])
def test_microobservatory_colors_resolve_to_standard_bands_without_prompt(
        monkeypatch, telescope, style, observed_filter, code, band, wl_min, wl_max):
    if style == 'lower':
        observed_filter = observed_filter.lower()
    elif style == 'padded':
        observed_filter = f' {observed_filter} '
    ld = exotic_module.LimbDarkening({})
    monkeypatch.setattr(ld, 'calculate_ld', lambda: None)
    monkeypatch.setattr(
        exotic_module, 'user_input',
        lambda *_args, **_kwargs: pytest.fail('Ben/Donald filters must resolve without a prompt'),
    )
    info = {
        'filter': observed_filter, 'wl_min': None, 'wl_max': None,
        'second_obs': 'MOBS',
    }
    exotic_module.apply_microobservatory_metadata(info, {
        'CREATOR': 'MicroObservatory2010-11-04', 'TELESCOP': telescope, 'FILTER': observed_filter,
    })
    exotic_module.nonlinear_ld(ld, info, non_interactive_run=True)

    assert info['filter'] == code
    assert info['filter_desc'] == band
    assert info['wl_min'] == wl_min
    assert info['wl_max'] == wl_max


@pytest.mark.parametrize('observed_filter', ['Red', 'Green', 'Blue'])
def test_other_telescope_color_label_is_preserved(observed_filter):
    info = {'filter': observed_filter, 'second_obs': 'MOBS'}
    exotic_module.apply_microobservatory_metadata(info, {
        'CREATOR': 'MicroObservatory2010-11-04', 'TELESCOP': 'Other', 'FILTER': observed_filter,
    })
    assert info['filter'] == observed_filter


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


def test_nonlinear_ld_non_interactive_uses_neutral_cbb_name_without_prompt(monkeypatch):
    ld = exotic_module.LimbDarkening({})
    monkeypatch.setattr(ld, "calculate_ld", lambda: None)
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda *_args, **_kwargs: pytest.fail("recognized CBB filter must not prompt"),
    )
    info_dict = {
        "filter": "CBB",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": "y",
    }

    exotic_module.nonlinear_ld(ld, info_dict, non_interactive_run=True)

    assert info_dict["filter"] == "CBB"
    assert info_dict["filter_desc"] == "CBB"
    assert info_dict["wl_min"] == 500.0
    assert info_dict["wl_max"] == 1000.0


def test_nonlinear_ld_non_interactive_treats_cv_as_clearv_without_prompt(monkeypatch):
    ld = exotic_module.LimbDarkening({})
    monkeypatch.setattr(ld, "calculate_ld", lambda: None)
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda *_args, **_kwargs: pytest.fail("recognized CV filter must not prompt"),
    )
    info_dict = {
        "filter": "CV",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": "y",
    }

    exotic_module.nonlinear_ld(ld, info_dict, non_interactive_run=True)

    assert info_dict["filter"] == "CV"
    assert info_dict["filter_desc"] == "CV"
    assert info_dict["wl_min"] == 350.0
    assert info_dict["wl_max"] == 1000.0


class UnrecognizedFilterLimbDarkening:
    fwhm_names_nonspecific = {}

    @staticmethod
    def check_fwhm(_observed_filter):
        return False

    @staticmethod
    def check_standard(_observed_filter):
        return False


class BooleanOptionLimbDarkening(UnrecognizedFilterLimbDarkening):
    filter_name = None
    filter_desc = None
    wl_min = None
    wl_max = None

    def calculate_ld(self):
        return None


def set_boolean_option_filter(ld, label):
    ld.filter_name = label
    ld.filter_desc = label
    ld.wl_min = 400.0
    ld.wl_max = 700.0


@pytest.mark.parametrize("config_value", [True, 1, "1", "y", "Y", "yes", "TRUE", "on"])
def test_nonlinear_ld_boolean_option_accepts_true_forms(monkeypatch, config_value):
    ld = BooleanOptionLimbDarkening()
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda prompt, **_kwargs: 1 if "enter 1" in prompt.lower() else pytest.fail(
            "valid true boolean must not trigger the y/n prompt"
        ),
    )
    monkeypatch.setattr(
        exotic_module,
        "standard_filter",
        lambda selected_ld, _observed_filter: set_boolean_option_filter(selected_ld, "standard"),
    )
    monkeypatch.setattr(
        exotic_module,
        "user_entered_ld",
        lambda *_args, **_kwargs: pytest.fail("true must select calculated limb darkening"),
    )
    info_dict = {
        "filter": "mystery-band",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": config_value,
    }

    exotic_module.nonlinear_ld(ld, info_dict)

    assert info_dict["filter"] == "standard"


@pytest.mark.parametrize("config_value", [False, 0, "0", "n", "N", "no", "FALSE", "off"])
def test_nonlinear_ld_boolean_option_accepts_false_forms(monkeypatch, config_value):
    ld = BooleanOptionLimbDarkening()
    monkeypatch.setattr(
        exotic_module,
        "user_input",
        lambda *_args, **_kwargs: pytest.fail("valid false boolean must not prompt"),
    )
    monkeypatch.setattr(
        exotic_module,
        "standard_filter",
        lambda *_args, **_kwargs: pytest.fail("false must select user-entered limb darkening"),
    )
    monkeypatch.setattr(
        exotic_module,
        "user_entered_ld",
        lambda selected_ld, _observed_filter: set_boolean_option_filter(selected_ld, "manual"),
    )
    info_dict = {
        "filter": "mystery-band",
        "wl_min": None,
        "wl_max": None,
        "ld_uncertainties": config_value,
    }

    exotic_module.nonlinear_ld(ld, info_dict)

    assert info_dict["filter"] == "manual"


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
