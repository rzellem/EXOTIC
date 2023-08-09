from exotic.api.ld import LimbDarkening

stellar_params = {
        'teff': 6001.0,
        'teffUncPos': 88.0,
        'teffUncNeg': -88.0,
        'met': -0.16,
        'metUncPos': 0.08,
        'metUncNeg': -0.08,
        'logg': 4.22,
        'loggUncPos': 0.04,
        'loggUncNeg': -0.04
    }

def setting_filter_values(observed_filter) -> None:
    ld_obj = LimbDarkening(stellar_params)
    ld_obj.check_standard(observed_filter)

def test_existing_standard_filter_abbreviation() -> None:
    observed_filter = {
        'filter': "SU",
        'name': None,
        'wl_min': None,
        'wl_max': None
    }

    expected_filter = {
        'filter': "Sloan u",
        'name': 'SU',
        'wl_min': '321.8',
        'wl_max': '386.8'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_existing_standard_filter_name() -> None:
    observed_filter = {
        'filter': "Optec Wing A",
        'name': None,
        'wl_min': None,
        'wl_max': None
    }

    expected_filter = {
        'filter': "Optec Wing A",
        'name': 'MA',
        'wl_min': '706.5',
        'wl_max': '717.5'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_existing_standard_filter_alias_name() -> None:
    observed_filter = {
        'filter': "LCO Bessell B",
        'name': None,
        'wl_min': None,
        'wl_max': None
    }

    expected_filter = {
        'filter': "Johnson B",
        'name': 'B',
        'wl_min': '391.6',
        'wl_max': '480.6'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_existing_mobs_standard_filter_name() -> None:
    observed_filter = {
        'filter': "MObs CV",
        'name': None,
        'wl_min': '300.0',
        'wl_max': '800.0'
    }

    expected_filter = {
        'filter': "MObs CV",
        'name': 'CV',
        'wl_min': '350.0',
        'wl_max': '850.0'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_custom_nonspecific_standard_filter_abbreviation_1() -> None:
    observed_filter = {
        'filter': "CV",
        'name': None,
        'wl_min': '300.0',
        'wl_max': '800.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_standard(observed_filter) == False

def test_custom_nonspecific_standard_filter_abbreviation_2() -> None:
    observed_filter = {
        'filter': "NA",
        'name': None,
        'wl_min': '500.0',
        'wl_max': '700.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_standard(observed_filter) == False

def test_nonexisting_standard_filter_abbreviation() -> None:
    observed_filter = {
        'filter': "TF",
        'name': None,
        'wl_min': '600.0',
        'wl_max': '700.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_standard(observed_filter) == False

def test_existing_standard_filter_fwhm() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '333.8',
        'wl_max': '398.8'
    }

    expected_filter = {
        'filter': "Johnson U",
        'name': 'U',
        'wl_min': '333.8',
        'wl_max': '398.8'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_existing_mobs_standard_filter_mobs() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '350.0',
        'wl_max': '850.0'
    }

    expected_filter = {
        'filter': "MObs CV",
        'name': 'CV',
        'wl_min': '350.0',
        'wl_max': '850.0'
    }

    setting_filter_values(observed_filter)

    assert observed_filter == expected_filter

def test_valid_fwhm_range() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '350.0',
        'wl_max': '850.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_fwhm(observed_filter) == True

def test_valid_fwhm_range_swapped_min_max() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '400.0',
        'wl_max': '200.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_fwhm(observed_filter) == True

def test_invalid_fwhm_range_1() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '100.0',
        'wl_max': '1000.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_fwhm(observed_filter) == False

def test_invalid_fwhm_range_2() -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': '-100.0',
        'wl_max': '3000.0'
    }

    ld_obj = LimbDarkening(stellar_params)

    assert ld_obj.check_fwhm(observed_filter) == False
