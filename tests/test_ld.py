import logging

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

def test_missing_fwhm_values_do_not_log_errors(caplog) -> None:
    observed_filter = {
        'filter': None,
        'name': None,
        'wl_min': None,
        'wl_max': None
    }

    ld_obj = LimbDarkening(stellar_params)

    with caplog.at_level(logging.ERROR, logger="exotic.api.ld"):
        assert ld_obj.check_fwhm(observed_filter) == False

    assert "FWHM matching failed" not in caplog.text

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


def test_photographic_filter_aliases_in_filter_column() -> None:
    alias_cases = [
        ("pb", "Photographic B", "PB", "391.6", "480.6"),
        ("pg", "Photographic G", "PG", "502.8", "586.8"),
        ("pr", "Photographic R", "PR", "590.0", "810.0"),
    ]

    for alias, expected_filter, expected_name, expected_min, expected_max in alias_cases:
        observed_filter = {'filter': alias, 'name': None, 'wl_min': None, 'wl_max': None}
        setting_filter_values(observed_filter)
        assert observed_filter == {
            'filter': expected_filter,
            'name': expected_name,
            'wl_min': expected_min,
            'wl_max': expected_max,
        }


def test_additional_standard_filter_aliases_in_filter_column() -> None:
    alias_cases = [
        ("bu", "Johnson U", "U", "333.8", "398.8"),
        ("bi", "Johnson I", "IJ", "780.0", "1020.0"),
        ("up", "Sloan u", "SU", "321.8", "386.8"),
        ("gp", "Sloan g", "SG", "402.5", "551.5"),
        ("rp", "Sloan r", "SR", "553.1", "693.1"),
        ("ip", "Sloan i", "SI", "697.5", "827.5"),
        ("zp", "Sloan z", "SZ", "841.2", "978.2"),
        ("su", "Stromgren u", "STU", "336.3", "367.7"),
        ("sv", "Stromgren v", "STV", "401.5", "418.5"),
        ("sb", "Stromgren b", "STB", "459.55", "478.05"),
        ("sy", "Stromgren y", "STY", "536.7", "559.3"),
        ("hb", "Stromgren Hbw", "STHBW", "481.5", "496.5"),
        ("zs", "PanSTARRS z-short", "ZS", "826.0", "920.0"),
        ("clearV", "MObs CV", "CV", "350.0", "850.0"),
        ("w", "MObs CV", "CV", "350.0", "850.0"),
        ("pl", "MObs CV", "CV", "350.0", "850.0"),
        ("exo", "Astrodon ExoPlanet-BB", "CBB", "500.0", "1000.0"),
        ("Astrodon-Exo", "Astrodon ExoPlanet-BB", "CBB", "500.0", "1000.0"),
    ]

    for alias, expected_filter, expected_name, expected_min, expected_max in alias_cases:
        observed_filter = {'filter': alias, 'name': None, 'wl_min': None, 'wl_max': None}
        setting_filter_values(observed_filter)
        assert observed_filter == {
            'filter': expected_filter,
            'name': expected_name,
            'wl_min': expected_min,
            'wl_max': expected_max,
        }


def test_osc_split_filter_aliases_in_filter_column() -> None:
    alias_cases = [
        ("B1", "Photographic B", "PB", "391.6", "480.6"),
        ("G1", "Photographic G", "PG", "502.8", "586.8"),
        ("G2", "Photographic G", "PG", "502.8", "586.8"),
        ("R1", "Photographic R", "PR", "590.0", "810.0"),
        ("R2", "Photographic R", "PR", "590.0", "810.0"),
    ]

    for alias, expected_filter, expected_name, expected_min, expected_max in alias_cases:
        observed_filter = {'filter': alias, 'name': None, 'wl_min': None, 'wl_max': None}
        setting_filter_values(observed_filter)
        assert observed_filter == {
            'filter': expected_filter,
            'name': expected_name,
            'wl_min': expected_min,
            'wl_max': expected_max,
        }
