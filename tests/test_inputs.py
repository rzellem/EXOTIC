import json

from exotic.inputs import Inputs, camera


def test_camera_accepts_cmos_as_ccd_without_prompt():
    assert camera("CMOS") == "CCD"


def test_camera_defaults_to_ccd_when_missing_or_unrecognized():
    assert camera(None) == "CCD"
    assert camera("") == "CCD"
    assert camera("mirrorless") == "CCD"


def test_camera_keeps_dslr_as_dslr():
    assert camera("DSLR") == "DSLR"
    assert camera("canon dslr") == "DSLR"


def test_comp_params_accepts_verbose_camera_key(tmp_path):
    init_data = {
        "user_info": {
            "Camera Type (e.g., CCD or DSLR; Note: if you are using a CMOS, please enter CCD here and then note your actual camera type in \"Observing Notes\")": "CCD"
        },
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["camera"] == "CCD"


def test_comp_params_defaults_require_comp_star_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["require_comp_star"] == "y"


def test_comp_params_defaults_ignore_header_wcs_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ignore_header_wcs"] == "n"


def test_comp_params_reads_observatory_full_title_from_user_info(tmp_path):
    init_data = {
        "user_info": {"Observatory Full Title": "Whipple Observatory"},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["obs_name"] == "Whipple Observatory"


def test_comp_params_reads_require_comp_star_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"require_comp_star": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["require_comp_star"] == "n"


def test_comp_params_reads_ignore_header_wcs_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"Ignore WCS in Header and Do Manual Alignment? (y/n)": "y"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ignore_header_wcs"] == "y"


def test_prereduced_mode_forces_aavso_comp_to_no(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text("time flux uncertainty\n")

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": {"ra": "", "dec": "", "x": "", "y": ""},
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["aavso_comp"] == "n"


def test_prereduced_allows_blank_observatory_location_for_bjd_tdb(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text("time flux uncertainty\n")

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "",
        "long": "",
        "elev": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["lat"] is None
    assert info_dict["long"] is None
    assert info_dict["elev"] is None


def test_comp_params_accepts_blank_if_none_phot_comp_star_key(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {
            "Comparison Star used in Photometry (blank if none)": {
                "ra": "",
                "dec": "",
                "x": "493",
                "y": "202",
            }
        },
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "493", "y": "202"}


def test_prereduced_uses_aavso_comp_star_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#COMP_STAR-XC={\"ra\": null, \"dec\": null, \"x\": \"493\", \"y\": \"202\"}\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "493", "y": "202"}


def test_prereduced_uses_aavso_observatory_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSLAT=+32.41638889\n"
        "#OBSLON=-110.73444444\n"
        "#OBSELEV=2616\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "",
        "long": "",
        "elev": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["lat"] == 32.41638889
    assert info_dict["long"] == -110.73444444
    assert info_dict["elev"] == 2616.0


def test_prereduced_uses_aavso_obsdate_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSDATE=2026-03-08\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2026-03-08"


def test_prereduced_prefers_aavso_obsdate_metadata_over_init_date(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSDATE=2026-03-08\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "1999-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2026-03-08"


def test_prereduced_derives_obsdate_from_first_data_row_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text(
        "time,flux,uncertainty\n"
        "2458849.5,0.979108,0.0386426\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "JD_UTC",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2020-01-01"


def test_prereduced_leaves_phot_comp_star_blank_when_missing_from_aavso_metadata(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "", "y": ""}
