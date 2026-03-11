from exotic.output_files import OutputFiles, save_comp_star_calibration_summary


class DummyFit:
    def __init__(self):
        self.parameters = {
            "tmid": 2450000.123456,
            "rprs": 0.1234,
            "inc": 88.5,
            "a1": 1.0,
            "a2": 0.0,
        }
        self.errors = {
            "tmid": 0.0001,
            "rprs": 0.001,
            "inc": 0.2,
            "a1": 0.1,
            "a2": 0.1,
        }
        self.time = [2450000.123456]
        self.data = [1.0]
        self.dataerr = [0.01]
        self.airmass_model = [1.0]


def test_aavso_output_includes_observatory_location_headers(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "Whipple Observatory",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "#OBSDATE=2020-01-01" in output_text
    assert "#OBSNAME=Whipple Observatory" in output_text
    assert "#OBSLAT=+32.41638889" in output_text
    assert "#OBSLON=-110.73444444" in output_text
    assert "#OBSELEV=2616" in output_text


def test_aavso_output_omits_obsname_header_when_blank(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "#OBSNAME=" not in output_text


def test_save_comp_star_calibration_summary_writes_selected_star(tmp_path):
    summary_path = save_comp_star_calibration_summary(
        tmp_path,
        "HAT-P-32 b",
        "2026-03-09",
        "PSF photometry",
        0.0012,
        [
            {
                "label": "Comp 1",
                "position": [101, 202],
                "selected": True,
                "aggregate_score": 0.0012,
                "ensemble_score": 0.0010,
                "pairwise_median_score": 0.0011,
                "pairwise_max_score": 0.0014,
                "self_score": 0.0009,
                "valid_pair_count": 2,
            },
            {
                "label": "Comp 2",
                "position": [303, 404],
                "selected": False,
                "aggregate_score": 0.0031,
                "ensemble_score": 0.0028,
                "pairwise_median_score": 0.0030,
                "pairwise_max_score": 0.0035,
                "self_score": 0.0012,
                "valid_pair_count": 2,
            },
        ],
        0,
    )

    text = summary_path.read_text()
    assert "# Selected comparison star,1" in text
    assert "Comp 1,101,202,true" in text
