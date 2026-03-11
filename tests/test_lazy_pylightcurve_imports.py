import subprocess
import sys
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_imports_eagerly_load_pylightcurve_without_noise():
    script = textwrap.dedent(
        """
        import tempfile
        import sys
        import types
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            package_dir = root / "pylightcurve"
            model_dir = package_dir / "models"
            model_dir.mkdir(parents=True)

            (package_dir / "__init__.py").write_text("", encoding="utf-8")
            (model_dir / "__init__.py").write_text("", encoding="utf-8")
            (model_dir / "exoplanet_lc.py").write_text(
                "import sys\\n"
                "print('LOUD-STDOUT')\\n"
                "print('LOUD-STDERR', file=sys.stderr)\\n"
                "def transit(*args, **kwargs): return 'stub-transit'\\n"
                "def eclipse_mid_time(*args, **kwargs): return 0.0\\n",
                encoding="utf-8",
            )

            sys.path.insert(0, str(root))

            for name in list(sys.modules):
                if name == "exotic.api.elca" or name == "exotic.api.joint_fitter" or name.startswith("pylightcurve"):
                    sys.modules.pop(name)

            fake_ultranest = types.ModuleType("ultranest")
            fake_ultranest.ReactiveNestedSampler = type("ReactiveNestedSampler", (), {})
            sys.modules["ultranest"] = fake_ultranest
            fake_plotting = types.ModuleType("plotting")
            fake_plotting.corner = lambda *args, **kwargs: None
            sys.modules["plotting"] = fake_plotting
            sys.modules["exotic.api.plotting"] = fake_plotting
            fake_ultranest_utils = types.ModuleType("ultranest_utils")
            fake_ultranest_utils.run_reactive_sampler = lambda *args, **kwargs: None
            sys.modules["ultranest_utils"] = fake_ultranest_utils
            sys.modules["exotic.api.ultranest_utils"] = fake_ultranest_utils

            import exotic.api.elca as elca
            import exotic.api.joint_fitter as joint_fitter

            assert "pylightcurve.models.exoplanet_lc" in sys.modules

            minimal_values = {
                "u0": 0.0,
                "u1": 0.0,
                "u2": 0.0,
                "u3": 0.0,
                "rprs": 0.1,
                "per": 1.0,
                "ars": 10.0,
                "ecc": 0.0,
                "inc": 89.0,
                "omega": 90.0,
                "tmid": 0.0,
            }

            assert elca.transit([0.0], minimal_values) == "stub-transit"
            assert joint_fitter.pytransit([0.0, 0.0, 0.0, 0.0], 0.1, 1.0, 10.0, 0.0, 89.0, 90.0, 0.0, [0.0]) == "stub-transit"
            print("imports-ok")
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr or result.stdout
    assert "imports-ok" in result.stdout
    assert result.stdout.count("Importing modules. Please wait.......") == 1
    assert "LOUD-STDOUT" not in result.stdout
    assert "LOUD-STDERR" not in result.stderr


def test_import_exotic_avoids_unused_astroquery_modules():
    script = textwrap.dedent(
        """
        import sys
        import tempfile
        import types
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            astroquery_dir = root / "astroquery"
            astroquery_dir.mkdir(parents=True)
            (astroquery_dir / "__init__.py").write_text("", encoding="utf-8")
            (astroquery_dir / "simbad.py").write_text(
                "import sys\\n"
                "print('LOUD-SIMBAD-STDOUT')\\n"
                "print('LOUD-SIMBAD-STDERR', file=sys.stderr)\\n"
                "class Simbad:\\n    pass\\n",
                encoding="utf-8",
            )
            (astroquery_dir / "gaia.py").write_text(
                "import sys\\n"
                "print('LOUD-GAIA-STDOUT')\\n"
                "print('LOUD-GAIA-STDERR', file=sys.stderr)\\n"
                "class Gaia:\\n    pass\\n",
                encoding="utf-8",
            )
            sys.path.insert(0, str(root))

            fake_barycorrpy = types.ModuleType("barycorrpy")
            fake_utc_tdb = types.ModuleType("barycorrpy.utc_tdb")
            fake_utc_tdb.JDUTC_to_BJDTDB = lambda *args, **kwargs: None
            fake_astroalign = types.ModuleType("astroalign")
            fake_astroalign.PIXEL_TOL = 1
            fake_imreg_dft = types.ModuleType("imreg_dft")
            fake_colour_demosaicing = types.ModuleType("colour_demosaicing")
            fake_colour_demosaicing.demosaicing_CFA_Bayer_bilinear = lambda *args, **kwargs: None
            fake_photutils = types.ModuleType("photutils")
            fake_photutils_aperture = types.ModuleType("photutils.aperture")
            fake_photutils_aperture.CircularAperture = type("CircularAperture", (), {})
            fake_photutils_detection = types.ModuleType("photutils.detection")
            fake_photutils_detection.DAOStarFinder = type("DAOStarFinder", (), {})
            fake_ldtk = types.ModuleType("ldtk")
            fake_ldtk.LDPSet = type("LDPSet", (), {})
            fake_ldtk.ldtk = types.SimpleNamespace(LDPSet=fake_ldtk.LDPSet)
            fake_ldtk_ldmodel = types.ModuleType("ldtk.ldmodel")
            fake_ldtk_ldmodel.LinearModel = type("LinearModel", (), {})
            fake_ldtk_ldmodel.QuadraticModel = type("QuadraticModel", (), {})
            fake_ldtk_ldmodel.NonlinearModel = type("NonlinearModel", (), {})
            fake_lmfit = types.ModuleType("lmfit")
            fake_pyvo = types.ModuleType("pyvo")
            fake_ultranest = types.ModuleType("ultranest")
            fake_ultranest.ReactiveNestedSampler = type("ReactiveNestedSampler", (), {})
            fake_elca = types.ModuleType("exotic.api.elca")
            fake_elca.lc_fitter = lambda *args, **kwargs: None
            fake_elca.binner = lambda *args, **kwargs: None
            fake_elca.transit = lambda *args, **kwargs: None
            fake_elca.get_phase = lambda *args, **kwargs: None
            fake_ld = types.ModuleType("exotic.api.ld")
            fake_ld.LimbDarkening = type("LimbDarkening", (), {})
            fake_ld.ld_re_punct_p = lambda *args, **kwargs: None

            sys.modules.setdefault("astroalign", fake_astroalign)
            sys.modules.setdefault("barycorrpy", fake_barycorrpy)
            sys.modules.setdefault("barycorrpy.utc_tdb", fake_utc_tdb)
            sys.modules.setdefault("imreg_dft", fake_imreg_dft)
            sys.modules.setdefault("colour_demosaicing", fake_colour_demosaicing)
            sys.modules.setdefault("photutils", fake_photutils)
            sys.modules.setdefault("photutils.aperture", fake_photutils_aperture)
            sys.modules.setdefault("photutils.detection", fake_photutils_detection)
            sys.modules.setdefault("ldtk", fake_ldtk)
            sys.modules.setdefault("ldtk.ldmodel", fake_ldtk_ldmodel)
            sys.modules.setdefault("lmfit", fake_lmfit)
            sys.modules.setdefault("pyvo", fake_pyvo)
            sys.modules.setdefault("ultranest", fake_ultranest)
            sys.modules.setdefault("exotic.api.elca", fake_elca)
            sys.modules.setdefault("exotic.api.ld", fake_ld)

            import exotic.exotic

            assert "astroquery.gaia" not in sys.modules
            assert "astroquery.simbad" not in sys.modules
            print("exotic-import-ok")
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr or result.stdout
    assert "exotic-import-ok" in result.stdout
    assert "LOUD-SIMBAD-STDOUT" not in result.stdout
    assert "LOUD-SIMBAD-STDERR" not in result.stderr
    assert "LOUD-GAIA-STDOUT" not in result.stdout
    assert "LOUD-GAIA-STDERR" not in result.stderr
