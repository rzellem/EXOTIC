import subprocess
import sys
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_normal_import_and_native_transit_do_not_initialize_pylightcurve():
    script = textwrap.dedent(
        """
        import builtins
        import sys
        import numpy as np

        original_import = builtins.__import__
        def forbid_catalogue_initialization(name, *args, **kwargs):
            if name == 'pylightcurve' or name.startswith('pylightcurve.'):
                raise AssertionError('C++ setup must not initialize PyLightcurve catalogues')
            return original_import(name, *args, **kwargs)
        builtins.__import__ = forbid_catalogue_initialization

        import exotic.exotic
        import exotic.api.elca as elca
        values = dict(rprs=.1, per=3., ars=10., ecc=0., inc=89., omega=90.,
                      tmid=0., u0=.3, u1=.2, u2=.1, u3=.05)
        flux = elca.transit(np.array([0., .01]), values)
        assert flux.shape == (2,) and np.isfinite(flux).all()
        assert 'pylightcurve' not in sys.modules
        print('native-startup-ok')
        """
    )
    result = subprocess.run([sys.executable, '-c', script], cwd=REPO_ROOT,
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr or result.stdout
    assert 'native-startup-ok' in result.stdout


def test_python_and_joint_models_need_no_pylightcurve_or_catalogues():
    script = textwrap.dedent(
        """
        import builtins
        import sys
        import numpy as np
        original_import = builtins.__import__
        def forbid_external_package(name, *args, **kwargs):
            if name == 'pylightcurve' or name.startswith('pylightcurve.'):
                raise AssertionError('Python reference must be data-free')
            return original_import(name, *args, **kwargs)
        builtins.__import__ = forbid_external_package
        from exotic.api import elca, joint_fitter
        times = np.linspace(-.1, .1, 31)
        values = dict(_claret_backend='python', rprs=.1, per=3., ars=10.,
                      ecc=.15, inc=89., omega=65., tmid=0.,
                      u0=.3, u1=.2, u2=.1, u3=.05)
        reference = elca.transit(times, values)
        np.testing.assert_array_equal(reference, joint_fitter.transit(times, values))
        assert np.isfinite(reference).all() and reference.min() < 1.
        assert joint_fitter.eclipse_mid_time(3.,10.,0.,89.,90.,0.) == 1.5
        assert not any(n == 'pylightcurve' or n.startswith('pylightcurve.') for n in sys.modules)
        print('local-models-ok')
        """
    )
    result = subprocess.run([sys.executable, '-c', script], cwd=REPO_ROOT,
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr or result.stdout
    assert 'local-models-ok' in result.stdout


def test_transit_supersamples_points_with_significant_exposure_smearing(monkeypatch):
    import numpy as np
    from exotic.api import elca

    def quadratic_transit(*args, **kwargs):
        return 1.0 + np.asarray(args[-1], dtype=float) ** 2

    monkeypatch.setattr(elca, 'pytransit', quadratic_transit)
    values = dict(_claret_backend='python', rprs=.1, per=1., ars=10.,
                  ecc=0., inc=90., omega=90., tmid=0., u0=0., u1=0., u2=0., u3=0.)
    values.update({elca.EXPOSURE_SMEARING_EXPOSURE_TIME_KEY: np.array([3600. / 86400.]),
                   elca.EXPOSURE_SMEARING_SUPERSAMPLE_KEY: 5,
                   elca.EXPOSURE_SMEARING_CHANGE_TOLERANCE_KEY: 0.})
    result = elca.transit(np.array([0.]), values)
    offsets = (np.arange(5, dtype=float) + .5) / 5. - .5
    expected = np.mean(1. + ((3600. / 86400.) * offsets) ** 2)
    np.testing.assert_allclose(result, [expected])
    assert result[0] > 1.


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
            fake_photutils_aperture.CircularAnnulus = type("CircularAnnulus", (), {})
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
