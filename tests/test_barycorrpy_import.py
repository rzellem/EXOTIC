"""Exercise the real startup dependency with an isolated Astropy cache."""
import subprocess
import sys
import textwrap
from pathlib import Path
from types import ModuleType

import pytest


@pytest.mark.parametrize('cache_api', ['legacy', 'modern'])
@pytest.mark.parametrize('cached', [False, True])
def test_bootstrap_supports_both_astropy_cache_apis(monkeypatch, cache_api, cached):
    from exotic import earth_orientation as orientation

    lookups = []
    imports = []
    converter = object()
    utc_tdb = ModuleType('barycorrpy.utc_tdb')
    utc_tdb.JDUTC_to_BJDTDB = converter
    monkeypatch.setitem(sys.modules, 'barycorrpy.utc_tdb', utc_tdb)

    # Strict signatures reproduce the pre-8.0 TypeError and require the quiet
    # missing-cache behavior on 8.0+, without masking any lookup exception.
    def legacy_lookup(url_key, pkgname='astropy'):
        lookups.append(url_key)
        return cached

    def modern_lookup(url_key, pkgname='astropy', *, on_missing):
        assert on_missing == 'ignore'
        lookups.append(url_key)
        return cached

    def import_reference(url_key, filename, *, replace):
        imports.append((url_key, filename, replace))

    monkeypatch.setattr(orientation, 'is_url_in_cache',
                        legacy_lookup if cache_api == 'legacy' else modern_lookup)
    monkeypatch.setattr(orientation, 'import_file_to_cache', import_reference)

    assert orientation.load_barycentric_converter() is converter
    assert lookups == [orientation.iers.IERS_B_URL]
    assert imports == ([] if cached else [
        (orientation.iers.IERS_B_URL, orientation.iers.IERS_B_FILE, False)])


def test_bootstrap_does_not_hide_cache_lookup_errors(monkeypatch):
    from exotic import earth_orientation as orientation

    def failed_lookup(url_key, *, on_missing):
        raise TypeError('cache lookup failed internally')

    monkeypatch.setattr(orientation, 'is_url_in_cache', failed_lookup)
    with pytest.raises(TypeError, match='cache lookup failed internally'):
        orientation.load_barycentric_converter()


@pytest.mark.parametrize('cached', [False, True])
@pytest.mark.parametrize('auto_download', [False, True])
def test_import_uses_local_iers_b_without_network(tmp_path, cached, auto_download):
    script = textwrap.dedent(
        """
        import logging
        import sys
        import warnings
        from pathlib import Path
        from unittest.mock import patch

        from astropy.config.paths import set_temp_cache
        from astropy.utils import data, iers

        cache, cached, auto_download = sys.argv[1:]
        cached = cached == 'True'
        auto_download = auto_download == 'True'
        logging.basicConfig(level=logging.INFO)

        with set_temp_cache(cache):
            if cached:
                cached_file = data.import_file_to_cache(iers.IERS_B_URL, iers.IERS_B_FILE)
                cached_mtime = Path(cached_file).stat().st_mtime_ns

            # A brand-new runtime has no download cache directory yet. Treat
            # that startup warning as a failure rather than overlooking it.
            warnings.simplefilter('error', data.CacheMissingWarning)

            def no_download(*args, **kwargs):
                raise OSError('IERS publisher unavailable')

            with iers.conf.set_temp('auto_download', auto_download), \
                    iers.conf.set_temp('auto_max_age', 30), \
                    patch.object(data, '_download_file_from_source', no_download):
                # Exercise the bootstrap before other runtime dependencies can
                # create the cache as an import side effect.
                from exotic.earth_orientation import load_barycentric_converter
                converter = load_barycentric_converter()
                import ldtk
                import exotic.exotic as reduction
                import barycorrpy.PINT_erfautils as pint

                assert reduction.JDUTC_to_BJDTDB is converter
                assert len(pint.iers_b) > 0
                assert iers.conf.auto_download is auto_download
                assert iers.conf.auto_max_age == 30
                if cached:
                    assert pint.iers_b_file == cached_file
                    assert Path(cached_file).stat().st_mtime_ns == cached_mtime
                print('startup-ok')
        """
    )
    result = subprocess.run(
        [sys.executable, '-c', script, str(tmp_path), str(cached), str(auto_download)],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr or result.stdout
    assert 'startup-ok' in result.stdout
    assert ('bundled IERS-B' in result.stderr) is (not cached)
