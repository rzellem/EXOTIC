"""Exercise the real startup dependency with an isolated Astropy cache."""
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest


@pytest.mark.parametrize('cached', [False, True])
@pytest.mark.parametrize('auto_download', [False, True])
def test_import_uses_local_iers_b_without_network(tmp_path, cached, auto_download):
    script = textwrap.dedent(
        """
        import logging
        import sys
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

            def no_download(*args, **kwargs):
                raise OSError('IERS publisher unavailable')

            with iers.conf.set_temp('auto_download', auto_download), \
                    iers.conf.set_temp('auto_max_age', 30), \
                    patch.object(data, '_download_file_from_source', no_download):
                import ldtk
                import exotic.exotic as reduction
                import barycorrpy.PINT_erfautils as pint

                assert callable(reduction.JDUTC_to_BJDTDB)
                assert len(pint.iers_b) > 0
                assert iers.conf.auto_download is auto_download
                assert iers.conf.auto_max_age == 30
                if cached:
                    assert pint.iers_b_file == cached_file
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
