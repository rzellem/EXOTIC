"""Install the tested Colab wheels without compiling EXOTIC or UltraNest.

The default URL becomes available when the prepared repository files are
published. This installer only reads that location; it never uploads files.
Use --wheelhouse PATH to install from an existing local wheelhouse instead.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from urllib.request import urlopen, urlretrieve


WHEELHOUSE_URL = (
    'https://raw.githubusercontent.com/rzellem/EXOTIC/'
    'michael_fitzgeralds_wonderful_branch_of_magic/wheelhouse/colab'
)


def install_prebuilt(wheelhouse=WHEELHOUSE_URL):
    source = str(wheelhouse).rstrip('/')
    with TemporaryDirectory(prefix='exotic-colab-') as temporary:
        if source.startswith(('https://', 'http://')):
            with urlopen(source + '/index.json') as response:
                manifest = json.load(response)
            directory = Path(temporary)
            for wheel in manifest['wheels']:
                urlretrieve(source + '/' + wheel, directory / wheel)
        else:
            directory = Path(source).resolve()
            manifest = json.loads((directory / 'index.json').read_text(encoding='utf-8'))
        subprocess.check_call([
            sys.executable, '-m', 'pip', 'install', '--upgrade',
            '--only-binary=exotic,ultranest,exoticcolab',
            '--find-links', str(directory), *manifest['requirements'],
        ])
    print('Prebuilt EXOTIC, UltraNest and notebook display helper installed.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheelhouse', default=WHEELHOUSE_URL,
                        help='Wheelhouse directory or URL; default: the EXOTIC Colab wheelhouse')
    args = parser.parse_args()
    install_prebuilt(args.wheelhouse)


if __name__ == '__main__':
    main()
