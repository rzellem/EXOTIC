"""Prepare the checkout C++ hybrid4 backend (use --force to rebuild)."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from exotic.api.native_build import main


if __name__ == '__main__':
    main()
