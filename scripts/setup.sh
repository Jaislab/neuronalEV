#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ ! -x .venv/bin/python ]]; then
  python3 -m venv --system-site-packages .venv
fi

# Reuse the exact tested versions in a restored cloud snapshot, without downloads.
if ! .venv/bin/python - <<'PY'
from importlib import metadata
from pathlib import Path
import sys
for line in Path('requirements.lock').read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    name, version = line.split('==')
    try:
        installed = metadata.version(name)
    except metadata.PackageNotFoundError:
        sys.exit(1)
    if installed != version:
        sys.exit(1)
PY
then
  .venv/bin/python -m pip install -r requirements.lock
fi
.venv/bin/python -m pip install --no-index --no-deps --no-build-isolation -e .
.venv/bin/facs-neuro --help
