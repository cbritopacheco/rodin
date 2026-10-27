"""Authoritative dedicated lane definitions, separate from generated scheduling."""
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]


def jobs():
    return yaml.safe_load((ROOT / '.github/ci/targeted-lanes.yml').read_text())['jobs']
