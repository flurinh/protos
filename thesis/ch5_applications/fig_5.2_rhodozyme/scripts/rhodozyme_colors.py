"""
Rhodozyme figure color loader.

Loads rhodozyme_colors.yaml and exposes colors as a simple dict.
Both hex strings and rgb_norm tuples are available per key.

Usage:
    from rhodozyme_colors import C, hex, rgb

    ax.plot(..., color=hex("designed"))
    ax.scatter(..., color=rgb("theozyme"))
    # or directly:
    C["retinal"]["hex"]        # '#bc6c25'
    C["retinal"]["rgb_norm"]   # [0.737, 0.424, 0.145]
"""

from pathlib import Path

import yaml

_COLORS_FILE = Path(__file__).resolve().parent / "rhodozyme_colors.yaml"

with open(_COLORS_FILE) as _f:
    C = yaml.safe_load(_f)


def hex(key: str) -> str:
    """Return hex color string for a given role key."""
    return C[key]["hex"]


def rgb(key: str) -> tuple:
    """Return (r, g, b) normalized tuple for a given role key."""
    return tuple(C[key]["rgb_norm"])
