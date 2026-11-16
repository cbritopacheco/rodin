"""Shared figure style for CMAME (Elsevier) submissions.

Column widths follow the Elsevier artwork guide (single 90 mm, 1.5 column
140 mm, double 190 mm); text is 8 pt serif with STIX mathematics (Times-like,
as in the journal body). Set CMAME_USETEX=1 to typeset with a local LaTeX.
"""
import os
import shutil

import matplotlib as mpl

MM = 1.0 / 25.4
SINGLE = 90 * MM
ONEHALF = 140 * MM
DOUBLE = 190 * MM

PA_PER_MMHG = 133.322368
M3S_TO_MLMIN = 6.0e7

# Okabe-Ito (colour-blind safe); scenario -> colour, geometry -> dash.
SCENARIO_COLOR = {
    "healthy": "#000000",
    "aortic_stenosis": "#D55E00",
    "hyperviscosity": "#0072B2",
    "lvh": "#009E73",
}
SCENARIO_LABEL = {
    "healthy": "Healthy",
    "aortic_stenosis": "Aortic stenosis",
    "hyperviscosity": "Hyperviscosity",
    "lvh": "LV hypertrophy",
}
GEOMETRY_DASH = {"nostenosis": "-", "stenosis80": (0, (4.0, 1.6))}
GEOMETRY_LABEL = {"nostenosis": "No stenosis", "stenosis80": "80% stenosis"}
SCENARIO_ORDER = ["healthy", "aortic_stenosis", "hyperviscosity", "lvh"]
GEOMETRY_ORDER = ["nostenosis", "stenosis80"]


def color(scenario):
    return SCENARIO_COLOR.get(scenario, "#777777")


def scenario_label(scenario):
    return SCENARIO_LABEL.get(scenario, scenario.replace("_", " ").capitalize())


def geometry_label(geometry):
    return GEOMETRY_LABEL.get(geometry, geometry)


def dash(geometry):
    return GEOMETRY_DASH.get(geometry, ":")


def apply():
    usetex = os.environ.get("CMAME_USETEX", "0") == "1" and shutil.which("latex")
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "text.usetex": bool(usetex),
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon": False,
        "legend.handlelength": 2.4,
        "axes.grid": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "image.interpolation": "lanczos",
    })


def panel_label(ax, letter, x=-0.02, y=1.02):
    """(a), (b), ... in the upper-left corner, outside the axes."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, fontweight="normal")


def save(fig, path_without_ext):
    os.makedirs(os.path.dirname(path_without_ext) or ".", exist_ok=True)
    fig.savefig(path_without_ext + ".pdf")
    fig.savefig(path_without_ext + ".png", dpi=600)
