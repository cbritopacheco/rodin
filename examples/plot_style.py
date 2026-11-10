"""Apply the documentation's Manrope typography and TeX-style mathematics."""
from pathlib import Path
import os
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt


def apply_style():
    examples = Path(__file__).resolve().parent
    fonts = examples.parent / 'brand/fonts/static'
    for name in ['Manrope-Regular.ttf', 'Manrope-SemiBold.ttf', 'Manrope-Bold.ttf']:
        font_manager.fontManager.addfont(fonts / name)
    plt.style.use(examples / 'rodin.mplstyle')
    ink = '#dee3eb' if os.environ.get('RODIN_PLOT_THEME') == 'dark' else '#303b49'
    for key in ['text.color', 'axes.labelcolor', 'axes.titlecolor', 'xtick.color', 'ytick.color']:
        plt.rcParams[key] = ink
    plt.rcParams['image.cmap'] = 'turbo'
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=turbo_colors(4))


def turbo_colors(count):
    """Sample curve colors across the Turbo colormap."""
    if count < 1:
        raise ValueError('At least one curve color is required')
    palette = plt.get_cmap('turbo')
    return [palette(.12 + .76 * index / max(1, count - 1))
            for index in range(count)]
