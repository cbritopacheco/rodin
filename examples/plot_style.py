"""Apply the documentation's Manrope typography and TeX-style mathematics."""
from pathlib import Path
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt


def apply_style():
    examples = Path(__file__).resolve().parent
    fonts = examples.parent / 'brand/fonts/static'
    for name in ['Manrope-Regular.ttf', 'Manrope-SemiBold.ttf', 'Manrope-Bold.ttf']:
        font_manager.fontManager.addfont(fonts / name)
    plt.style.use(examples / 'rodin.mplstyle')
