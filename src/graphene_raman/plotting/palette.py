"""Palette and matplotlib style of the thesis figures (2026-09-21; in the package since 2026-09-30, ex scripts/_palette.py and
figures/memoire.mplstyle). Main colour: navy #000080 (= \\definecolor{darkblue}{rgb}{0,0,0.5}, the hyperlink colour of the PDF).
The cycler of memoire.mplstyle repeats CYCLE in the same order. One hue per category, never reused within a figure.
Usage: `from graphene_raman.plotting.palette import use_style, NAVY, ...; use_style()` (the style is found from __file__, never from the cwd)."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

STYLE = Path(__file__).with_name("memoire.mplstyle")   # single copy of the thesis style


def use_style():
    """Load the thesis style (memoire.mplstyle next to this file)."""
    plt.style.use(str(STYLE))


# metadata that would otherwise make every PDF/PNG differ at each run (creation date, software version): removed so that a
# regenerated figure is byte-identical to the tracked one (recette de la copie vers le dépôt du mémoire, 2026-09-30)
PDF_METADATA = {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}
PNG_METADATA = {"Software": None}


def save(fig, stem, outdir=None):
    """Save `fig` as <outdir>/<stem>.pdf (vector, for LaTeX) and .png (preview), reproducible byte for byte."""
    stem = Path(outdir) / stem if outdir is not None else Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), metadata=PDF_METADATA)
    fig.savefig(stem.with_suffix(".png"), metadata=PNG_METADATA)


NAVY = "#000080"      # principale : résultat de production / méthode validée
ORANGE = "#eb6834"    # contraste principal (Born, référence DFPT directe, Γ^ed dans le chapitre 5)
GREEN = "#1baf7a"
GOLD = "#eda100"
PINK = "#e87ba4"
SKY = "#56b4e9"       # « version claire » de la principale (même grandeur, paramètre non retenu)
REF = "#52514e"       # référence (DFT, grille grossière)
INK = "#0b0b0b"
MUTED = "#8a8984"
LIGHT = "#b5b4b0"
CYCLE = [NAVY, ORANGE, GREEN, GOLD, PINK, SKY]

# tailles de super-cellule : la référence de production (9×9) porte la principale
COL = {"9x9": NAVY, "7x7": ORANGE, "8x8": GREEN, "12x12": GOLD, "6x6": PINK, "5x5": SKY}

# cartes : séquentielle blanc → marine ; divergente marine ↔ blanc ↔ orange foncé
CMAP_SEQ = LinearSegmentedColormap.from_list("navy_seq", ["#ffffff", "#c6cfee", "#7f8fd8", "#3d4fb5", NAVY])
CMAP_DIV = LinearSegmentedColormap.from_list("navy_div", [NAVY, "#4f63c4", "#b9c2ea", "#ffffff", "#f5b899", "#e8783f", "#a83a0d"])
