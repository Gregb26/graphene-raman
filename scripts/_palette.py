"""Palette des figures du mémoire (2026-09-21). Couleur principale : bleu marine #000080 (= \\definecolor{darkblue}{rgb}{0,0,0.5},
couleur des hyperliens du PDF). Importée par make_figures.py, make_figures_memoire.py, make_figures_epw.py ; le cycler de
figures/memoire.mplstyle reprend CYCLE dans le même ordre. Une teinte par catégorie, jamais recyclée dans une même figure."""
from matplotlib.colors import LinearSegmentedColormap

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
