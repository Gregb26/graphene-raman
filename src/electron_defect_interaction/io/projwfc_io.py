"""
projwfc_io.py
    Readers for projwfc.x outputs (R4, 2026-09-25):
        read_projwfc_states -- the atomic-wavefunction list ("state #") and, per k-point block and band, the
                               decomposition "psi = w*[# i] + ..." with |psi|^2, from projwfc.out
        read_pdos_m         -- m-resolved pdos_atm#N(X)_wfc#j(l) files (E, ldos, pdos per m and spin)
    projwfc prints only the components above its internal threshold (0.001), so the listed weights sum
    to less than |psi|^2. QE real-harmonic order for l = 1: m = 1 -> z, m = 2 -> x, m = 3 -> y.
"""
import re

import numpy as np

_STATE = re.compile(r"state #\s*(\d+): atom\s+(\d+)\s*\((\S+)\s*\), wfc\s+(\d+)\s*\(l=(\d+) m=\s*(\d+)\)")
_KLINE = re.compile(r"^\s*k =\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")
_ELINE = re.compile(r"==== e\(\s*(\d+)\) =\s*([-\d.]+) eV ====")
_COMP = re.compile(r"([\d.]+)\*\[#\s*(\d+)\]")
_PSI2 = re.compile(r"\|psi\|\^2 =\s*([\d.]+)")
M_LABEL = {(0, 1): "s", (1, 1): "pz", (1, 2): "px", (1, 3): "py"}


def read_projwfc_states(path):
    """
    Returns (states, bands):
        states: dict state_index (1-based) -> dict(atom (1-based), species, wfc, l, m, label)
        bands : list of dict(kblock (0-based index of the "k =" block; LSDA: 0 = up, 1 = down),
                             k (3,), band (1-based), e_eV, comps {state_index: weight}, psi2)
    """
    states = {}; bands = []; kblock = -1; k = None; cur = None
    with open(path) as f:
        for line in f:
            m = _STATE.search(line)
            if m:
                i, atom, sp, wfc, l, mm = m.groups()
                states[int(i)] = dict(atom=int(atom), species=sp, wfc=int(wfc), l=int(l), m=int(mm),
                                      label=M_LABEL.get((int(l), int(mm)), f"l{l}m{mm}"))
                continue
            m = _KLINE.match(line)
            if m:
                kblock += 1; k = np.array([float(x) for x in m.groups()]); cur = None
                continue
            m = _ELINE.search(line)
            if m:
                cur = dict(kblock=kblock, k=k, band=int(m.group(1)), e_eV=float(m.group(2)), comps={}, psi2=np.nan)
                bands.append(cur)
                continue
            if cur is not None:
                for w, i in _COMP.findall(line):
                    cur["comps"][int(i)] = cur["comps"].get(int(i), 0.0) + float(w)
                m = _PSI2.search(line)
                if m:
                    cur["psi2"] = float(m.group(1)); cur = None
    return states, bands


def read_pdos_m(path):
    """
    m-resolved pdos file. Returns (E (nE,), ldos (nE, nspin), pdos (nE, nspin, nm)).
    Column layout (QE): E, ldos[spin...], then for each m: pdos[spin...].
    """
    with open(path) as f:
        header = f.readline()
    nspin = 2 if "ldosdw" in header else 1
    data = np.loadtxt(path)
    E = data[:, 0]; ldos = data[:, 1:1 + nspin]
    rest = data[:, 1 + nspin:]
    nm = rest.shape[1] // nspin
    pdos = rest.reshape(len(E), nm, nspin).transpose(0, 2, 1)
    return E, ldos, pdos


def weights_by_group(band, states, groups):
    """
    Sum of the listed projwfc weights of one band over atom/orbital groups.
    groups: dict name -> (set of atom indices (1-based), set of labels among 's', 'pz', 'px', 'py').
    Returns dict name -> weight.
    """
    out = {g: 0.0 for g in groups}
    for i, w in band["comps"].items():
        st = states.get(i)
        if st is None:
            continue
        for g, (atoms, labels) in groups.items():
            if st["atom"] in atoms and st["label"] in labels:
                out[g] += w
    return out
