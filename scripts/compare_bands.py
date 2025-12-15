"""
compare_bands.py
    Python module to compare the DFT electronic bands computed by Abinit and the Wannier interpolated bands computed by Wannier90.
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def read_eig_nc(filepath):
    """
    Reads and extracts the eigenvalues written in a Abinit *EIG.nc output file.
    Inputs:
        filepath: str, path to EIG.nc file containing eigenvalues
    Returns:
        dft_eig: (nspol, nkpt, nband) array of floats: eigenvalues in eV computed by Abinit DFT code. nspol in number of spin polarizations, nkpt is number of kpoints, band is number of DFT bands
        kpt: (nkpt, 3) array of floats, kpoint path in reduced coordinates
        fermie: float, Fermie energy in eV
    """
    with Dataset(filepath, mode='r') as nc:
        dft_eig = np.array(nc.variables["Eigenvalues"][:]) # (nspol, nkpt, nband), eigenvalues in Ha on kpoint path
        kpt = np.array(nc.variables["Kptns"][:]) # (nkpt, 3), kpoint path in reduced coordinates on which the band structure was computed
        fermie = np.array(nc.variables["fermie"][:]) # float, Fermi energy in Ha

        Ha_to_eV = 27.211386245988
        # convert to eV and shift by Fermi energy
        dft_eig = Ha_to_eV * (dft_eig - fermie)
        fermie *= Ha_to_eV

        return dft_eig, kpt, fermie

def check_hermicity_HR(HR, R):
    """
    Checks the Hermicity of the Hamiltonian written in the Wannier tight binding basis. The condition is
    H(R) = H^\\dag(-R).
    Inputs
        HR: (nrpts, nw, nw) array of complex
            Hamiltonien written in Wannier tight binding basis
        R: (nrpts, 3) array of floats
            R vectors computed by Wannier90 used in the tight binding basis
    """

    # Map R = (n1, n2, n3) to index in array R 
    R_tuples = [tuple(r) for r in R]
    R_index  = {r: i for i, r in enumerate(R_tuples)}

    for i, Rs in enumerate(R):
        R_ = tuple(-Rs) # -R

        # some R may be unpaired, skip those
        j = R_index.get(R_)
        if j is None:
            continue
        H = HR[i]
        H_h = HR[j].conj().T

        if not np.allclose(H, H_h, atol=1e-10):
            print(f"Hermiticity fails for R = {R_}, -R index = {j}")
            raise TypeError('Hamiltonian not Hermitian')
    
    return True

def read_tb_dat(w90_path):
    """
    Reads H(R) and R from Wannier90 tb.dat output files.
    Inputs:
        w90_path: str, path to tb.dat Wannier90 output file
    Returns:
        H(R): (nR, nwann, nwann) array of complex, Hamilontian in Wannier basis computed by Wannier90
        R: (nR, 3) array of ints: lattice vectors used by Wannier90 to construct H(R)
    """
    from pathlib import Path
    with Path(w90_path).open() as f:

        date = f.readline()

        # Lattice vectors in angstrom
        a1 = np.fromstring(f.readline(), sep=' ')
        a2 = np.fromstring(f.readline(), sep=' ')
        a3 = np.fromstring(f.readline(), sep=' ')

        A = np.column_stack((a1, a2, a3)) # A[:,i] = a_i

        nwann = int(f.readline()) # number of wannier functions

        nrpts = int(f.readline()) # number of inequivalent R points used by Wannier90 

        # skip the R point degeneracy list
        count = 0
        while count < nrpts:
            line = f.readline().split()
            count += len(line)

        HR = []
        R = []
        for _ in range(nrpts):

            f.readline() # blank space

            r = np.array([int(n) for n in f.readline().split()]) # R vectors
            hr = np.zeros((nwann, nwann), dtype=complex) # allocate matrix

            for _ in range(nwann*nwann):
                m_str, n_str, Re_str, Im_str = f.readline().split()
                m = int(m_str) - 1 # convert to 0 based
                n = int(n_str) - 1
                hr[m,n] = float(Re_str) + 1j*float(Im_str)
            
            R.append(r)
            HR.append(hr)

        HR = np.array(HR)
        R = np.array(R) 

        if check_hermicity_HR(HR,R):
            return HR, R

def main():
    dft_eig_path = "EIG.nc"
    wannier_path = "tb.dat"

    # get dft eigenvalues and kpoint path used
    dft_eig, kpt, fermie = read_eig_nc(dft_eig_path)

    # get H(R) and R
    HR, R = read_tb_dat(wannier_path) # (nR, nwann, nwann) and (nR, 3)

    # Interpolate HR on kpoint path using a Fourier transform
    phase = np.exp(2j*np.pi * (kpt @ R.T)) # (nkpt, nR)
    Hk = np.einsum("kr, rmn -> kmn", phase, HR) # (nkpt, nwann, nwann), summed over R

    # Diagonalize Hk to get Wannier eigenvalues on kpoint path
    wannier_eig, _ = np.linalg.eigh(Hk) # (nkpt, nwann)
    wannier_eig -= fermie

    # plotting
    nspol, nkpt, nband = dft_eig.shape
    _, nwann = wannier_eig.shape

    if nspol==1:
        fig, ax = plt.subplots(figsize=(5,6), dpi=120)
        for ib in range(nband):
            if ib==0:
                ax.plot(dft_eig[0,:,ib], "k", lw=1.0, label='Abinit')
            else:
                ax.plot(dft_eig[0,:,ib], "k", lw=1.0)
        for iw in range(nwann):
            if iw==0:
                ax.plot(wannier_eig[:,iw], "--r", lw=1.0, label="Wannier90")
            else:
                ax.plot(wannier_eig[:,iw], "--r", lw=1.0)
        ax.set_xlabel("k-path", fontsize=12)
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_title("Band Structure Comparison", fontsize=13)
        ax.grid(True, alpha=0.3)
        ax.tick_params(direction='in', length=4)
        ax.set_xlim(0,nkpt)
        ax.set_ylim(np.min(dft_eig), np.max(dft_eig))
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig("comparison.png")
        plt.show()

    else:
        raise ValueError(f"Code not written for spin polarization greater than one. You have nspol={nspol}")

if __name__ == "__main__":
    main()