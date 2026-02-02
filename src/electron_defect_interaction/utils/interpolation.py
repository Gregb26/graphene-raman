"""
interpolation.py

Python module containing functions to do interpolation of all sorts.
"""

import numpy as np

def trilinear_periodic(F,x,y,z,dx=1.0,dy=1.0,dz=1.0):
    """
    Periodc trilinear interpolation on a uniform 3D grid.
    Inputs:
        F: (Nx,Ny,Nz) array of complex, function to interpolate.
        x,y,z: floats, coordinates of the query points.
        dx,dy,dz: floats, grid spacing.
    Returns:
        Interpolated values F(x,y,z).
    """
    F = np.asarray(F)
    Nx, Ny, Nz = F.shape

    # fractional indices
    u = np.asarray(x) / dx; v = np.asarray(y) / dy; w = np.asarray(z) / dz

    # wrap periodically to [0,Ni)
    u = np.mod(u, Nx); v = np.mod(v, Ny); w = np.mod(w, Nz)

    # grid point to the left of query point
    i0 = np.floor(u, dtype=int); j0 = np.floor(v, dtype=int); k0 = np.floor(w, dtype=int)
    # grid point to the right of query point
    i1 = (i0 + 1) % Nx; j1 = (j0 + 1) % Ny; k1 = (k0 + 1) % Ny;  

    # trinlinear interpolation
    t = u - i0; s = v - j0; r = w - k0
    # gather 8 corners
    c000 = F[i0, j0, k0]
    c100 = F[i1, j0, k0]
    c010 = F[i0, j1, k0]
    c110 = F[i1, j1, k0]
    c001 = F[i0, j0, k1]
    c101 = F[i1, j0, k1]
    c011 = F[i0, j1, k1]
    c111 = F[i1, j1, k1]
    # trilinear blend
    one_t = 1.0 - t
    one_s = 1.0 - s
    one_r = 1.0 - r

    return (
        c000 * one_t * one_s * one_r +
        c100 * t     * one_s * one_r +
        c010 * one_t * s     * one_r +
        c110 * t     * s     * one_r +
        c001 * one_t * one_s * r     +
        c101 * t     * one_s * r     +
        c011 * one_t * s     * r     +
        c111 * t     * s     * r
    )


 