"""
Adapted from a code snippet for making an interactive 3D visualization
of a 2D map on a sphere written by Andrea Zonca (https://github.com/zonca).
See his demo here:

http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html

The original code snippet can be found here:

https://gist.github.com/zonca/5146356#file-3d_map-py

I had to install some new modules to make this work:

    conda install vtk
    conda install mayavi

I then added a `mlab.draw()` statement to try and make the plot render,
but no joy - seems as though mayavi/vtk and the ipython notebook just don't
play well together. Seems like others are having this problem too.
Some hope here: http://nbviewer.ipython.org/github/empet/Math/blob/master/DomainColoring.ipynb
but otherwise, this is staying in the attic.
"""

from mayavi import mlab
import numpy as np
import healpy as hp

def zoncaview(m):
    """
    m is a healpix sky map, such as provided by WMAP or Planck.
    """

    nside = hp.npix2nside(len(m))
    vmin = -1e3; vmax = 1e3

    # Set up some grids:
    xsize = ysize = 1000
    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    # Project the map to a rectangular matrix xsize x ysize:
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)
    grid_map = m[grid_pix]

    # Create a sphere:
    r = 0.3
    x = r*np.sin(THETA)*np.cos(PHI)
    y = r*np.sin(THETA)*np.sin(PHI)
    z = r*np.cos(THETA)

    # The figure:
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    mlab.clf()

    mlab.mesh(x, y, z, scalars=grid_map, colormap="jet", vmin=vmin, vmax=vmax)

    mlab.draw()

    return
