import numpy as np
import matplotlib.pylab as plt
import healpy as hp

import beatbox

# ====================================================================

class Universe(object):
    """
    A simple model universe in a box.
    """
    def __init__(self):
        self.size = 100
        self.phi = np.zeros([3,self.size])
        return

    def __str__(self):
        return "a model universe, containing potential map phi"

    def read_in_CMB_T_map(self,from_this=None):
        if from_this is None:
            print "No CMB T map file supplied."
            self.Tmap = None
            self.NSIDE = None
        else:
            self.Tmap = hp.read_map(from_this)
            self.NSIDE = hp.npix2nside(len(self.Tmap))
        return

    def show_CMB_T_map(self,from_perspective_of="observer"):
        if from_perspective_of == "observer":
            # Sky map:
            hp.mollview(self.Tmap)
        else:
            # "External" view ([like this](http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html))            pass
            #   beatbox.zoncaview(self.Tmap)
            # This did not work, sadly. Maybe we can find a 3D
            # spherical surface plot routine using matplotlib?
            print "Ahem - we can't visualize maps on the surface of the sphere yet, sorry."
            pass
        return

    def decompose_T_map_into_spherical_harmonics(self,lmax=None):
        """
        See healpy documentation at https://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.map2alm.html
        self.alm is a 1D numpy array of type=complex128.
        Indexing is described at https://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.Alm.html
        """
        if lmax is None:
            self.lmax = 3*self.NSIDE - 1
        else:
            self.lmax = lmax
        self.mmax = self.lmax

        self.alm = hp.sphtfunc.map2alm(self.Tmap,lmax=self.lmax,mmax=self.mmax)

        return

    def show_one_spherical_harmonic_of_CMB_T_map(self,l=1,m=1,max=20):
        i = hp.Alm.getidx(self.lmax, l, m)
        projected_alm = self.alm * 0.0
        projected_alm[i] = self.alm[i]
        projected_map = hp.alm2map(projected_alm,self.NSIDE)
        hp.mollview(projected_map,min=-max,max=max)
        return

    def show_lowest_spherical_harmonics_of_CMB_T_map(self,lmax=3,max=20):
        i = []
        for l in range(lmax):
            for m in range(-lmax,lmax+1):
                i.append(hp.Alm.getidx(self.lmax, l, m))
        print "Displaying sky map of the ",len(i)," lowest spherical harmonics only..."
        truncated_alm = self.alm * 0.0
        truncated_alm[i] = self.alm[i]
        truncated_map = hp.alm2map(truncated_alm,self.NSIDE)
        hp.mollview(truncated_map,min=-max,max=max)
        return

# ====================================================================
