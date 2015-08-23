import numpy as np
import matplotlib.pylab as plt
import healpy as hp

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
        else:
            self.Tmap = hp.read_map(from_this)
        return

    def show_CMB_T_map(self,from_perspective_of="observer"):
        if from_perspective_of == "observer":
            # Sky map:
            hp.mollview(self.Tmap)
        else:
            #  "External" view ([like this](http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html))            pass
            # TO BE CODED!
        return

# ====================================================================
