import os.path
import numpy as np
import matplotlib
import matplotlib.pylab as plt
import healpy as hp
import string
import yt
import os
import glob
from PIL import Image as PIL_Image

import beatbox.universe

# ====================================================================

class Multiverse(object):
    """
    Universe sampler.
    """
    
    def __init__(self, truncated_nmax=None, truncated_nmin=1, truncated_lmax=None, truncated_lmin=1):
     
        self.truncated_nmax=truncated_nmax
        self.truncated_nmin=truncated_nmin
        self.truncated_lmax=truncated_lmax
        self.truncated_lmin=truncated_lmin
    
    def initiate_data_universe(self):
        '''
        Makes an instance of the universe.
        '''
        
        We=beatbox.Universe()
        if We.kfilter is None:
            We.set_Universe_k_filter()
        if We.R is None:
            We.populate_Universe_R()
        self.We=We
    
    
    
    def read_Planck_samples(self):
        
        datadir = "data/commander_32band_Clsamples100/"
        Tmapfiles = glob.glob(datadir+"cmb_Cl_c000*.fits")
        Nmaps = len(Tmapfiles)
        self.Tmap = np.array([],ndim=Nmapss) # We are now legion!

        for k in range(Nmaps):
            self.Tmap[k]=self.read_in_CMB_T_map(from_this=Tmapfiles[k])
        
            self.Tmapfile = from_this
            self.Tmap = hp.read_map(from_this)
            self.NSIDE = hp.npix2nside(len(self.Tmap))
        return
    
    