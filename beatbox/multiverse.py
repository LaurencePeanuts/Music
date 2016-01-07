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

#import beatbox.universe
import beatbox

# ====================================================================

class Multiverse(object):
    """
    Universe sampler.
    """
    
    truncated_nmax=None
    truncated_nmin=None
    truncated_lmax=None
    truncated_lmin=None
    
    def __init__(self, truncated_nmax=None, truncated_nmin=1, truncated_lmax=None, truncated_lmin=1):
     
        beatbox.Multiverse.truncated_nmax=truncated_nmax
        beatbox.Multiverse.truncated_nmin=truncated_nmin
        Multiverse.truncated_lmax=truncated_lmax
        Multiverse.truncated_lmin=truncated_lmin
    
        self.all_data_universes=np.array([])
        self.all_simulated_universes=np.array([])
        
        #import beatbox.universe
        #from beatbox.universe import Universe   
        
    def create_original_Universe(self):
        #import beatbox.universe
        #from beatbox.universe import Universe
        We_first=beatbox.Universe()
        We_first.set_Universe_k_filter()
        We_first.populate_Universe_R()
        
    def initiate_data_universe(self):
        '''
        Makes an instance of the universe.
        '''
        We=beatbox.Universe()
        
        self.all_data_universes = np.append(self.all_data_universes,We)
    
    
    
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
    
    