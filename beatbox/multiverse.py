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
    # ====================================================================

    #Initialize the class variables
   
    truncated_nmax=None
    truncated_nmin=None
    truncated_lmax=None
    truncated_lmin=None
 
    # ====================================================================

    def __init__(self, truncated_nmax=None, truncated_nmin=1, truncated_lmax=None, truncated_lmin=1):
     
        Multiverse.truncated_nmax=truncated_nmax
        Multiverse.truncated_nmin=truncated_nmin
        Multiverse.truncated_lmax=truncated_lmax
        Multiverse.truncated_lmin=truncated_lmin
    
        self.all_data_universes=np.array([])
        self.all_simulated_universes=np.array([]) 

    # ====================================================================
    
    def __str__(self):
        return "a multiverse, capable of containing legions of Universes"

    # ====================================================================
        
    def create_original_Universe(self):
        '''
        Create a first universe in order to initiate the R matrix and
        the k filter for all subsequent universes. After this all 
        attibutes of the Universe class have their proper value.
        '''
        
        We_first=beatbox.Universe()
        We_first.set_Universe_k_filter()
        We_first.populate_Universe_R()
        
    def initiate_data_universe(self):
        '''
        Makes an instance of the universe, which is made of one Planck
        data realization.
        '''
        We=beatbox.Universe()
        
        self.all_data_universes = np.append(self.all_data_universes,We)

    def initiate_simulated_universe(self, truncated_nmax=None, truncated_nmin=None, truncated_lmax=None, truncated_lmin=None, n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=np.pi,Pvar=0.0):
        '''
        Makes an instance of the universe with containing a random
        realization of the gravitational field phi.
        '''
        We=beatbox.Universe()
        usedefault=1
        if truncated_nmax is not None:
            We.truncated_nmax=truncated_nmax
            usedefault=usedefault+1
        if truncated_nmin is not None:
            We.truncated_nmin=truncated_nmin
            usedefault=usedefault+1
        if truncated_lmax is not None:
            We.truncated_lmax=truncated_lmax   
            usedefault=usedefault+1
        if truncated_lmin is not None:
            We.truncated_lmin=truncated_lmin    
            usedefault=usedefault+1
        
            
        We.generate_a_random_potential_field(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin, n_s=n_s,kstar=kstar,PSnorm=PSnorm,Pdist=Pdist,Pmax=Pmax,Pvar=Pvar)
        We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)

        
        self.all_simulated_universes = np.append(self.all_simulated_universes,We)
    
    
    def read_Planck_samples(self):
        '''
        Read the 100 Planck samples into 100 instances of Universe 
        '''
        
        # download the tarball containing 100 posterior sample "COMMANDER-Ruler"
        #    low resolution maps, if not there already
        # tarball = "commander_32band_Clsamples100.tar.gz"
        datadir = "data/commander_32band_Clsamples100/"

        # if not os.path.isfile(tarball):
        #    URL = "http://folk.uio.no/ingunnkw/planck/32band/"+tarball
        #    !wget -O "$tarball" "$URL"
        #    !tar xvfz "$tarball"
        #    !mkdir -p "$datadir"
        #    !mv cmb_Cl_c000*.fits "$datadir"
        
        
        Tmapfiles = glob.glob(datadir+"cmb_Cl_c000*.fits")
        Nmaps = len(Tmapfiles) 
        
        self.all_data_universes = np.append(self.all_data_universes, [beatbox.Universe() for i in range(Nmaps)])
        
        
        for k in range(Nmaps):
            self.all_data_universes[-1-k].read_in_CMB_T_map(from_this=Tmapfiles[k])
            self.all_data_universes[-1-k].decompose_T_map_into_spherical_harmonics()
        print "Read in",Nmaps,"maps into",len(We),"beatbox universe objects."
        
        return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    