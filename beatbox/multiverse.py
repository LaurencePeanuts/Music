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
        
        return

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
        
        return
    
    
    def initiate_data_universe(self):
        '''
        Makes an instance of the universe, which is made of one Planck
        data realization.
        '''
        We=beatbox.Universe()
        
        self.all_data_universes = np.append(self.all_data_universes,We)

        return 
    
    
    def initiate_simulated_universe(self, truncated_nmax=None, truncated_nmin=None, truncated_lmax=None, truncated_lmin=None, n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=np.pi,Pvar=0.0, fngrid=None):
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
        
        if fngrid is None: 
            We.generate_a_random_potential_field(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin, n_s=n_s,kstar=kstar,PSnorm=PSnorm,Pdist=Pdist,Pmax=Pmax,Pvar=Pvar)
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)

        else:
            We.fngrid=fngrid
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)
        
        
        
        self.all_simulated_universes = np.append(self.all_simulated_universes,We)
    
        return
    
    def read_Planck_samples(self):
        '''
        Read the 100 Planck posterior samples into 100 instances of Universe 
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
        print "Read in",Nmaps,"maps into",Nmaps,"beatbox universe objects."
        
        return
    
    
    def calculate_covariance_matrix(self, filename = 'covCyy.txt'):
        '''
        Calculate the a_y covariance matric from the
        100 Planck posterior samples, and save it to 
        disk
        '''
        Nmaps=100
        Planck_a_y=np.zeros((Nmaps, len(beatbox.Universe.lms)), dtype=np.float)
        for k in range(Nmaps):
            values=self.all_data_universes[-1-k].alm2ay()
            Planck_a_y[-1-k,:]=self.all_data_universes[-1-k].ay2ayreal_for_inference(values)
        
        r=0
        meanPlanck=np.zeros(len(beatbox.Universe.lms))
        meanPlanck=np.mean(Planck_a_y, axis=0)
        C_yy=np.zeros((len(self.all_data_universes[-1].lms), len(self.all_data_universes[-1].lms)))
        for r in range(Nmaps):
            C_yy=C_yy+np.outer(Planck_a_y[r,:]- meanPlanck, Planck_a_y[r,:]- meanPlanck)
        self.C_yy=C_yy/Nmaps
        
        np.savetxt( "data/"+filename, self.C_yy)
        
        return
    
    def load_covariance_matrix(self, filename = 'covCyy.txt'):
        '''
        Load the previously calculated a_y covariance matrix
        '''
        
        self.C_yy = np.loadtxt("data/"+filename)
        
    
    def calculate_sdv_Cyy_inverse(self):
        '''
        Calculate the inverse of the covariance matrix using singular
        value decomposition.
        '''
        
        U, s, V_star = np.linalg.svd(self.C_yy)
        
        #S = np.zeros((U.shape[1], V_star.shape[0]), dtype=complex)
        #S=np.diag(s)
        #S_cross = np.transpose(1./S)
        S_cross=np.diag(1./s)
        
        V = V_star.conj().T
        U_star=U.conj().T
        self.inv_Cyy=np.dot(V, np.dot(S_cross,U_star))
        
        return
    
    
    def reconstruct_3D_potential(self, datamap ,inv_Cyy=None):
        '''
        Given an observed (or mock) map of the CMB sky, reconstruct the most 
        probable value of the 3D interior sphere using inv_Cyy.
        '''
        
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy
        
        import scipy.stats as st
        
        naccept = 0
        niters = 10000
        self.samples = np.zeros(niters+1)
        self.samples = self.samples.tolist()
        self.logpost = np.zeros(niters+1)
        
        #self.initiate_simulated_universe()
        ay = self.all_simulated_universes[0].ay2ayreal_for_inference(self.all_simulated_universes[0].ay)
        fngrid = self.all_simulated_universes[0].fngrid
        
        self.samples[0] = self.all_simulated_universes[0].fngrid
        sigma = self.all_simulated_universes[0].Power_Spectrum
        
        self.inv_Cf=1./(self.all_simulated_universes[-1].Power_Spectrum)
        
        self.logpost[0] = self.get_logpost(datamap, fngrid, inv_Cyy)
        
        for i in range(niters):
            fngrid_p = fngrid + np.random.normal(0, np.sqrt(sigma) )
            
            logpost_p = self.get_logpost(datamap, fngrid_p, inv_Cyy)
            rho = min( 1, np.exp(logpost_p) / np.exp(self.logpost[-1]) )
            u = np.random.uniform()
            if u < rho:
                naccept += 1
                fngrid = fngrid_p
        
            self.samples[i+1] = fngrid
            self.logpost[i+1] = logpost_p
            
        nmcmc = len(self.samples)//2
        print "Efficiency = ", naccept/niters
         
        
        return
        
        
    def get_logpost(self, datamap, fngrid_p, inv_Cyy):
        
        self.initiate_simulated_universe(fngrid=fngrid_p)
        values = self.all_simulated_universes[-1].ay
        ay = self.all_simulated_universes[-1].ay2ayreal_for_inference(values)
        
        log_likelihood = -0.5 * np.dot( (datamap.T - ay.T), np.dot( inv_Cyy, (datamap-ay) ) )
        log_prior = -0.5 * np.sum( fngrid_p**2/self.inv_Cf )
        
        return log_likelihood+log_prior
    
    
    def solve_for_3D_potential(self, datamap, inv_Cyy=None):
        
        #Initiate the inverse of the a_y covariance matrix
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy
        
        #Initiate the inverse covariance matrix of the prior
        ind = np.where(beatbox.Universe.kfilter>0)
        PS = np.zeros(2*len(ind[1]))
        PS[:len(ind[1])] = (self.all_simulated_universes[-1].Power_Spectrum[ind]/2)
        PS[len(ind[1]):] = (self.all_simulated_universes[-1].Power_Spectrum[ind]/2)
        
        inv_Cf=np.diag(1./PS)
        
        #Deal with the R matrix to make it real
        # Select the m values out the the lms tupples
        m = np.array([m[1] for m in beatbox.Universe.lms])
        # Find the indices of the positive ms
        pos_ind = (m>0)
        # Find the indices of the m=0
        zero_ind = (m==0)
        # Find the indices of the negative ms
        neg_ind = (m<0)
        
        R_real = np.zeros((len(beatbox.Universe.lms), len(beatbox.You.all_simulated_universes[0].fn)), dtype=np.float)
        
        R_real[pos_ind,:] = beatbox.Universe.R[pos_ind,:].real
        R_real[neg_ind,:] = beatbox.Universe.R[pos_ind,:].imag
        R_real[zero_ind,:] = beatbox.Universe.R[zero_ind,:].astype(np.float)
        
        A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
        
        U, s, V_star = np.linalg.svd(A)
        inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
        
        from numpy.linalg import inv
        self.reconstrunct_fn = np.dot( inv_A , np.dot(R_real.T , np.dot (inv_Cyy , datamap) ) )
        
        return
    
    def generate_one_realization_of_noise(self):
        '''
        Generate one realization of the noise given by C_yy
        '''
        
        mean=np.zeros(self.C_yy.shape[0])
        noise = np.random.multivariate_normal(mean, self.C_yy, 1)
        
        return noise
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    