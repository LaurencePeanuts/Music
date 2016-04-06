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
from scipy.stats import chi2

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
        self.all_reconstructed_universes=np.array([])
        
        try:
            #print beatbox.Universe.truncated_nmax
        
            i = 0
            reslms = 0
            if Multiverse.truncated_nmax is not beatbox.Universe.truncated_nmax:
                beatbox.Universe.truncated_nmax = Multiverse.truncated_nmax
                i=i+1
            if Multiverse.truncated_nmin is not beatbox.Universe.truncated_nmin:
                beatbox.Universe.truncated_nmin = Multiverse.truncated_nmin
                i=i+1
            if Multiverse.truncated_lmax is not beatbox.Universe.truncated_lmax:
                beatbox.Universe.truncated_lmax = Multiverse.truncated_lmax
                i=i+1
                reslms = reslms + 1
            if Multiverse.truncated_lmin is not beatbox.Universe.truncated_lmin:
                beatbox.Universe.truncated_lmin = Multiverse.truncated_lmin
                i=i+1
                reslms = reslms + 1
            if reslms is not 0:
                beatbox.Universe.lms = [(l, m) for l in range(beatbox.Universe.truncated_lmin,beatbox.Universe.truncated_lmax+1) for m in range(-l, l+1)]
            if i is not 0:
                We_first=beatbox.Universe()
                We_first.set_Universe_k_filter()
                We_first.populate_Universe_R()
                print 'Restarted new Universes from scratch with new properties. '
            print 'This is the end of the universe as we know it...'
        
        except AttributeError:
            i=None
        
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
    
    
    def initiate_simulated_universe(self, truncated_nmax=None, truncated_nmin=None, truncated_lmax=None, truncated_lmin=None, n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=np.pi,Pvar=0.0, fngrid=None, printout=1):
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
            We.generate_a_random_potential_field(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin, n_s=n_s,kstar=kstar,PSnorm=PSnorm, Pdist=Pdist, Pmax=Pmax,Pvar=Pvar, printout=printout)
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)

        else:
            We.fngrid=fngrid
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)
        
        
        
        self.all_simulated_universes = np.append(self.all_simulated_universes,We)
    
        return
    
    def read_Planck_samples(self, Nmaps=100, from_this = 'from_commander_32band_100', URL=None, fits_name=None):
        '''
        Read the Nmaps Planck posterior samples into Nmaps instances
        of Universe 
        - fits_name is a string containing the common name of the .fits
          files inside the tarball
        - URL is a string containing the path from which the tarball 
          must be downloaded
        - from_this is a string containing the name of the tarball 
          without the extension .tar.gz
        '''
        
        # If we are using the 100 Commander samples, load the Tmaps 
        if from_this is 'from_commander_32band_100':
            # download the tarball containing 100 posterior sample 
            #    "COMMANDER-Ruler" low resolution maps, if not there already
            tarball = "commander_32band_Clsamples100.tar.gz"
            datadir = "../data/commander_32band_Clsamples100/"
            
#            test1 = "testdirectory"
#            testcomd = "mkdir " + test1
#            print testcomd
#            os.system(testcomd)
            
            if not os.path.isfile(datadir+tarball):
                URL = "http://folk.uio.no/ingunnkw/planck/32band/"+tarball
                command = "wget -O " + tarball + " " +URL
                os.system(command)
                tarcmd = "tar xvfz " + tarball
                os.system(tarcmd)
                mkdircmd = "mkdir -p " + datadir
                os.system(mkdircmd)
                mvcmd = "cmb_Cl_c000*.fits" + datadir
                os.system(mvcmd)
        
        
            Tmapfiles = glob.glob(datadir+"cmb_Cl_c000*.fits")
            Nmaps = len(Tmapfiles) 
        
        else:
            datadir = "../data/%s" % from_this
            tarball = "%s.tar.gz" % from_this
            
            if not os.path.isfile(datadir+tarball):
                URL = URL+tarball
                command = "wget -O " + tarball + " " +URL
                os.system(command)
                tarcmd = "tar xvfz " + tarball
                os.system(tarcmd)
                mkdircmd = "mkdir -p " + datadir
                os.system(mkdircmd)
                mvcmd = fits_name+"*.fits" + datadir
                os.system(mvcmd)
               
        
            Tmapfiles = glob.glob(datadir+fits_name+"*.fits")
            Nmaps = len(Tmapfiles) 
        
        self.all_data_universes = np.append(self.all_data_universes, [beatbox.Universe() for i in range(Nmaps)])
        
        
        for k in range(Nmaps):
            self.all_data_universes[-1-k].read_in_CMB_T_map(from_this=Tmapfiles[k])
            self.all_data_universes[-1-k].decompose_T_map_into_spherical_harmonics()
        print "Read in",Nmaps,"maps into",Nmaps,"beatbox universe objects."
        
        return
    
    
    def calculate_covariance_matrix(self, filename = 'lmax8'):
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
        
        np.savetxt( "../data/covCyy_"+filename+".txt", self.C_yy)
        
        return
    
    def load_covariance_matrix(self, filename = 'covCyy.txt'):
        '''
        Load the previously calculated a_y covariance matrix
        '''
        
        self.C_yy = np.loadtxt("../data/"+filename)
        
    
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
        
        
        
        naccept = 0
        niters = 10000
        self.samples = np.zeros(niters+1)
        self.samples = self.samples.tolist()
        self.logpost = np.zeros(niters+1)
        
        
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
    
    
    def solve_for_3D_potential(self, datamap, inv_Cyy=None, print_alpha=0):
        
        #Initiate the inverse of the a_y covariance matrix
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy
        
        #Initiate the inverse covariance matrix of the prior
        ind = np.where(beatbox.Universe.kfilter>0)
        PS = np.zeros(2*len(ind[1]))
        PS[:len(ind[1])] = (self.all_simulated_universes[-1].Power_Spectrum[ind])
        PS[len(ind[1]):] = (self.all_simulated_universes[-1].Power_Spectrum[ind])
        
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
        
        self.R_real = R_real
        
        # Define the A matrix and its inverse 
        A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
        
        U, s, V_star = np.linalg.svd(A)
        inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
        
        #Solve for the normalization of the prior
        inv_Cf = self.solve_for_prior_normalization(inv_Cyy, inv_Cf, A, inv_A, R_real, datamap, print_alpha)
        
        # Redefine A with the properly normalized prior
        A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
        
        # Find the inverse of A that has been the properly normalized Cf
        U, s, V_star = np.linalg.svd(A)
        inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
        self.inv_A=inv_A
        
        # Use linear algebra to solve the A*f_n=b linear equation
        b =  np.dot(R_real.T , np.dot (inv_Cyy , datamap) )
        #self.reconstrunct_fn = np.linalg.solve(A, b)
    
        #from numpy.linalg import inv
        self.reconstrunct_fn = np.dot( inv_A , np.dot(R_real.T , np.dot (inv_Cyy , datamap) ) )
        
        return
    
    def solve_for_prior_normalization(self, inv_Cyy, inv_Cf, A, inv_A, R_real, datamap, print_alpha):
        
        N = R_real.shape[1]
        
        #from numpy.linalg import inv
        alpha = N / ( np.trace( np.dot(inv_A, inv_Cf) ) + np.dot( datamap.T  , np.dot( inv_Cyy, np.dot(R_real, np.dot( inv_A , np.dot(inv_Cf , np.dot( inv_A , np.dot( R_real.T ,np.dot(inv_Cyy, datamap)))))))) )
        if print_alpha is 1:
            outfile = 'alpha_lmax' + str(Multiverse.truncated_lmax) + 'lmin' + str(Multiverse.truncated_lmin) + 'nmax' + str(Multiverse.truncated_nmax) + 'nmin' + str(Multiverse.truncated_nmin) + '.txt'
            f = open('RobustnessAnalysis/' + outfile, 'a')
            towrite = str(alpha) + "\n"
            f.write(towrite)
            f.close()
            
        return alpha*inv_Cf
    
    def generate_one_realization_of_noise(self):
        '''
        Generate one realization of the noise given by C_yy
        '''
        
        #np.random.seed(3)
        
        mean=np.zeros(self.C_yy.shape[0])
        noise = np.random.multivariate_normal(mean, self.C_yy, 1)
        
        return noise
    
    
    def generate_realizations_from_posterior(self, mean, number_of_realizations=100):
        '''
        Generate realizations of the posterior
        '''
        
        if self.inv_A is None:
            print 'A matrix not initialized'
            return
            
        samples = np.random.multivariate_normal(mean.reshape(len(mean)), self.inv_A, number_of_realizations)
        
        
        self.all_simulated_universes = np.append(self.all_simulated_universes, [beatbox.Universe() for i in range(number_of_realizations)])
        
        
        
        for k in range(number_of_realizations):
            self.all_simulated_universes[-1-k].fn = samples[k,:]
            self.all_simulated_universes[-1-k].transform_3D_potential_into_alm(usedefault=1, fn=1)
        
        return 
    
    
    
    def calculate_chi2_in_posterior(self, fn_true, fn_rec):
        '''
        Calculate the chi2 and the p-value distance between the most probable
        values of fn's and their true, noiseless value.
        '''
        
        if self.inv_A is None:
            print 'A matrix not initialized'
            return
        
        Delta_fn = fn_true-fn_rec
        chi2value = np.dot (Delta_fn.T , np.dot( self.inv_A, Delta_fn  ))
        
        p_value = 1-chi2.cdf(chi2value, len(fn_true)) 
        
        return p_value, chi2value
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    