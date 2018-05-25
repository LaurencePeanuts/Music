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

import time

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
    
    
    def initiate_simulated_universe(self, truncated_nmax=None, truncated_nmin=None, truncated_lmax=None, truncated_lmin=None, n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=2*np.pi,Pvar=0.0, fngrid=None, printout=1,dofft=0):
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
            We.generate_a_random_potential_field(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin, n_s=n_s,kstar=kstar,PSnorm=PSnorm, Pdist=Pdist, Pmax=Pmax,Pvar=Pvar, printout=printout, do_fft=dofft)
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault)

        else:
            We.fngrid=fngrid
            We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=usedefault, do_fft=dofft)
        
        
        
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
        
        #U, s, V_star = np.linalg.svd(self.C_yy)
        
        #S = np.zeros((U.shape[1], V_star.shape[0]), dtype=complex)
        #S=np.diag(s)
        #S_cross = np.transpose(1./S)
        #S_cross=np.diag(1./s)
        
        #V = V_star.conj().T
        #U_star=U.conj().T
        #self.inv_Cyy=np.dot(V, np.dot(S_cross,U_star))
        
        self.inv_Cyy = np.linalg.inv(self.C_yy)
        
        return

    def calculate_A_matrix(self, inv_Cyy = None, inv_Cf = None, R_real = None, filename = 'nmax_unspecified'):
        '''
        Calculate the A matrix, and save it to disk
        '''
        #Initiate the inverse of the a_y covariance matrix
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy

        if inv_Cf is None:
#            #Initiate the inverse covariance matrix of the prior
#            ind = np.where(beatbox.Universe.kfilter>0)
#            PS_long = np.zeros(2*len(ind[1]))
#            PS_long[:len(ind[1])] = (self.all_simulated_universes[-1].Power_Spectrum[ind])
#            PS_long[len(ind[1]):] = (self.all_simulated_universes[-1].Power_Spectrum[ind])
#        
#            PS = np.zeros(len(ind[1]))
#            PS[:len(ind[1])/2] = PS_long[:len(ind[1])/2]
#            PS[len(ind[1])/2:] = PS_long[len(ind[1]):3*len(ind[1])/2]
#        
#            inv_Cf=np.diag(1./PS)
            print 'Prior matrix not provided, using default Universe'
            
            ind = np.where(beatbox.Universe.kfilter>0)
            fn_length = len(ind[1])
            self.inv_Cf = calculate_prior_Cf(ind, fn_length)

        if R_real is None:
            print 'R_real matrix not provided, using default Universe'
            #Deal with the R matrix to make it real
            # Select the m values out the the lms tupples
            m = np.array([m[1] for m in beatbox.Universe.lms])
            # Find the indices of the positive ms
            pos_ind = (m>0)
            # Find the indices of the m=0
            zero_ind = (m==0)
            # Find the indices of the negative ms
            neg_ind = (m<0)
        
            R_real = np.zeros((len(beatbox.Universe.lms), len(beatbox.You.all_simulated_universes[-1].fn)), dtype=np.float)
        
            R_real[pos_ind,:] = beatbox.Universe.R[pos_ind,:].real
            R_real[neg_ind,:] = beatbox.Universe.R[pos_ind,:].imag
            R_real[zero_ind,:] = beatbox.Universe.R[zero_ind,:].astype(np.float)
        
            self.R_real = R_real
        
        A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
        self.A=A
        
        inv_A = self.calculate_pure_A_inverse()
        
#        np.savetxt( "../data/A_"+filename+".txt", A)
        np.save("../data/A_"+filename+".npy", A)
        np.save("../data/inv_A_"+filename+".npy", inv_A)
        return inv_A
    
    def load_A_matrix(self,  filename, inv = 'n'):
        '''
        Load the previously calculated a_y covariance matrix
        '''
        
#        self.A = np.loadtxt("../data/"+filename)
        if inv == 'n':
            self.A = np.load("../data/A_"+filename+".npy")
        else:
            self.inv_A = np.load("../data/inv_A_"+filename+".npy")
        return
    
    def calculate_sdv_A_inverse(self):
        '''
        Calculate the inverse of the covariance matrix using singular
        value decomposition.
        '''
        
        U, s, V_star = np.linalg.svd(self.A)
        
        #S = np.zeros((U.shape[1], V_star.shape[0]), dtype=complex)
        #S=np.diag(s)
        #S_cross = np.transpose(1./S)
        S_cross=np.diag(1./s)
        
        V = V_star.conj().T
        U_star=U.conj().T
        self.inv_A=np.dot(V, np.dot(S_cross,U_star))
        
        return inv_A
    
    def calculate_pure_A_inverse(self):
        '''
        Calculate the inverse of the covariance matrix using singular
        value decomposition.
        '''
        
        from numpy.linalg import inv
        
        inv_A=inv(self.A)
        
        return inv_A
    
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
    
    

    def calculate_prior_Cf(self, ind, fn_length, n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=2*np.pi,Pvar=0.0):
        '''
        Initiate the inverse covariance matrix of the prior.
        '''
        
        # Define the constants that go in the power spectrum
        # Change units of the pivot scale kstar from Mpc^-1 to normalize the smallest k
        #    mode to 1 (i.e. the radius of the CMB photosphere at 13.94Gpc)
        kstar = kstar*1.394e4

        # Draw Gaussian random Fourier coefficients with a k^{-3+(n_s-1)} power spectrum:
        Power_Spectrum = PSnorm*10000*np.power((beatbox.Universe.k/kstar) ,(-3+(n_s-1)))
        Power_Spectrum[np.isinf(Power_Spectrum)] = 10**-9

        
        #ind = np.where(beatbox.Universe.kfilter>0)
        PS_long = np.zeros(2*fn_length)
        PS_long[:fn_length] = (self.all_simulated_universes[-1].Power_Spectrum[ind]) 
        PS_long[fn_length:] = (self.all_simulated_universes[-1].Power_Spectrum[ind]) 
        
        PS = np.zeros(fn_length)
        PS[:fn_length/2] = PS_long[:fn_length/2]
        PS[fn_length/2:] = PS_long[fn_length:3*fn_length/2]
        
        inv_Cf = np.diag(1./PS) * 2 # factor of 2 is to account for real and imag covariance from power-spectrum covariance
            
        self.inv_Cf = inv_Cf
        
        return inv_Cf
        
    def solve_for_3D_potential(self, datamap, inv_Cyy=None, A=None, print_alpha=0):
        
        start = time.time()
        
        #Initiate the inverse of the a_y covariance matrix
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy
        
#        end1 = time.time()
#        print end1-start1
        
        
        
#        end2 = time.time()
#        print end2-end1
        
        #Deal with the R matrix to make it real
        # Select the m values out the the lms tupples
        m = np.array([m[1] for m in beatbox.Universe.lms])
        # Find the indices of the positive ms
        pos_ind = (m>0)
        # Find the indices of the m=0
        zero_ind = (m==0)
        # Find the indices of the negative ms
        neg_ind = (m<0)
        
        #Calculate the number of fns for the recontruction
        ind = np.where(beatbox.Universe.kfilter>0)
        fn_length = len(ind[1])
        
        
        R_real = np.zeros((len(beatbox.Universe.lms), fn_length), dtype=np.float)
        
        R_real[pos_ind,:] = beatbox.Universe.R[pos_ind,:].real
        R_real[neg_ind,:] = beatbox.Universe.R[pos_ind,:].imag
        R_real[zero_ind,:] = beatbox.Universe.R[zero_ind,:].astype(np.float)
        
        self.R_real = R_real
        
        end2 = time.time()
        print end2-start
#        print 'begining of long part'
        
        # Define the A matrix and its inverse 
#        A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
#        self.A1= A
        
#        end4 = time.time()
#        print 'A'
#        print end4-end3
        
#        U, s, V_star = np.linalg.svd(A)
        
#        end5 = time.time()
#        print 'B'
#        print end5-end4
        
#        inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
        
        
        #Solve for the normalization of the prior
        #inv_Cf = self.solve_for_prior_normalization(inv_Cyy, inv_Cf, A, inv_A, R_real, datamap, print_alpha)
        
        # Redefine A with the properly normalized prior
        
        if A is None:
            
            inv_Cf = self.calculate_prior_Cf(ind=ind, fn_length=fn_length)
            
#            #Initiate the inverse covariance matrix of the prior
#            ind = np.where(beatbox.Universe.kfilter>0)
#            PS_long = np.zeros(2*len(ind[1]))
#            PS_long[:len(ind[1])] = (self.all_simulated_universes[-1].Power_Spectrum[ind]) 
#            PS_long[len(ind[1]):] = (self.all_simulated_universes[-1].Power_Spectrum[ind]) 
#        
#            PS = np.zeros(len(ind[1]))
#            PS[:len(ind[1])/2] = PS_long[:len(ind[1])/2]
#            PS[len(ind[1])/2:] = PS_long[len(ind[1]):3*len(ind[1])/2]
#        
#            inv_Cf = np.diag(1./PS) * 2 # factor of 2 is to account for real and imag covariance from power-spectrum covariance
#            
#            self.inv_Cf = inv_Cf
            
            filename = 'nmax'+str(beatbox.Universe.truncated_nmax)
            
            inv_A = self.calculate_A_matrix(inv_Cyy = inv_Cyy, R_real = R_real, inv_Cf = inv_Cf, filename = filename)
            #A = self.A
        else:
            filename = 'nmax'+str(beatbox.Universe.truncated_nmax)
            
            end4 = time.time()

            self.load_A_matrix(filename, inv = 'y')
            
            end5 = time.time()
            print "time for loading is:"
            print end5-end4

            inv_A = self.inv_A
            print "successfully loaded A from the disk."
        #self.A = A
        
        
        # Find the inverse of A that has been the properly normalized Cf
#        U, s, V_star = np.linalg.svd(A)
#        inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
#        self.inv_A=inv_A
        

        
#        U2, s2, V_star2 = np.linalg.svd(inv_A)
#        inv_inv_A = np.dot(V_star2.T, np.dot(np.diag(1./s2),U2.T))
#        self.A=inv_inv_A
        

        
        
        # Use linear algebra to solve the A*f_n=b linear equation
        b =  np.dot(R_real.T , np.dot (inv_Cyy , datamap) )
        #self.reconstrunct_fn = np.linalg.solve(A, b)
    
        #from numpy.linalg import inv
#        self.reconstrunct_fn = np.linalg.solve(A, b)
        self.reconstrunct_fn = np.dot( inv_A , np.dot(R_real.T , np.dot (inv_Cyy , datamap) ) )
        
        end = time.time()
        print 'total time is:'
        print end-start
        
        return
    
    def solve_for_prior_normalization(self, inv_Cyy, inv_Cf, A, inv_A, R_real, datamap, print_alpha):
        
        N = R_real.shape[1]
        
        #from numpy.linalg import inv
        alpha = N / ( np.trace( np.dot(inv_A, inv_Cf) ) + np.dot( datamap.T  , np.dot( inv_Cyy, np.dot(R_real, np.dot( inv_A , np.dot(inv_Cf , np.dot( inv_A , np.dot( R_real.T ,np.dot(inv_Cyy, datamap)))))))) )
        
        if print_alpha is 1:
            alpha_towrite=np.flatten(alpha)[0]
            outfile = 'alpha_lmax' + str(Multiverse.truncated_lmax) + 'lmin' + str(Multiverse.truncated_lmin) + 'nmax' + str(Multiverse.truncated_nmax) + 'nmin' + str(Multiverse.truncated_nmin) + '.txt'
            f = open('RobustnessAnalysis/' + outfile, 'a')
            towrite = str(alpha_towrite) + "\n"
            f.write(towrite)
            f.close()
            
        return alpha*inv_Cf
    
    def solve_for_3D_potential_marginalized_over_large_n(self, datamap, inv_Cyy=None, print_alpha=0, truncated_nmax=None, nmax=None, NewAcov = None, inv_NewAcov = None, count=0):
        
        #Initiate the inverse of the a_y covariance matrix
        if inv_Cyy is None:
            inv_Cyy = self.inv_Cyy
        if truncated_nmax is None:
            truncated_nmax = beatbox.Universe.truncated_nmax
        if nmax is None:
            nmax = beatbox.Universe.truncated_nmax
        
        #Get the indices of the f_n's ordered by k norm
        ind = np.where(beatbox.Universe.kfilter>0)
        
        k, theta, phi = beatbox.Universe.k[ind], np.arctan2(beatbox.Universe.ky[ind],beatbox.Universe.kx[ind]), np.arccos(beatbox.Universe.kz[ind]/beatbox.Universe.k[ind])
        
        
        kvec_long = np.zeros(2*len(ind[1]))
        kvec_long[:len(ind[1])] = k
        kvec_long[len(ind[1]):] = k
        
        kvec = np.zeros(len(ind[1]))
        kvec[:len(ind[1])/2] = kvec_long[:len(ind[1])/2]
        kvec[len(ind[1])/2:] = kvec_long[len(ind[1]):3*len(ind[1])/2]
        
        ind_for_A_marginalization = np.argsort(kvec)
        
        kvec1_ind = ind_for_A_marginalization[np.in1d(ind_for_A_marginalization, np.where(kvec <= truncated_nmax*beatbox.Universe.Deltak), assume_unique=False)]
        
        kvec2_ind = ind_for_A_marginalization[np.in1d(ind_for_A_marginalization, np.where(kvec > truncated_nmax*beatbox.Universe.Deltak), assume_unique=False)]
        
        
        #print kvec1_ind.shape
        #print kvec2_ind.shape
        
        #Reindex the data
        #datamap1 = datamap[kvec1_ind]
        #datamap2 = datamap[kvec2_ind]
        #datamap_reordformarg = datamap[ind_for_A_marginalization]
        
        #Initiate the inverse covariance matrix of the prior
        #ind = np.where(beatbox.Universe.kfilter>0)
        PS_long = np.zeros(2*len(ind[1]))
        PS_long[:len(ind[1])] = (self.all_simulated_universes[0].Power_Spectrum[ind])
        PS_long[len(ind[1]):] = (self.all_simulated_universes[0].Power_Spectrum[ind])
        
        PS = np.zeros(len(ind[1]))
        PS[:len(ind[1])/2] = PS_long[:len(ind[1])/2]
        PS[len(ind[1])/2:] = PS_long[len(ind[1]):3*len(ind[1])/2]
        
        inv_Cf = np.diag(1./PS)
        
        inv_Cf_reordformarg = inv_Cf[ind_for_A_marginalization, ind_for_A_marginalization] 
        
        #Deal with the R matrix to make it real
        # Select the m values out the the lms tupples
        m = np.array([m[1] for m in beatbox.Universe.lms])
        # Find the indices of the positive ms
        pos_ind = (m>0)
        # Find the indices of the m=0
        zero_ind = (m==0)
        # Find the indices of the negative ms
        neg_ind = (m<0)
        
        R_real = np.zeros((len(beatbox.Universe.lms), len(self.all_simulated_universes[-1].fn)), dtype=np.float)
        
        R_real[pos_ind,:] = beatbox.Universe.R[pos_ind,:].real
        R_real[neg_ind,:] = beatbox.Universe.R[pos_ind,:].imag
        R_real[zero_ind,:] = beatbox.Universe.R[zero_ind,:].astype(np.float)
        
        self.R_real = R_real
        #R_real_reordformarg = R_real[ind_for_A_marginalization]
        #self.R_real_reordformarg = R_real_reordformarg
        
        # Define the A matrix and it's four components, and the needed inverses 
        if NewAcov is None:
            A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
        
            A1 = np.zeros((len(kvec1_ind), len(kvec1_ind)))
            A2 = np.zeros((len(kvec1_ind), len(kvec2_ind)))
            A3 = np.zeros((len(kvec2_ind), len(kvec1_ind)))
            A4 = np.zeros((len(kvec2_ind), len(kvec2_ind)))
            print A1.shape 
            print A2.shape
            print A3.shape
            print A4.shape 
        
        
            for i in range(len(kvec1_ind)):
                A1[i,:] = A[kvec1_ind[i], kvec1_ind]
            for i in range(len(kvec1_ind)):
                A2[i,:] = A[kvec1_ind[i], kvec2_ind]
            print A2.shape
            print  A2[:5,:5]
            for i in range(len(kvec1_ind)):
                A3[:,i] = A[kvec2_ind, kvec1_ind[i]]
            for i in range(len(kvec2_ind)):
                A4[:,i] = A[kvec2_ind, kvec2_ind[i]]
            print A1.shape 
            print A2.shape
            print A3.shape
            print A4.shape 
            
            U, s, V_star = np.linalg.svd(A4)
            inv_A4 = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
            
            #Solve for the normalization of the prior
            #inv_Cf = self.solve_for_prior_normalization(inv_Cyy, inv_Cf,   A, inv_A, R_real, datamap, print_alpha)
        
            # Redefine A with the properly normalized prior
            #A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
            #self.A = A
            # Find the inverse of A that has been the properly normalized  Cf
            #U, s, V_star = np.linalg.svd(A)
            #inv_A = np.dot(V_star.T, np.dot(np.diag(1./s),U.T))
            self.inv_A4=inv_A4
            
            NewAcov = (A1-np.dot(np.dot(A2, inv_A4), A2.T ) )
            self.NewAcov = NewAcov
            count = count+1
        else:
            invA4dotA2t=self.invA4dotA2t
#            A = np.dot(R_real.T , np.dot( inv_Cyy , R_real)) + inv_Cf
#        
#            A2 = np.zeros((len(kvec1_ind), len(kvec2_ind)))
#            for i in range(len(kvec1_ind)):
#                A2[i,:] = A[kvec1_ind[i], kvec2_ind]
#            print A2.shape
            
            
        if inv_NewAcov is None:
            U2, s2, V_star2 = np.linalg.svd(NewAcov)
            inv_NewAcov = np.dot(V_star2.T, np.dot(np.diag(1./s2),U2.T))
            self.inv_NewAcov=inv_NewAcov
        
        #This should not be done:
        #inv_Cyy_reordformarg = np.zeros((len(ind_for_A_marginalization), len(ind_for_A_marginalization)))
        
        #print len(ind_for_A_marginalization)
        #print inv_Cyy_reordformarg.shape
        
        #for i in range(len(ind_for_A_marginalization)):
        #    inv_Cyy_reordformarg[i,:] = inv_Cyy[ind_for_A_marginalization[i], ind_for_A_marginalization]
        
        # Use linear algebra to solve the A*f_n=b linear equation
        #b_reordformarg =  np.dot(R_real_reordformarg.T , np.dot (inv_Cyy_reordformarg , datamap_reordformarg) )
        b =  np.dot(R_real.T , np.dot (inv_Cyy , datamap) )
        b1 = b[kvec1_ind]
        b2 = b[kvec2_ind]
        #self.reconstrunct_fn = np.linalg.solve(A, b)
    
        #from numpy.linalg import inv
        if count is not 0:
            self.invA4dotA2t = np.dot(inv_A4 , A2.T )
            invA4dotA2t = self.invA4dotA2t
            self.reconstrunct_fn_ordformarg = np.dot( inv_NewAcov , b1.T - np.dot( b2.T,  self.invA4dotA2t ) )
            
        else:
            self.reconstrunct_fn_ordformarg = np.dot( inv_NewAcov , b1.T - np.dot( b2.T,  invA4dotA2t ) )
        
        self.reconstrunct_fn = np.zeros(len(kvec))
        
        self.reconstrunct_fn[kvec1_ind] = self.reconstrunct_fn_ordformarg
        self.reconstrunct_fn[kvec2_ind] = np.zeros(kvec2_ind.shape)
        
        return

    
    
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
        
        Delta_fn = fn_true.reshape(len(fn_true),1)-fn_rec.reshape(len(fn_rec),1)
        chi2value = np.dot (Delta_fn.T , np.dot( self.A, Delta_fn  ))
        
        p_value = 1-chi2.cdf(chi2value, len(fn_true)) 
        
        return p_value, chi2value
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    