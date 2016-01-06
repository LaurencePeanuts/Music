import numpy as np
import matplotlib
import matplotlib.pylab as plt
import healpy as hp
import string
import yt
import os
import glob
from PIL import Image as PIL_Image
from images2gif import writeGif
from scipy.special import sph_harm,sph_jn


import beatbox

# ===================================================================

def set_k_filter(self):
    """
    Define a filter over the k space for the modes between kmin and kmax 
    """
    #Define lower & upper bounds for the filter
    Universe.high_k_cutoff=Universe.truncated_nmax*Universe.Deltak
    Universe.low_k_cutoff=Universe.truncated_nmin*Universe.Deltak
    
    
    # Define the filter
    #low_k_filter = (~(Universe.n < Universe.low_k_cutoff)).astype(int)
    #high_k_filter = (~(Universe.n > Universe.high_k_cutoff)).astype(int)
    #Universe.kfilter = high_k_filter*low_k_filter
    
    # Define the filter
    low_k_filter = (~(Universe.n < Universe.truncated_nmin)).astype(int)
    high_k_filter = (~(Universe.n > Universe.truncated_nmax)).astype(int)
    Universe.kfilter = high_k_filter*low_k_filter
    return

def populate_response_matrix(self):
    """
    Populate the R matrix for the default range of l and n, or
    or over the range specified above
    """


    truncated_nmax=Universe.truncated_nmax
    truncated_nmin=Universe.truncated_nmin
    truncated_lmax=Universe.truncated_lmax
    truncated_lmin=Universe.truncated_lmin
    lms=Universe.lms
    kfilter=Universe.kfilter
    
    # Initialize R matrix:
    NY = (truncated_lmax + 1)**2-(truncated_lmin)**2
    # Find the indices of the non-zero elements of the filter
    ind=np.where(Universe.kfilter>0)
    # The n index spans 2x that length, 1st half for the cos coefficients, 2nd half
    #    for the sin coefficients
    NN = 2*len(ind[1])
    Universe.R = np.zeros([NY,NN], dtype=np.complex128)

    k, theta, phi = Universe.k[ind], np.arctan2(Universe.ky[ind],Universe.kx[ind]), np.arccos(Universe.kz[ind]/Universe.k[ind])
    # We need to fix the 'nan' theta element that came from having ky=0
    theta[np.isnan(theta)] = np.pi/2.0
    
    # Get ready to loop over y
    y=0
    A=[sph_jn(truncated_lmax,ki)[0] for ki in k]        
    # Loop over y, computing elements of R_yn 
    for i in lms:        
        l=i[0]
        m=i[1]
        
        trigpart = np.cos(np.pi*l/2.0)
        B=np.asarray([A[ki][l] for ki in range(len(k))])
        Universe.R[y,:NN/2] = 4.0 * np.pi * sph_harm(m,l,theta,phi).reshape(NN/2)*B.reshape(NN/2) * trigpart
        trigpart = np.sin(np.pi*l/2.0)
        Universe.R[y,NN/2:] = 4.0 * np.pi * sph_harm(m,l,theta,phi).reshape(NN/2)*B.reshape(NN/2)* trigpart
                
        y=y+1    
        return


# ====================================================================

class Universe(object):
    """
    A simple model universe in a box.
    """
    
    # ====================================================================

    #Initialize the class variables
    PIXSCALE = 0.1
    BOXSIZE = 4.0

    # Real space: define a coordinate grid:
    NPIX = int(BOXSIZE/PIXSCALE) + 1
    Nj = np.complex(0.0,NPIX)
    x, y, z = np.mgrid[-BOXSIZE/2.0+BOXSIZE/(2*float(NPIX)):BOXSIZE/2.0-BOXSIZE/(2*float(NPIX)):Nj, -BOXSIZE/2.0+BOXSIZE/(2*float(NPIX)):BOXSIZE/2.0-BOXSIZE/(2*float(NPIX)):Nj, -BOXSIZE/2.0+BOXSIZE/(2*float(NPIX)):BOXSIZE/2.0-BOXSIZE/(2*float(NPIX)):Nj]
    
    # Define the truncatad range of modes (in n and l) we want in our Universe:
    try:
        truncated_nmax=Multiverse.truncated_nmax
        truncated_nmin=Multiverse.truncated_nmin
        truncated_lmax=Multiverse.truncated_lmax
        truncated_lmin=Multiverse.truncated_lmin
    except NameError:
        truncated_nmax=2
        truncated_nmin=1
        truncated_lmax=8
        truncated_lmin=0
    

    # If only truncated_lmax is provided, calculated the largest truncated_nmax we can reconstruct
    if (truncated_lmax is not None) and (truncated_nmax is None):
        truncated_nmax=int(np.floor((3.0*(truncated_lmax+1)**2.0/(4.0*np.pi))**(1.0/3.0)))
    # Else define a default value for truncated_nmax if not already done
    elif truncated_nmax is None:
        truncated_nmax=6
    # If only truncated_nmax is provided, calculated the truncated_lmax needed for no information
    #    from the 3D map to be lost
    if (truncated_nmax is not None) and (truncated_lmax is None):
        truncated_lmax=int(np.ceil(-0.5+2.0*truncated_nmax**(3.0/2.0)*np.sqrt(np.pi/3.0)))
    
    # Make a y_max-long tupple of l and m pairs
    if None not in (truncated_lmin, truncated_lmax):
        lms=[(l, m) for l in range(truncated_lmin,truncated_lmax+1) for m in range(-l, l+1)]
    
    
    # Fourier space: define a coordinate grid:

    # The nmax we need for the resolution we want in our Universe is:
    nmax = int(BOXSIZE/(2*PIXSCALE))
    Deltak = 2.0*np.pi/BOXSIZE
    kmax = nmax*Deltak
    kx, ky, kz = np.meshgrid(np.linspace(-kmax,kmax,NPIX),np.linspace(-kmax,kmax,NPIX),np.linspace(-kmax,kmax,NPIX), indexing='ij')
    k = np.sqrt(np.power(kx, 2)+np.power(ky,2)+np.power(kz,2))
    nx, ny, nz = np.meshgrid(np.linspace(-nmax,nmax,NPIX),np.linspace(-nmax,nmax,NPIX),np.linspace(-nmax,nmax,NPIX), indexing='ij');
    n = np.sqrt(np.power(nx, 2)+np.power(ny,2)+np.power(nz,2));
    
    # Define filter in k-space, that contains the modes we want:
    kfilter = None
    set_Universe_k_filter = set_k_filter
    
    #Define and populate the R matrix:
    R = None
    populate_Universe_R = populate_response_matrix

    #==========================================================
    
    def __init__(self):

        # The potential map (pure real):
        self.phi = self.x * 0.0
        # The CMB temperature map:
        self.Tmap = None
        self.NSIDE = None
 
        return

    def __str__(self):
        return "a model universe, containing potential map phi"

    # ----------------------------------------------------------------
    def read_in_CMB_T_map(self,from_this=None):
        if from_this is None:
            print "No CMB T map file supplied."
            self.Tmapfile = None
        else:
            self.Tmapfile = from_this
            self.Tmap = hp.read_map(from_this)
            self.NSIDE = hp.npix2nside(len(self.Tmap))
        return


    def show_CMB_T_map(self,Tmap=None,from_perspective_of="observer"):
        if Tmap is None:
            self.NSIDE = 256
            self.Tmap = hp.alm2map(self.alm,self.NSIDE)
        else:
            self.Tmap=Tmap
            
        if from_perspective_of == "observer":
            # Sky map:
            hp.mollview(self.Tmap, coord='G', rot=(-90,0,0),title="CMB temperature fluctuations as seen from inside the LSS")
        else:
            # Interactive "external" view ([like this](http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html))            pass
            #   beatbox.zoncaview(self.Tmap)
            # This did not work, sadly. Maybe we can find a 3D
            # spherical surface plot routine using matplotlib? For
            # now, just use the healpix vis.
            R = (0.0,0.0,0.0) # (lon,lat,psi) to specify center of map and rotation to apply
            hp.orthview(self.Tmap,rot=R,flip='geo',half_sky=True,title="CMB temperature fluctuations as seen from outside the LSS")
            # print "Ahem - we can't visualize maps on the surface of the sphere yet, sorry."
        return


    def decompose_T_map_into_spherical_harmonics(self,lmax=None):
        """
        See healpy documentation at https://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.map2alm.html
        self.alm is a 1D numpy array of type=complexx128.
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
        """
        To do this we need to make a healpy-format alm array, with
        just one non-zero complex value in it, which we extract
        from the parent alm array. Since healpy only returns positive
        m coefficients, we just ask to see component with that |m|.
        """
        projected_alm = self.alm * 0.0
        i = hp.Alm.getidx(self.lmax, l, np.abs(m)) # Note |m| here
        projected_alm[i] = self.alm[i]
        projected_map = hp.alm2map(projected_alm,self.NSIDE)
        hp.mollview(projected_map)
        return


    def show_lowest_spherical_harmonics_of_CMB_T_map(self,lmax=10,max=20):
        """
        To do this, we construct a healpy-formatted alm array based on
        a subset of the parent one, again observing the positive m-only
        convention.
        """
        truncated_alm = self.alm * 0.0
        i = []
        for l in range(lmax+1):
            for m in range(l+1):
                i.append(hp.Alm.getidx(self.lmax, l, m))
        print "Displaying sky map of the l = ",l," and lower spherical harmonics only..."
        truncated_alm[i] = self.alm[i]
        truncated_map = hp.alm2map(truncated_alm,self.NSIDE)
        hp.mollview(truncated_map,min=-max,max=max)
        return


    def get_alm(self,l=None,m=None,lms=None):
        """
        hp.map2alm only returns the positive m coefficients - we need
        to derive the negative ones ourselves if we are going to
        do anything with them outside healpy. See
        http://stackoverflow.com/questions/30888908/healpy-map2alm-function-does-not-return-expected-number-of-alm-values?lq=1
        for discussion.
        """
        if (l is None or m is None) and lms is None:
            return None
        
        elif l is None and m is None:
            ay=np.zeros(len(lms),dtype=np.complex128)
            for i in lms:
                if i[1] >= 0:
                    index = hp.Alm.getidx(self.lmax, i[0], i[1])
                    prefactor = 1.0
                    value = self.alm[index]
                else:
                    index = hp.Alm.getidx(self.lmax, i[0], -i[1])
                    prefactor = (-1.0)**i[1]
                    value = np.conjugate(self.alm[index])
                ay[i[0]**2+i[0]+i[1]-(lms[0][0])**2]=prefactor * value
            return ay

        elif m >= 0:
            index = hp.Alm.getidx(self.lmax, l, m)
            prefactor = 1.0
            value = self.alm[index]
        else:
            index = hp.Alm.getidx(self.lmax, l, -m)
            prefactor = (-1.0)**m
            value = np.conjugate(self.alm[index])

        return prefactor * value

    def put_alm(self,value,l=None,m=None,lms=None):
        '''
        Re-arranges the value or vector of a_y values into the
        correct order to be used by healpy as a_lm.
        If lms is given, len(lms) must equal len(value), while
        if l and m are specified, value must be a scalar.
        '''
        print 'made it to put_alm!!'
        if (l is None or m is None) and lms is None:
            return None
        elif l is None and m is None:
            if len(lms) != len(value):
                print 'a_y and (l, m) are of unequal lenghts, cannot proceed'
                return
            index=np.zeros(len(lms), dtype=int)
            count=0
            for i in lms:
                index[count] = hp.Alm.getidx(max(lms)[0], i[0], i[1])
                count=count+1
            lmax=max(lms)[0]
            mmax=max(lms)[1]
            self.alm=np.zeros(mmax*(2*lmax+1-mmax)/2+lmax+1, dtype=np.complex128)
            # Throw away the negative indices (which correspond to the negative m's)
            #     since the maps are real, negative m coefficients can be deduced
            #     from the positive ones.
            index_positive=index[~(index<0)]
            ind1=np.arange(len(value))
            self.alm[index_positive] = value[ind1[~(index<0)]]
            return
        index = hp.Alm.getidx(self.truncated_lmax, l, m)
        self.alm[index] = value
        
        return


    def alm2ay(self, truncated_lmax=None, truncated_lmin=None, usedefault=1):
        """
        Read its own a_lm array, and return the corresponding
        a_y array (in the correct order).
        The conversion between y index and l_max is:
        (l+1)**2-(2l+1)/2 +1/2 +m = l**2+2*l+1-l-1/2+1/2+m = l**2+l+1+m
        and the first element has index 0 so subtract 1, so 
        y=l**2+l+m is the index, need to subtract the elements before lmin
        so y=l**2+l+m-(lmin+1)**2
        """
        if usedefault==1:
            truncated_lmax=self.truncated_lmax
            truncated_lmin=self.truncated_lmin
            lms=self.lms
        # Make a y_max-long tupple of l and m pairs
        else:
            lms=[(l, m) for l in range(truncated_lmin,truncated_lmax+1) for m in range(-l, l+1)]
        
        ay = np.zeros((truncated_lmax+1)**2-(truncated_lmin)**2,dtype=np.complex128)
        ay=self.get_alm(lms=lms)
        
        #self.ay=ay
        return ay


    def ay2alm(self, ay,truncated_lmax=None, truncated_lmin=None, usedefault=1):
        """
        Repackage the a_y array into healpy-readable a_lm's
        """
        if usedefault==1:
            truncated_lmax=self.truncated_lmin
            truncated_lmin=self.truncated_lmax
            lms=self.lms
        
        # Make a y_max-long tupple of l and m pairs
        else:
            lms=[(l, m) for l in range(truncated_lmin,truncated_lmax+1) for m in range(-l, l+1)]
        
        self.put_alm(ay, lms=lms)
        
        return
    
    def ay2ayreal_for_inference(self,value):
        '''
        Reorganize the ays so that only independent measurements are kept,
        and split the real and imaginary values into different elements.
        The negative m values ara dependent on the positive m, so they must 
        be discarted, and all but m=0 values are complex.
        Therefore, we replace the positive m values by their respective real
        part, and the negative m values by the imaginary part of the
        corresponding positive m. This way, each l retains 2l+1 independent
        real degrees of freedom.
        '''
        
        #Select the m values out the the lms tupples
        m=np.array([m[1] for m in self.lms])
        #Find the indices of the positive ms
        pos_ind=(m>0)
        #Find the indices of the m=0
        zero_ind=(m==0)
        #Find the indices of the negative ms
        neg_ind=(m<0)
        
        ay_real=np.zeros(len(self.lms), dtype=np.float)
        
        ay_real[pos_ind]=value[pos_ind].real
        ay_real[neg_ind]=value[pos_ind].imag
        ay_real[zero_ind]=value[zero_ind].astype(np.float)
        
        
        return ay_real

    def write_out_spherical_harmonic_coefficients(self,lmax=10):
        outfile = string.join(string.split(self.Tmapfile,'.')[0:-1],'.') + '_alm_lmax' + str(lmax) + '.txt'
        f = open(outfile, 'w')
        f.write("#    l    m    alm_real   alm_imag\n")
        count = 0
        for l in range(lmax+1):
            for m in range(-l,l+1):
                alm = self.get_alm(l,m)
                line = " {0:d}  {1:d}  {2:g}  {3:g}\n".format(l,m,float(np.real(alm)),float(np.imag(alm)))
                f.write(line)
                count += 1
        f.close()
        print count,"alm's (lmax =",lmax,") written to",outfile
        return

    # ----------------------------------------------------------------

    def populate_instance_response_matrix(self,truncated_nmax=None, truncated_nmin=None,truncated_lmax=None, truncated_lmin=None, usedefault=1):
        """
        Populate the R matrix for the default range of l and n, or
        or over the range specified above
        """

        if usedefault==1:
            truncated_nmax=self.truncated_nmax
            truncated_nmin=self.truncated_nmin
            truncated_lmax=self.truncated_lmax
            truncated_lmin=self.truncated_lmin
            lms=self.lms
            kfilter=self.kfilter    
        else:
            low_k_cutoff=truncated_nmin*self.Deltak
            high_k_cutoff=truncated_nmax*self.Deltak
            self.set_instance_k_filter(truncated_nmax=truncated_nmax,truncated_nmin=truncated_nmin)    
            lms=[(l, m) for l in range(truncated_lmin,truncated_lmax+1) for m in range(-l, l+1)]
        
        
        # Initialize R matrix:
        NY = (truncated_lmax + 1)**2-(truncated_lmin)**2
        # Find the indices of the non-zero elements of the filter
        ind=np.where(self.kfilter>0)
        # The n index spans 2x that length, 1st half for the cos coefficients, 2nd half
        #    for the sin coefficients
        NN = 2*len(ind[1])
        self.R = np.zeros([NY,NN], dtype=np.complex128)

        # In case we need n1, n2, n3 at some point...:
        #    n1, n2, n3 = self.kx[ind]/self.Deltak , self.ky[ind]/self.Deltak, self.kz[ind]/self.Deltak
        k, theta, phi = self.k[ind], np.arctan2(self.ky[ind],self.kx[ind]), np.arccos(self.kz[ind]/self.k[ind])
        # We need to fix the 'nan' theta element that came from having ky=0
        theta[np.isnan(theta)] = np.pi/2.0
        
        # Get ready to loop over y
        y=0
        A=[sph_jn(truncated_lmax,ki)[0] for ki in k]        
        # Loop over y, computing elements of R_yn 
        for i in lms:        
            l=i[0]
            m=i[1]
            trigpart = np.cos(np.pi*l/2.0)
            B=np.asarray([A[ki][l] for ki in range(len(k))])
            self.R[y,:NN/2] = 4.0 * np.pi * sph_harm(m,l,theta,phi).reshape(NN/2)*B.reshape(NN/2) * trigpart

            trigpart = np.sin(np.pi*l/2.0)
            self.R[y,NN/2:] = 4.0 * np.pi * sph_harm(m,l,theta,phi).reshape(NN/2)*B.reshape(NN/2)* trigpart
                
            y=y+1   
        return

    # ----------------------------------------------------------------
    
    def load_mathematica_data(self):
        
        f= open("data/f_ns.txt", 'r')
        data = f.read()
        f.close()
        columns = data.split()
        f_n=np.zeros(len(columns))
        for count in range(int(len(columns))):
            f_n[count] = float(columns[count])
        
        g= open("data/fncoordinates.txt", 'r')
        data2 = g.read()
        g.close()
        columns2 = data2.split()
        k_vec=np.zeros(len(columns2))
        for count2 in range(int(len(columns2))):
            k_vec[count2] = float(columns2[count2])
        
        k_x=k_vec[0::3]
        k_y=k_vec[1::3]
        k_z=k_vec[2::3]

        return f_n, k_x, k_y, k_z
        
    # ----------------------------------------------------------------

    def set_instance_k_filter(self,truncated_nmax=None,truncated_nmin=None):
        """
        Define a filter over the k space for the modes between kmin and kmax 
        """
        #Make sure we have lower & upper bounds for the filter
        if truncated_nmax is None:
            self.high_k_cutoff=self.truncated_nmax*self.Deltak
        else:
            self.truncated_nmax = truncated_nmax
            self.high_k_cutoff=truncated_nmax*self.Deltak
        if truncated_nmin is None:
            self.low_k_cutoff=self.truncated_nmin*self.Deltak
        else:
            self.truncated_nmin = truncated_nmin
            self.low_k_cutoff = truncated_nmin*self.Deltak
        
        # Define the filter
        low_k_filter = (~(self.n < self.truncated_nmin)).astype(int)
        high_k_filter = (~(self.n > self.truncated_nmax)).astype(int)
        self.kfilter = high_k_filter*low_k_filter
        return


    def generate_a_random_potential_field(self,truncated_nmax=6,truncated_nmin=2,n_s=0.97,kstar=0.02,PSnorm=2.43e-9,Pdist=1,Pmax=np.pi,Pvar=0.0):

        #is this realy necessary since filter def moved up in __init__ function??
        # Set the k filter:
        if (beatbox.Universe.kfilter is None) or (truncated_nmax != self.truncated_nmax) or (truncated_nmin != self.truncated_nmin):
            self.set_instance_k_filter(truncated_nmax=truncated_nmax,truncated_nmin=truncated_nmin)

        # Define the constants that go in the power spectrum
        #    scalar spectral index
        self.n_s = n_s
        #   power spectrum normalization
        self.PSnorm = PSnorm
        # Change units of the pivot scale kstar from Mpc^-1 to normalize the smallest k
        #    mode to 1 (i.e. the radius of the CMB photosphere at 13.94Gpc)
        self.kstar = kstar*1.394e4

        # Draw Gaussian random Fourier coefficients with a k^{-3+(n_s-1)} power spectrum:
        self.Power_Spectrum = self.PSnorm*np.power((self.k/self.kstar) ,(-3+(self.n_s-1)))
        self.Power_Spectrum[np.isinf(self.Power_Spectrum)] = 10**-9

        fn_Norm = np.random.normal(0, np.sqrt(self.Power_Spectrum) )*np.power(self.kfilter,2)
        # Draw the phases for the modes: use p=1 for a uniform distribution in [0,Pmax],
        #    and p=0 for a Gaussian distribution with mean Pmax and variance Pvar
        self.Pdist=Pdist
        self.Pvar=Pvar
        self.Pmax=Pmax

        if Pdist==1:
            fn_Phase = np.random.uniform(0, Pmax*np.ones(self.k.shape,dtype=np.float_) )*np.power(self.kfilter,2)
        else:
            fn_Phase = np.random.normal(Pmax, np.sqrt(Pvar)*np.ones(self.k.shape,dtype=np.float_) )*np.power(self.kfilter,2)

        # Need to ensure that f_-k = f^*_k
        # FT = fn_R + fn_I*1j
        FT = fn_Norm*np.cos(fn_Phase)+fn_Norm*np.sin(fn_Phase)*1j

        X = np.concatenate((np.append(FT[:self.nmax, self.nmax ,self.nmax ], 0), np.conjugate(FT[:self.nmax, self.nmax ,self.nmax ])[::-1]), axis=0)
        Z = np.concatenate( ( FT[:, :self.nmax ,self.nmax ], X.reshape(2*self.nmax+1,1), np.conjugate(FT[:, :self.nmax ,self.nmax ])[::-1,::-1]), axis=1 )
        self.fngrid = np.concatenate( (FT[:,:,:self.nmax], Z.reshape(2*self.nmax+1,2*self.nmax+1,1), np.conjugate( FT[:,:,:self.nmax])[::-1,::-1,::-1] ), axis=2  )

        print "Generated ",self.fngrid[~(self.fngrid[:,:,:] == 0)].size," potential Fourier coefficients"

        if Pdist==1:
            print " with phases uniformly distributed between 0 and ", Pmax
        else:
            print " with phases sampled from a Gaussian distribution with mean ", Pmax," and variance ", Pvar

        # Evaluate it on our Phi grid:
        self.evaluate_potential_given_fourier_coefficients()

        return

    
    def evaluate_potential_given_fourier_coefficients(self):

        self.phi = np.zeros(self.x.shape,dtype=np.float_)
        ComplexPhi = np.zeros(self.x.shape,dtype=np.complex128)

        #THIS PART DID THE iFFT MANUALLY
        # for i in range((2*self.nmax+1)**3):
        #    phase = self.kx.reshape((2*self.nmax+1)**3,1)[i] * self.x + self.ky.reshape((2*self.nmax+1)**3,1)[i] * self.y + self.kz.reshape((2*self.nmax+1)**3,1)[i] * self.z
        #    ComplexPhi += self.fn.reshape((2*self.nmax+1)**3,1)[i] * (np.cos(phase)+np.sin(phase)*1j)

        #Now use iFFT to invert the Fourier coefficients f_n to a real space potential
        ComplexPhi = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(self.fngrid)))

        # Throw out the residual imaginary part of the potential [< O(10^-16)]
        self.phi = ComplexPhi.real

        print "Built potential grid, with dimensions ",self.phi.shape,\
              " and mean value ", round(np.mean(self.phi),4),"+/-",round(np.std(self.phi),7)

        return

    
    def rearrange_fn_from_grid_to_vector(self):
        '''
        It's easiest to generate a potential from the prior on a 3D
        grid, so we can use the iFFT. For the linear algebra in the
        inference, we need the fourier coefficients arranged in a
        vector.
        '''
        ind=np.where(self.kfilter>0)
        self.fn=np.zeros(2*len(ind[1]))
        self.fn[:len(ind[1])] = (self.fngrid[ind]).real
        self.fn[len(ind[1]):] = (self.fngrid[ind]).imag
        return

    def transform_3D_potential_into_alm(self, truncated_nmax=None, truncated_nmin=None,truncated_lmax=None, truncated_lmin=None, usedefault=1):
        '''
        From the f_n on a 3D grid, rearrange the Fourier coefficients 
        in a vector and generate the R matrix. From these, calculate the a_y 
        and finally rearrange them in a a_lm vector useable by healpy to 
        make a T map.
        The method can do this either for the harmonics correcponding to the
        full range of n values of the 3D potential (if usedefault=1), or else 
        for the specified values. If truncated_nmax is too large for the
        specified truncated_lmax, some information will be lost.        
        '''
        # Make a vector out of the fn grid of Fourier coefficients
        self.rearrange_fn_from_grid_to_vector()
        if usedefault==1:
            # Populate the R matrix
            if beatbox.Universe.R is None:
                self.populate_instance_response_matrix(usedefault=usedefault)
            # Calculate the a_y matrix
            ay=np.dot(self.R,self.fn)
            self.ay=ay
            # Reorganize a_y into a_lm
            self.ay2alm(ay, usedefault=usedefault)
        else:
            # Populate the R matrix
            self.populate_instance_response_matrix(truncated_nmax=truncated_nmax, truncated_nmin=truncated_nmin,truncated_lmax=truncated_lmax, truncated_lmin=truncated_lmin, usedefault=0)        
            # Calculate the a_y matrix
            ay=np.dot(self.R,self.fn)
            self.ay=ay
            # Reorganize a_y into a_lm
            self.ay2alm(ay,truncated_lmax=truncated_lmax, truncated_lmin=truncated_lmin, usedefault=0)

        return


    def show_potential_with_yt(self,output='',angle=np.pi/4.0, N_layer=5, alpha_norm=1.0, cmap='BrBG', Proj=0, Slice=0, gifmaking=0, show3D=0):
        """
        Visualize the gravitational potential using yt. We're after something
        like http://yt-project.org/doc/_images/vr_sample.jpg - described
        at http://yt-project.org/doc/visualizing/volume_rendering.html
        """

        # Load the potential field into a yt data structure,
        #   offsetting such that minimum value is zero.
        # First get extrema of phi array:
        mi = np.min(self.phi)
        ma = np.max(self.phi)
        # Symmetrize to put zero at center of range:
        ma = np.max(np.abs([mi,ma]))
        mi = -ma
        # Offset to make minimum value zero:
        offset = ma
        ma = 2.0*ma
        mi = 0.0

        # Size of the box containing the phi
        # Physical -2 to 2 box
        # bbox = np.array([[-2, 2], [-2, 2], [-2, 2]])
        # Physical box from the iFFT 
        bbox = np.array([[np.min(self.x), np.max(self.x)], [np.min(self.y), np.max(self.y)], [np.min(self.z), np.max(self.z)]])

        # Apply offset and store phi array in a yt data structure,
        #    I'm putting some random density units here 
        #    (seems to be needed to display properly):
        ds = yt.load_uniform_grid(dict(density=(self.phi+offset, 'g/cm**3')), self.phi.shape, bbox=bbox,  nprocs=1)
        field='density'
        #Check that the loaded field is recognized by yt
        #    print ds.field_list


        # Here's Sam's gist, from https://gist.github.com/samskillman/0e574d1a4f67d3a3b1b1
        #   im, sc = yt.volume_render(ds, field='phi')
        #   sc.annotate_domain(ds)
        #   sc.annotate_axes()
        #   im = sc.render()
        #   im.write_png(output, background='white')
        # volume_render is not yet available, though.

        # Following the example at http://yt-project.org/doc/visualizing/volume_rendering.html

        # Set minimum and maximum of plotting range (in proper yt units):
        dd = ds.all_data()
        mi2, ma2 = dd.quantities.extrema(field)
        #print "Extrema of ds phi:",mi,ma, mi2, ma2

        use_log = False

        # Instantiate the ColorTransferFunction.
        #        tf = yt.ColorTransferFunction((mi2, ma2))
        #        tf.grey_opacity=True
        #    Add some isopotential surface layers:
        #        tf.add_layers(N_layer, 0.0000005*(ma2 - mi2) / N_layer, alpha=alpha_norm*np.ones(N_layer,dtype='float64'), colormap = cmap)

        # Instantiate the ColorTransferFunction using the transfer function helper.
        from IPython.core.display import Image
        from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper

        tfh = yt.TransferFunctionHelper(ds)
        tfh.set_field('density')
        tfh.set_log(False)
        tfh.set_bounds()
        tfh.build_transfer_function()
        tfh.tf.grey_opacity=False
        #For small units, wide Gaussians:
        #tfh.tf.add_layers(N_layer,  w=10*(ma2 - mi2) /N_layer, mi=0.2*ma, ma=ma-0.2*ma, alpha=alpha_norm*np.ones(N_layer,dtype='float64'), col_bounds=[0.2*ma,ma-0.3*ma] , colormap=cmap)
        #For big units, small Gaussians
        tfh.tf.add_layers(N_layer,  w=0.0000001*(ma2 - mi2) /N_layer, mi=0.2*ma, ma=ma-0.2*ma, alpha=alpha_norm*np.ones(N_layer,dtype='float64'), col_bounds=[0.2*ma,ma-0.3*ma] , colormap=cmap)
        # Check if the transfer function captures the data properly:
        densityplot1 = tfh.plot('densityplot1')
        densityplot2 = tfh.plot('densityplot2', profile_field='cell_mass')

        # Set up the camera parameters: center, looking direction, width, resolution
        c = (np.max(self.x)+np.min(self.x))/2.0
        Lx = np.sqrt(2.0)*np.cos(angle)
        Ly = np.sqrt(2.0)*np.sin(angle)
        Lz = 1.0
        L = np.array([Lx, Ly, Lz])
        W = ds.quan(1.6, 'unitary')
        N = 512

#        # Create a camera object
#        cam = ds.camera(c, L, W, N, tfh.tf, fields=[field], log_fields = [use_log],  no_ghost = False)
#        if show3D==1:
#            cam.show()
#        self.cam=cam
#
#        if self.Pdist==1:
#        	im1 = cam.snapshot('opac_phi3D_Uniform_phases_0-'+str(self.Pmax)+'.png', clip_ratio=5)
#        else:
#            im1 = cam.snapshot('opac_phi3D_Gauss_phases_mean'+str(self.Pmax)+'_var'+str(self.Pvar)+'.png', clip_ratio=5)

        if gifmaking==1:
        	# Add the domain box to the image:
        	nim = cam.draw_grids(im1)

        	# Save the image to a file:
        	nim.write_png(output)

        if Proj==1:
            s = yt.ProjectionPlot(ds, "z", "density")
            s.show()
            s.save('phi')

        if Slice==1:
            w = yt.SlicePlot(ds, "x", "density", center="c")
            w.show()
            w.save('phitest')

        return


    def show_potential_from_all_angles_with_yt(self,output='phi.gif'):

        # Create 36 frames for the animated gif, one for each angle:
        steps = 36
        angles = np.arange(steps)*np.pi/np.float(steps)/2.0+np.pi/4
        # current bug: the frames jump at pi/4, 3pi/4 etc..

        # book-keeping:
        folder = 'frames/'
        os.system("rm -rf "+folder)
        os.system("mkdir -p "+folder)

        # Now create the individual frames:
        for k,angle in enumerate(angles):
            framefile = folder+str(k).zfill(3)
            print "Making frame",k,": ",framefile,"at viewing angle",angle
            self.show_potential_with_yt(output=framefile,angle=angle, N_layer=5, alpha_norm=1.0, cmap='BrBG', Proj=0, Slice=0, gifmaking=1)

        # Create an animated gif of all the frames:
        images = [PIL_Image.open(framefile) for framefile in glob.glob(folder+'*.png')]
        writeGif(output, images, duration=0.2)

        return


# ====================================================================

"""
Response matrix from Roger's mathematica notebook:

# Construct the klst:
nmax = 6;
klst = {};
Do[
  If[0 < n1^2 + n2^2 + n3^2 <= nmax^2, klst = Append[klst, {n1, n2, n3}]],
      {n1, -nmax, nmax}, {n2, -nmax, nmax}, {n3, -nmax, nmax}
  ];
NN = Length[klst];

# Set size of box, via separation in k space:
[CapitalDelta]k = .5 [Pi];

# Construct llst, an array of l's and m's for use in Spherical Harmonics:
# Note that the monopole and dipole are ignored!
lmax = 10;
llst = {};
Do[
  If[1 < l <= lmax, llst = Append[llst, {l, m}]], {l, 2, lmax}, {m, -l, l}
  ];
llst; # Not sure what this line does.
L = Length[llst];

# Construct R matrix:
R = Chop[ # Clean out rounding errors (esp in imaginary parts)
     Table[4. [Pi] I^llst[[y, 1]] # i^l - imaginary i!

      SphericalHarmonicY[llst[[y, 1]],
                        llst[[y, 2]],
                        ArcCos[klst[[n, 3]]/Norm[klst[[n]]]], # theta'
                        If[klst[[n, 1]] == klst[[n, 2]] == 0, 0, ArcTan[klst[[n, 1]], klst[[n, 2]]]]] # phi'
                      [Conjugate] # Take complex conjugate of the Ylm

      SphericalBesselJ[llst[[y, 1]], [CapitalDelta]k Norm[klst[[n]]]], # Norm gives the length of the k vector

      {y, 1, L}, # for y in range 1 to L
      {n, 1, NN} # for n in range 1 to NN
     ] # End of Table command
  ];

# Write it out:
(*Export["myn.txt",R]*)

"""
