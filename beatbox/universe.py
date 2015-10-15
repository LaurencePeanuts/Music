import numpy as np
import matplotlib.pylab as plt
import healpy as hp
import string
import yt
import os
import glob
from PIL import Image as PIL_Image
from images2gif import writeGif

import beatbox

# ====================================================================

class Universe(object):
    """
    A simple model universe in a box.
    """
    def __init__(self):
        self.PIXSCALE = 0.1
        self.BOXSIZE = 4.0
        self.NPIX = int(self.BOXSIZE/self.PIXSCALE) + 1
        Nj = np.complex(0.0,self.NPIX)
        self.x, self.y, self.z = np.mgrid[-self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj, -self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj, -self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj]
        # self.phi = np.zeros([3,self.NPIX])
        self.phi = self.x * 0.0
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


    def show_CMB_T_map(self,from_perspective_of="observer"):

        if self.Tmap is None:
            self.NSIDE = 16
            self.Tmap = hp.alm2map(self.alm,self.NSIDE)

        if from_perspective_of == "observer":
            # Sky map:
            hp.mollview(self.Tmap,title="CMB temperature fluctuations as seen from inside the LSS")
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
        hp.mollview(projected_map,min=-max,max=max)
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


    def get_alm(self,l=None,m=None):
        """
        hp.map2alm only returns the positive m coefficients - we need
        to derive the negative ones ourselves if we are going to
        do anything with them outside healpy. See
        http://stackoverflow.com/questions/30888908/healpy-map2alm-function-does-not-return-expected-number-of-alm-values?lq=1
        for discussion.
        """
        if l is None or m is None:
            return None

        elif m >= 0:
            index = hp.Alm.getidx(self.lmax, l, m)
            prefactor = 1.0
            value = self.alm[index]
        else:
            index = hp.Alm.getidx(self.lmax, l, -m)
            prefactor = (-1.0)**m
            value = np.conjugate(self.alm[index])

        return prefactor * value


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

    def generate_a_random_potential_field(self, high_k_cutoff=6, low_k_cutoff=2, nmax=10, n_s=0.96, kstar=0.02, PSnorm=2.43e-9):

        self.nmax = nmax

		self.high_k_cutoff=high_k_cutoff
		self.low_k_cutoff=low_k_cutoff

		self.Deltak=2*np.pi/self.BOXSIZE;
		self.kmax=self.nmax*self.Deltak;
		self.N=2*self.nmax;
		
		self.kx, self.ky, self.kz = np.meshgrid(np.linspace(-kmax,kmax,N+1),np.linspace(-kmax,kmax,N+1),np.linspace(-kmax,kmax,N+1), indexing='ij');
		self.k = np.sqrt(np.power(self.kx, 2)+np.power(self.ky,2)+np.power(self.kz,2));
		
		# Define filter in k-space
		low_k_filter = (~(self.k <= self.low_k_cutoff)).astype(int)
		high_k_filter = (~(self.k >= self.high_k_cutoff)).astype(int)
		self.filter = high_k_filter*low_k_filter
		
		# Define the constants that go in the power spectrum
		# scalar spectral index
		self.n_s=n_s
		# Power spectrum normalization
		self.PSnorm=PSnorm
		# Change units of the pivot scale kstar from Mpc^-1 to normalize the smallest k
		# mode to 1 (i.e. the radius of the CMB photosphere at 13.94Gpc)
		self.kstar=kstar*1.394e4
		
        # Draw Gaussian random Fourier coefficients with a k^{-3+(n_s-1)} power spectrum:
		self.Power_Spectrum = self.PSnorm*np.power((self.k/self.kstar) ,(-3+(self.n_s-1)))
		self.Power_Spectrum[np.isinf(self.Power_Spectrum)]=10**-9

		fn_R = np.random.normal(0, np.sqrt(self.Power_Spectrum/2) )*np.power(self.filter,2)
		fn_I = np.random.normal(0, np.sqrt(self.Power_Spectrum/2) )*np.power(self.filter,2)

		#Need to ensure that f_-k=f^*_k!!!!!
		FT=fn_R+fn_I*1j
		
		X=np.concatenate((np.append(FT[:nmax, nmax+1 ,nmax+1 ], 0), np.conjugate(np.flipud(FT[:nmax, nmax+1 ,nmax+1 ]))), axis=0)
		Z=np.concatenate( ( FT[:, :nmax ,nmax ], X.reshape(21,1), np.conjugate(np.fliplr(np.flipud(FT[:, :nmax ,nmax ])))), axis=1 )
		self.fn= np.concatenate( (FT[:,:,:nmax], Z.reshape(21,21,1), np.conjugate( np.fliplr(np.flipud(FT[:,:,:nmax])))[:,:,::-1] ), axis=2  )

        print "Generated ",self.fn[~(self.fn[:,:,:] == 0)].size," potential Fourier coefficients"


        # Evaluate it on our Phi grid:
        self.evaluate_potential_given_fourier_coefficients()

        return


    def evaluate_potential_given_fourier_coefficients(self):

        
        i = 0

        #for nx in range(-self.nmax,self.nmax+1):
        #   for ny in range(-self.nmax,self.nmax+1):
        #       for nz in range(-self.nmax,self.nmax+1):
        #           if (nx*nx+ny*ny+nz*nz <= self.nmax*self.nmax and not np.all(nx,ny,nz)):
        #               kx,ky,kz = self.Deltak*self.klst[i]
        phase = self.kx * self.x + self.ky * self.y + self.kz * self.z
        self.phi += self.fn * np.complex(np.cos(phase),np.sin(phase))
                        # Bug: this potential is not Real. Need to make coeffs obey f_k = f_-k*
                        #i += 1


        print " Built potential grid, with dimensions ",self.phi.shape,\
              " and mean value ", round(np.mean(self.phi)),"+/-",round(np.std(self.phi))

        return


    def show_potential_with_yt(self,output='phi.png',angle=np.pi/4.0):
        """
        Visualize the gravitational potential using yt. We're after something
        like http://yt-project.org/doc/_images/vr_sample.jpg - described
        at http://yt-project.org/doc/visualizing/volume_rendering.html
        """

        # Load the potential field into a yt data structure,
        # offsetting such that minimum value is zero:
        offset = np.abs(np.min(self.phi))
        ds = yt.load_uniform_grid(dict(phi=self.phi+offset), self.phi.shape)

        # Here's Sam's gist, from https://gist.github.com/samskillman/0e574d1a4f67d3a3b1b1
        #   im, sc = yt.volume_render(ds, field='phi')
        #   sc.annotate_domain(ds)
        #   sc.annotate_axes()
        #   im = sc.render()
        #   im.write_png(output, background='white')
        # volume_render is not yet available, though.

        # Following the example at http://yt-project.org/doc/visualizing/volume_rendering.html

        # Set minimum and maximum of plotting range, taking offset
        # into account so that zero appears in center of color map:
        mi, ma = ds.all_data().quantities.extrema('phi')
        print "Extrema of ds phi:",mi,ma
        mi -= offset
        ma -= offset
        ma = np.max(np.abs([mi,ma]))
        mi = -ma + offset
        ma += offset
        print "Extrema after symmetrizing:", mi,ma

        # Instantiate the ColorTransferFunction.
        tf = yt.ColorTransferFunction((mi, ma))

        # Set up the camera parameters: center, looking direction, width, resolution
        c = (ds.domain_right_edge + ds.domain_left_edge)/2.0

        Lx = np.sqrt(2.0)*np.cos(angle)
        Ly = np.sqrt(2.0)*np.sin(angle)
        Lz = 1.0
        L = np.array([Lx, Ly, Lz])
        W = ds.quan(1.6, 'unitary')
        N = 256

        # Create a camera object
        cam = ds.camera(c, L, W, N, tf, fields=['phi'])

        # Now let's add some isocontours, and take a snapshot:
        tf.add_layers(21, colormap='BrBG')
        im = cam.snapshot(output)

        # Add the domain box to the image:
        nim = cam.draw_domain(im)

        # Save the image to a file:
        nim.write_png(output)

        return


    def show_potential_from_all_angles_with_yt(self,output='phi.gif'):

        # Create 36 frames for the animated gif, one for each angle:
        steps = 36
        angles = np.arange(steps)*2.0*np.pi/np.float(steps)

        # book-keeping:
        folder = 'frames/'
        os.system("rm -rf "+folder)
        os.system("mkdir -p "+folder)

        # Now create the individual frames:
        for k,angle in enumerate(angles):
            framefile = folder+str(k).zfill(3)+'.png'
            print "Making frame",k,": ",framefile,"at viewing angle",angle
            self.show_potential_with_yt(output=framefile,angle=angle)

        # Create an animated gif of all the frames:
        images = [PIL_Image.open(framefile) for framefile in glob.glob(folder+'*.png')]
        writeGif(output, images, duration=0.2)

        return


    def transform_3D_potential_fourier_series_into_T_spherical_harmonics(self,lmax=3):

        self.lmax = lmax
        Nalm = 0
        for l in range(self.lmax+1):
            for m in range(-l,l+1):
                Nalm += 1

        # Dummy code until I figure out Roger's matrix:
        self.R = np.zeros([Nalm,len(self.fn)])

        Re = np.multiply(self.R,self.fn)
        Im = Re.copy()
        self.alm = np.flatten(Re + 1.0j*Im)
        # This does not work - healpy book-keeping is wrong:
        #   TypeError: The a_lm must be a 1D array.

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
