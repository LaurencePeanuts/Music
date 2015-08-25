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
        self.PIXSCALE = 0.1
        self.BOXSIZE = 4.0
        self.NPIX = int(self.BOXSIZE/self.PIXSCALE) + 1
        Nj = np.complex(0.0,self.NPIX)
        self.x, self.y, self.z = np.mgrid[-self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj, -self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj, -self.BOXSIZE/2.0:self.BOXSIZE/2.0:Nj]
        # self.phi = np.zeros([3,self.NPIX])
        self.phi = self.x * 0.0
        return

    def __str__(self):
        return "a model universe, containing potential map phi"

    # ----------------------------------------------------------------

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
            hp.mollview(self.Tmap,title="CMB temperature fluctuations as seen from inside the LSS")
        else:
            # Interactive "external" view ([like this](http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html))            pass
            #   beatbox.zoncaview(self.Tmap)
            # This did not work, sadly. Maybe we can find a 3D
            # spherical surface plot routine using matplotlib? For
            # now, just use the healpix vis.
            R = (0.0,0.0,0.0) # (lon,lat,psi) to specify center of map and rotation to apply
            hp.orthview(self.Tmap,rot=R,half_sky=True,title="CMB temperature fluctuations as seen from outside the LSS")
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

    # ----------------------------------------------------------------

    def generate_a_random_potential_field(self,nmax=6):

        # Draw unit Gaussian random Fourier coefficients:
        self.nmax = nmax
        self.fn = np.array([])
        self.n = []
        for nx in range(-self.nmax,self.nmax+1):
            for ny in range(-self.nmax,self.nmax+1):
                for nz in range(-self.nmax,self.nmax+1):
                    if (nx*nx+ny*ny+nz*nz <= self.nmax*self.nmax):
                        self.n.append([nx,ny,nz])
                        self.fn = np.append(self.fn,np.random.randn())

        self.n = np.array(self.n)
        print "Generated ",len(self.fn)," potential Fourier coefficients"

        # Evaluate it on our Phi grid:
        self.evaluate_potential_given_fourier_coefficients()

        return


    def evaluate_potential_given_fourier_coefficients(self):

        deltak = 2.0*np.pi/self.BOXSIZE
        i = 0
        for nx in range(-self.nmax,self.nmax+1):
            for ny in range(-self.nmax,self.nmax+1):
                for nz in range(-self.nmax,self.nmax+1):
                    if (nx*nx+ny*ny+nz*nz <= self.nmax*self.nmax):
                        kx,ky,kz = deltak*self.n[i]
                        phase = kx * self.x + ky * self.y + kz * self.z
                        self.phi += self.fn[i]*np.cos(phase)
                        i += 1

        print " Built potential grid, with dimensions ",self.phi.shape,\
              " and mean value ", round(np.mean(self.phi)),"+/-",round(np.std(self.phi))

        return


    def show_potential_with_yt(self):
        """
        Vaporware to make a beautiful visualization of the gravitational
        potential of the entire universe, using yt. We're after something
        like http://yt-project.org/doc/_images/vr_sample.jpg - described
        at http://yt-project.org/doc/visualizing/volume_rendering.html
        Hopefully someone out there can help us out!
        """

        # make_amazing_yt_3D_plot(self.phi)

        # Until then, hang our heads:
        print "We need someone with yt skillz to code this."

        return


    def transform_3D_potential_fourier_series_into_T_spherical_harmonics(self,lmax=3):

        self.lmax = lmax
        Nalm = 0
        for l in range(self.lmax+1):
            for m in range(-l,l+1):
                Nalm += 1

        # Dummy code until I figure out Roger's matrix:
        self.R = np.zeros([Nalm,len(self.fn)])

        self.alm = np.multiply(self.R,self.fn)

        return

# ====================================================================

"""
Response matrix from Roger's mathematica notebook:

nmax = 6; klst = {};
Do[If[0 < n1^2 + n2^2 + n3^2 <= nmax^2,
  klst = Append[klst, {n1, n2, n3}]], {n1, -nmax, nmax}, {n2, -nmax,
  nmax}, {n3, -nmax, nmax}]; NN =
 Length[klst]; \[CapitalDelta]k = .5 \[Pi];
lmax = 10; llst = {}; Do[
 If[1 < l <= lmax, llst = Append[llst, {l, m}]], {l, 2, lmax}, {m, -l,
   l}]; llst; L = Length[llst];
myn = Chop[Table[4. \[Pi] I^llst[[y, 1]]
     SphericalHarmonicY[llst[[y, 1]], llst[[y, 2]],
       ArcCos[klst[[n, 3]]/Norm[klst[[n]]]],
       If[klst[[n, 1]] == klst[[n, 2]] == 0, 0,
        ArcTan[klst[[n, 1]],
         klst[[n, 2]]]]]\[Conjugate] SphericalBesselJ[
      llst[[y, 1]], \[CapitalDelta]k Norm[klst[[n]]]], {y, 1, L}, {n,
     1, NN}]];
(*Export["myn.txt",myn]*)

"""
