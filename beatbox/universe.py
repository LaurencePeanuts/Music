import numpy as np
import matplotlib.pylab as plt
import ipdb

plt.ion()
np.random.seed(3) # for reproduceability
r_cmb_mpc = 14.0
cmap = 'gray'
vscale = 15.

def demo():
     # Generate fake CMB data on a Healpix sphere.
     f = FakeHealpixData()
     f.show()

     # Define a 2d slice through our universe.
     s = SliceSurface()

     # Define an Inference object, then infer and visualize 
     # the minimum-variance phi field on the slice, given data
     # on the sphere.
     inf = Inference(f, s)
     inf.calculate_mv_phi()
     inf.view_phi_mv_slice()

     # Make a bunch of realizations and analyze/visualize them.
     '''

     # RK: realizations have lower variance around the CMB ring (good), but
     # have too-high variance in center of ring.  I think it's an artifact
     # of how I've defined the correlation/covariance function, namely as 
     # an integral that starts at k_min = 2*pi/(2*r_cmb).  Not sure where 
     # to go from here.

     slice_realizations = []
     for i in range(20):
          print i
          this_slice_realization = inf.calculate_phi_realization()
          slice_realizations.append(this_slice_realization)
     slice_realizations = np.array(slice_realizations)
     ipdb.set_trace()
     '''


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


class CartesianCoordinates(object):
    def __init__(self):
        pass

    def update_xyz(self):
        self.xyz = np.vstack([self.x, self.y, self.z]).T

    def make_distance_array(self, other_cart_coord):
        #print '...making distance array...'
        # Fast pairwise distances, see 
        # https://jakevdp.github.io/blog/2013/06/15/numba-vs-cython-take-2/
        from scipy.spatial.distance import cdist
        return cdist(self.xyz, other_cart_coord.xyz)

    def make_auto_distance_array(self):
        return self.make_distance_array(self)


class SliceSurface(CartesianCoordinates):
    def __init__(self, position=0., side_mpc=30., reso_mpc=0.8):
        self.side_mpc = side_mpc
        self.reso_mpc = reso_mpc
        n_side = np.ceil(side_mpc/reso_mpc)
        self.n_side = n_side
        x_2d = self.reso_mpc*np.tile(np.arange(n_side),n_side).reshape(n_side, n_side)
        x_2d -= x_2d.mean()
        y_2d = self.reso_mpc*np.tile(np.arange(n_side),n_side).reshape(n_side, n_side).T
        y_2d -= y_2d.mean()
        z_2d = self.reso_mpc*np.zeros_like(x_2d) + position
        z_2d -= z_2d.mean()
        self.n_side = n_side
        self.x = x_2d.ravel()
        self.y = y_2d.ravel()
        self.z = z_2d.ravel()
        self.update_xyz()



class HealpixSphericalSurface(CartesianCoordinates):
    def __init__(self, radius_mpc=r_cmb_mpc, n_side=2**4):
        # FYI: n_side = 2**4 corresponds to 
        # 0.064 radians resolution = ~0.9 Mpc at z~1000.
        from healpy import nside2npix, pix2vec
        self.n_pix = nside2npix(n_side)
        x, y, z = pix2vec(n_side, range(self.n_pix))
        self.radius_mpc = radius_mpc
        self.x = self.radius_mpc*x
        self.y = self.radius_mpc*y
        self.z = self.radius_mpc*z
        self.update_xyz()
        

class FakeHealpixData(HealpixSphericalSurface):
    def __init__(self, sigma=1e-10):
        HealpixSphericalSurface.__init__(self)
        self.sigma = sigma
        self.data = np.zeros(self.n_pix)
        self.add_truth() 
        self.add_noise()

    def add_truth(self):
        #print '...adding truth...'
        distance = self.make_auto_distance_array()
        delta = distance[distance!=0].min()
        cov = large_scale_phi_covariance(distance)
        from numpy.random import multivariate_normal
        self.data += multivariate_normal(np.zeros(self.n_pix), cov)

    def add_noise(self):
        #print '...adding noise...'
        from numpy.random import randn
        self.data += self.sigma*randn(self.n_pix)
        pass

    def show(self):
        from healpy import mollview
        mollview(self.data)#, cmap=cmap, min=-vscale, max=+vscale)


class HealpixPlusSlice(CartesianCoordinates):
     def __init__(self):
          healpix = HealpixSphericalSurface()
          slice = SliceSurface()
          self.n_healpix = len(healpix.x)
          self.n_slice = len(slice.x)
          self.n_total = self.n_healpix + self.n_slice
          self.ind_healpix = range(0, self.n_healpix)
          self.ind_slice = range(self.n_healpix, self.n_total)
          self.x = np.hstack([healpix.x, slice.x])
          self.y = np.hstack([healpix.y, slice.y])
          self.z = np.hstack([healpix.z, slice.z])
          self.update_xyz()


def large_scale_phi_covariance(distance):
    # should be something like 
    # cov(r) ~ Int(dk * sin(k*r)/(k**2 * r) )
    # see Equation 9.32 from Dodelson's Cosmology.

    # The integral will diverge unless we put in this k_min.
    k_min = 2.*np.pi / (2. * r_cmb_mpc) # hack
    k_max = 2.*np.pi / (2. * 0.25) # hack
     
    # Evaluate covariance on 1d grid.
    k_vec = np.arange(k_min, k_max, k_min/4.)
    d_vec = np.arange(0., 1.01*distance.max(), 0.1)
    pk_phi = k_vec**(-3.)
    kd_vec = k_vec * d_vec[:,np.newaxis]
    from scipy.special import jv
    cov_vec = np.sum(pk_phi / k_vec * k_vec**3. * jv(0, kd_vec), axis=1)
    #plt.plot(d_vec, cov_vec)

    # Now interpolate onto 2d grid.
    from scipy import interpolate
    f = interpolate.interp1d(d_vec, cov_vec)
    cov = f(distance)

    # Let's force the covariance to be unity along the diagonal.
    # I.e. let's define the variance of each point to be 1.0.
    #cov_diag = cov.diagonal().copy()
    #cov /= np.sqrt(cov_diag)
    #cov /= np.sqrt(cov_diag.T)
    
    return cov



class Inference(object):
     def __init__(self, data_object, test_object):
          # DATA_OBJECT is e.g. a FakeHealpixData object.
          #    It's where you have data.
          # TEST_OBJECT is e.g. a SliceSurface object.
          #     It's where you want to make inferences.
          self.data = data_object
          self.test = test_object

     def calculate_phi_realization(self):
          ###############################################################
          # Coded up from Equation 18 in Roland's note, 
          # https://www.dropbox.com/s/hsq44r7cs1rwkuq/MusicofSphere.pdf
          # Is there a faster algorithm than this?
          ###############################################################
          # Ryan's understanding of this:

          # Define a coordinate object that includes points on the sphere 
          # and on a 2d slice.
          joint = HealpixPlusSlice()
          
          # Do some preparatory work.
          # We only do this once when making multiple realizations.
          if not(hasattr(self, 'cov_joint')):
               dist = joint.make_auto_distance_array()
               cov_joint = large_scale_phi_covariance(dist)
               self.cov_joint = cov_joint
          if not(hasattr(self, 'phi_mv')):
               self.calculate_mv_phi()

          # Generate noise-free truth *simultaneously* on Sphere and Slice.
          from numpy.random import multivariate_normal
          realization_truth = multivariate_normal(np.zeros(joint.n_total), self.cov_joint)
          sphere_truth = realization_truth[joint.ind_healpix]
          slice_truth = realization_truth[joint.ind_slice]

          # Add noise to Sphere points.
          noise = self.data.sigma*np.random.randn(joint.n_healpix)
          sphere_data = sphere_truth + noise

          # Generate MV estimate on Slice.
          tmp = np.dot(self.inv_cov_data_data , sphere_data)
          this_phi_mv_slice = np.dot(self.cov_data_test.T, tmp)

          # Get the difference of the MV estimate on Slice and the truth on Slice.
          diff_mv = this_phi_mv_slice - slice_truth

          # Add that difference to your *original* MV estimate on Slice.
          # Now you have a sample/realization of the posterior on the Slice, given original data.
          this_realization = self.phi_mv + diff_mv

          return this_realization



     def calculate_mv_phi(self):
          self.get_data_data_covariance()
          self.get_data_test_covariance()
          tmp = np.dot(self.inv_cov_data_data , self.data.data)
          self.phi_mv = np.dot(self.cov_data_test.T, tmp)

     def get_data_data_covariance(self):
          # Get phi covariance between data space and data space.
          from numpy.linalg import inv
          dist_data_data = self.data.make_auto_distance_array()
          cov_data_data = large_scale_phi_covariance(dist_data_data)
          self.inv_cov_data_data = inv(cov_data_data)

     def get_data_test_covariance(self):
          # Get phi covariance between data space and test space.
          dist_data_test = self.data.make_distance_array(self.test)
          cov_data_test = large_scale_phi_covariance(dist_data_test)
          self.cov_data_test = cov_data_test

     def view_phi_mv_slice(self):
          self.view_slice(self.phi_mv)

     def view_slice(self, slice_1d):
          slice_2d = slice_1d.reshape(self.test.n_side, self.test.n_side)
          plt.figure(figsize=(7,7))
          plt.imshow(slice_2d, cmap=cmap, vmin=-vscale, vmax=+vscale)

     



     
