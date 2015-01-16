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
     for i in range(10):
          this_realization = inf.calculate_phi_realization()


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
    def __init__(self, position=0., side_mpc=35., reso_mpc=0.5):
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
          # Need to code up Equation 18 in Roland's note, 
          # https://www.dropbox.com/s/hsq44r7cs1rwkuq/MusicofSphere.pdf
          # Ryan's understanding of this:
          # Generate noise-free truth *simultaneously* on Sphere and Slice.
          # Add noise to Sphere points.
          # Generate MV estimate on Slice.
          # Get the difference of the MV estimate on Slice and the truth on Slice.
          # Add that difference to your *original* MV estimate on Slice.
          # Now you have a sample/realization of the posterior on the Slice, given original data.
          pass

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

     



     
