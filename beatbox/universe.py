import numpy as np
import matplotlib.pylab as plt
plt.ion()
import ipdb

r_cmb_mpc = 14.0
cmap = 'gray'

def demo():
     f = FakeHealpixData()
     f.show()
     s = SliceSurface()

     # need to put the below code into objects/methods.

     # 
     autodist = f.make_auto_distance_array()
     autocov = large_scale_phi_covariance(autodist)
     from numpy.linalg import inv
     inv_autocov = inv(autocov)

     crossdist = f.make_distance_array(s)
     crosscov = large_scale_phi_covariance(crossdist)

     w = np.dot(crosscov.T, np.dot(inv_autocov , f.data))
     w_2d = w.reshape(s.n_side, s.n_side)
     plt.clf(); plt.imshow(w_2d, cmap=cmap)


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
        print '...making distance array...'
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


class SphericalSurface(CartesianCoordinates):
    def __init__(self, radius_mpc=r_cmb_mpc, reso_mpc=1.):
        self.epsilon = 0.5*reso_mpc
        pass
        

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
        print '...adding truth...'
        distance = self.make_auto_distance_array()
        delta = distance[distance!=0].min()
        cov = large_scale_phi_covariance(distance)
        from numpy.random import multivariate_normal
        self.data += multivariate_normal(np.zeros(self.n_pix), cov)

    def add_noise(self):
        print '...adding noise...'
        from numpy.random import randn
        self.data += self.sigma*randn(self.n_pix)
        pass

    def show(self):
        from healpy import mollview
        mollview(self.data)


def large_scale_phi_covariance(distance):
    # should be something like 
    # cov(r) ~ Int(dk * sin(k*r)/(k**2 * r) )
    # see Equation 9.32 from Dodelson's Cosmology.

    # The integral will diverge unless we put in this k_min.
    k_min = 2.*np.pi / (2. * r_cmb_mpc) # hack
    k_max = 2.*np.pi / (2. * 0.5) # hack
     
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
