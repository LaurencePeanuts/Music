import numpy as np
import matplotlib.pylab as plt
plt.ion()
import ipdb


def demo():
     f = FakeHealpixData()
     f.show()

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


class SliceCoordinates(CartesianCoordinates):
    def __init__(self, position=0., side_mpc=25.,reso_mpc=1.):
        self.side_mpc = side_mpc
        self.reso_mpc = reso_mpc
        n_side = int(np.ceil(side_mpc/reso_mpc))
        self.n_side = n_side
        x_2d = np.tile(np.arange(n_side),n_side).reshape(n_side, n_side)
        y_2d = np.tile(np.arange(n_side),n_side).reshape(n_side, n_side).T
        z_2d = np.zeros_like(x_2d) + position
        self.n_side = n_side
        self.x = x_2d.ravel()
        self.y = y_2d.ravel()
        self.z = z_2d.ravel()
        self.update_xyz()


class SphericalSurace(CartesianCoordinates):
    def __init__(self, radius_mpc=14., reso_mpc=1.):
        self.epsilon = 0.5*reso_mpc
        pass
        

class HealpixSphericalSurface(CartesianCoordinates):
    def __init__(self, radius_mpc=14., n_side=2**4):
        # FYI: n_side = 2**4 corresponds to 
        # 0.064 radians resolution = ~0.9 Mpc at z~1000.
        from healpy import nside2npix, pix2vec
        self.n_pix = nside2npix(n_side)
        x, y, z = pix2vec(n_side, range(self.n_pix))
        self.x = radius_mpc*x
        self.y = radius_mpc*y
        self.z = radius_mpc*z
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
    cov = distance**(3.)
    # above is totally wrong!  should be something like 
    # cov(r) ~ Int(dk * sin(k*r)/(k**2 * r) )
    # see Equation 9.32 from Dodelson's Cosmology.
    return cov
