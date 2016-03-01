# First run the _ImportScript.py so that these don't have to be imported:
#import numpy as np
#import math
#import beatbox
#import os.path
#import healpy as hp

# declaring initial objects
beatbox.You.create_original_Universe()
#  Number of noise realizations to make
numreal=10000

MockUniverse = np.array([])
MockUniverse = np.append(MockUniverse, [beatbox.Universe() for i in range(numreal)])
beatbox.You.all_reconstructed_universes = np.append(beatbox.You.all_reconstructed_universes, [beatbox.Universe() for i in range(numreal)])

# Calculate C_yy from the 100 posterior sample Commander Planck CMB temperature maps 
#    or load the C_yy matrix if already calculated
if not os.path.isfile('../data/covCyy_lmax%d_lmin%d.txt' % (beatbox.Multiverse.truncated_lmax, beatbox.Multiverse.truncated_lmin)):
    beatbox.You.read_Planck_samples()
    beatbox.You.calculate_covariance_matrix(filename='lmax%d_lmin%d' % (beatbox.Multiverse.truncated_lmax, beatbox.Multiverse.truncated_lmin))
else:
    beatbox.You.load_covariance_matrix(filename='covCyy_lmax%d_lmin%d.txt' % (beatbox.Multiverse.truncated_lmax, beatbox.Multiverse.truncated_lmin))
    
# Calculate the inverse of the a_y covariance matrix
beatbox.You.calculate_sdv_Cyy_inverse()

# Make numreal simulated Universes to be used as a mocks:
for i in range(numreal):
    beatbox.You.initiate_simulated_universe(printout=0)
    noise = beatbox.You.generate_one_realization_of_noise()
    # Add the noise to the realization of the sky map
    datamap = beatbox.You.all_simulated_universes[-1].ay2ayreal_for_inference(beatbox.You.all_simulated_universes[-1].ay)+noise
    datamap = datamap.T
    
    MockUniverse[i].ayreal2ay_for_mapping(datamap)
    MockUniverse[i].ay2alm(MockUniverse[i].ay)
    #MockUniverse[i].show_CMB_T_map(title="Mock CMB graviational potential     fluctuations as seen from inside the LSS, l_max=%d" % MockUniverse[i].truncated_lmax, from_perspective_of="observer")
    
    # Reconstruct the potential
    beatbox.You.solve_for_3D_potential(datamap, print_alpha = 1)
    
    # Give those reconstructed f_n's to one of the reconstructed Universes
    beatbox.You.all_reconstructed_universes[i].fn = beatbox.You.reconstrunct_fn
#We.fn = beatbox.You.all_simulated_universes[-1].fn * 0.5
    beatbox.You.all_reconstructed_universes[i].transform_3D_potential_into_alm( truncated_nmax=beatbox.You.all_reconstructed_universes[i].truncated_nmax, truncated_nmin=beatbox.You.all_reconstructed_universes[i].truncated_nmin, truncated_lmax=beatbox.You.all_reconstructed_universes[i].truncated_lmax, truncated_lmin=beatbox.You.all_reconstructed_universes[i].truncated_lmin,usedefault=1, fn=1)
    #beatbox.You.all_reconstructed_universes[i].fn.show_CMB_T_map( from_perspective_of="observer")
    
#-----------------------------------------------

# Make a histogram with the alpha values for those numreal reconstructions

# Make the histogram
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import multivariate_normal

n, bins, patches = plt.hist(alphas, 40, normed=1, facecolor='green', alpha=0.75)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# bug: the Tex interpretor doesn't work properly...
plt.xlabel('$\displaystyle\alpha$')
plt.ylabel('Probability')
plt.title('Histogram of $\displaystyle\alpha$ values')


# Find the best fit mu and sigma
sigma, mu = np.mgrid[0.1:0.2:.01, 0.9:1.2:.01]
posterior = np.zeros(sigma.shape)

for i in range(sigma.shape[0]):
    for j in range(sigma.shape[1]):
        posterior[i, j] = multivariate_normal.pdf(alphas, mu[i, j]*np.ones(alphas.shape), sigma[i, j]*np.diag(np.ones(alphas.shape)))

max_ind = np.argmax(posterior)



x = np.linspace(0,1.6, 100)
rv = norm()
plt.plot(x, norm.pdf(x, loc=1.02, scale=0.16), 'k-', lw=2, label='frozen pdf')


#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()

#-------------------------------------------

# Calculate the 'recovery success' for the alpha parameter:
sum_alpha = 0
for i in range(len(alphas)):
    sum_alpha = sum_alpha+(alphas[i]-1)**2
    
sum_alpha = sum_alpha/len(alphas)

print 'The alpha recovery success is ' + sum_alpha + '.' 