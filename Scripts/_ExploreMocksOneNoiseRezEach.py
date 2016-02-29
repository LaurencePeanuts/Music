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
    beatbox.You.initiate_simulated_universe()
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