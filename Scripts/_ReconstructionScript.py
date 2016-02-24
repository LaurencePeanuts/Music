#%load_ext autoreload
#%autoreload 2

import numpy as np
import math
import beatbox
import os.path
import healpy as hp

np.random.seed(1)

# declaring initial objects
#You=beatbox.Multiverse(truncated_nmax=3, truncated_nmin=1, truncated_lmax=10, truncated_lmin=1)
beatbox.You.create_original_Universe()

#make a realization of the sky to be used as a mock:
beatbox.You.initiate_simulated_universe()

# Calculate C_yy from the 100 posterior sample Commander Planck CMB temperature maps 
#    or load the C_yy matrix if already calculated
if not os.path.isfile('../data/covCyy_lmax%d.txt' % beatbox.Universe.truncated_lmax):
    beatbox.You.read_Planck_samples()
    beatbox.You.calculate_covariance_matrix(filename='lmax%d' % beatbox.Universe.truncated_lmax)
else:
    beatbox.You.load_covariance_matrix(filename='covCyy_lmax%d.txt' % beatbox.Universe.truncated_lmax)
    
# Calculate the inverse of the a_y covariance matrix
beatbox.You.calculate_sdv_Cyy_inverse()

# Fake Cyy matrix
#beatbox.You.inv_Cyy = np.diag(1000*np.ones(beatbox.You.inv_Cyy.shape[1]))
#print "this is a fake Cyy inverse"

if MOCK == 1:

    # Use the simulated Universe as mock data
    # First, generate one realization of noise
    noise = beatbox.You.generate_one_realization_of_noise()
    # Add the noise to the realization of the sky map
    datamap = beatbox.You.all_simulated_universes[-1].ay2ayreal_for_inference(beatbox.You.all_simulated_universes[-1].ay)+noise
    datamap = datamap.T

    # Plot the generated T map
    beatbox.You.all_simulated_universes[-1].show_CMB_T_map(from_perspective_of="observer")

    # Plot the Mock Univers with the noise
    MockUniverse=beatbox.Universe()
    MockUniverse.ayreal2ay_for_mapping(datamap)
    MockUniverse.ay2alm(MockUniverse.ay)
    MockUniverse.show_CMB_T_map( from_perspective_of="observer")
    
    
    
    #_______________________________________________________________

else:
    beatbox.You.all_data_universes = np.append(beatbox.You.all_data_universes, beatbox.Universe())
    beatbox.You.all_data_universes[-1].read_in_CMB_T_map(from_this = 'data/commander_32band_Clsamples100/cmb_Cl_c0001_k00031.fits')
    beatbox.You.all_data_universes[-1].decompose_T_map_into_spherical_harmonics()
    beatbox.You.all_data_universes[-1].alm2ay()
    beatbox.You.all_data_universes[-1].ay2alm(beatbox.You.all_data_universes[-1].ay)
    #beatbox.You.all_data_universes[-1].show_CMB_T_map(from_perspective_of="observer")
    datamap = beatbox.You.all_data_universes[-1].ay2ayreal_for_inference(beatbox.You.all_data_universes[-1].ay)

    # Plot the generated T map
    beatbox.You.all_data_universes[-1].show_CMB_T_map(from_perspective_of="observer")




    #________________________________________________________________





# Reconstruct the potential
beatbox.You.solve_for_3D_potential(datamap)


# Give those reconstructed f_n's to one realization of the Universe
We = beatbox.Universe()
We.fn = beatbox.You.reconstrunct_fn
#We.fn = beatbox.You.all_simulated_universes[-1].fn * 0.5
We.transform_3D_potential_into_alm(truncated_nmax=We.truncated_nmax, truncated_nmin=We.truncated_nmin,truncated_lmax=We.truncated_lmax, truncated_lmin=We.truncated_lmin,usedefault=1, fn=1)
We.show_CMB_T_map(from_perspective_of="observer")

if MOCK == 1:
    # Plot the residuals:
    hp.mollview(MockUniverse.Tmap-We.Tmap,  rot=(-90,0,0),title="CMB graviational potential     fluctuations as seen from inside the LSS, l_max=%d, Tmap diff" % We.truncated_lmax)


    # residuals of the Tmap:
    WeRes = beatbox.Universe()
    WeRes.alm = (beatbox.You.all_simulated_universes[-1].alm-We.alm)

    #WeRes.ay2alm(WeRes.ay)
    WeRes.NSIDE = 256
    WeRes.Tmap = hp.alm2map(WeRes.alm,WeRes.NSIDE)

    hp.mollview(WeRes.Tmap,  rot=(-90,0,0),title="CMB graviational potential fluctuations as seen from inside the LSS, l_max=%d, alms diff" % We.truncated_lmax)
    #WeRes.show_CMB_T_map( from_perspective_of="observer")

else:
    # Plot the residuals:
    hp.mollview(beatbox.You.all_data_universes[-1].Tmap-We.Tmap,  rot=(-90,0,0),title="CMB graviational potential     fluctuations as seen from inside the LSS, l_max=%d, Tmap diff" % We.truncated_lmax)


    # residuals of the Tmap:
    WeRes = beatbox.Universe()
    WeRes.alm = (beatbox.You.all_data_universes[-1].alm-We.alm)

    #WeRes.ay2alm(WeRes.ay)
    WeRes.NSIDE = 256
    WeRes.Tmap = hp.alm2map(WeRes.alm,WeRes.NSIDE)

    hp.mollview(WeRes.Tmap,  rot=(-90,0,0),title="CMB graviational potential fluctuations as seen from inside the LSS, l_max=%d, alms diff" % We.truncated_lmax)
    #WeRes.show_CMB_T_map( from_perspective_of="observer")
