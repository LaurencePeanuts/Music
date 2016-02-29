# First run the _ImportScript.py so that these don't have to be imported:
#import numpy as np
#import math
#import beatbox
#import os.path
#import healpy as hp

# declaring initial objects
beatbox.You.create_original_Universe()
#  Number of noise realizations to make
numreal=1000

#Make the mock Universe
beatbox.You.all_simulated_universes[0].show_CMB_T_map(from_perspective_of="observer")


beatbox.You.all_reconstructed_universes = np.append(beatbox.You.all_reconstructed_universes, [beatbox.Universe() for i in range(numreal)])

for i in range(numreal):
    # First, generate one realization of noise
    noise = beatbox.You.generate_one_realization_of_noise()
    # Add the noise to the initial realization of the sky map
    datamap = beatbox.You.all_simulated_universes[0].ay2ayreal_for_inference(beatbox.You.all_simulated_universes[0].ay)+noise
    datamap = datamap.T
    
    # Reconstruct the fn's
    beatbox.You.solve_for_3D_potential(datamap)
    beatbox.You.all_reconstructed_universes[i].fn=beatbox.You.reconstrunct_fn
    beatbox.You.all_reconstructed_universes[i].transform_3D_potential_into_alm( usedefault=1, fn=1)
    

# ---------------------------------------

# Now look at the distribution of reconstructed f_n's
# Which f_n to look at in each reconstruction:    
n = 30

smarts=np.zeros(numreal)
for i in range(numreal):
    smarts[i] = beatbox.You.all_reconstructed_universes[i].fn[n] 
print smarts[:].shape

# Make a histogram
nbins, bins, patches = plt.hist(smarts, 100, normed=1, facecolor='green', alpha=0.75)
plt.xlabel('Value of f_n for n = '+n)
plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.axvline(beatbox.You.all_simulated_universes[0].fn[n])

#plt.axvline(We.fn[5])

plt.show()

# -------------------------------------------

# Now look at the distribution of reconstructed a_y's

# what value of y to look at
yval = 74
ayval=np.zeros(numreal)
for i in range(numreal):
    ayval[i] = beatbox.You.all_reconstructed_universes[i].ay[yval][0]


n, bins, patches = plt.hist(ayval, 40, normed=1, facecolor='yellow', alpha=0.75)

plt.xlabel('Value of a_y for y = ' + yval)
plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.axvline(beatbox.You.all_simulated_universes[0].ay[yval])

#plt.axvline(We.ay[5])

plt.show()

# -----------------------------------------------

#Sample the posterior distributon

beatbox.You.generate_realizations_from_posterior(We.fn, number_of_realizations=1000)

for k in range(10):
    beatbox.You.all_simulated_universes[-1-k].show_CMB_T_map(from_perspective_of =  "observer")
    
#for k in range(10):
#    beatbox.You.all_simulated_universes[-1-k].show_CMB_T_map(from_perspective_of = s"observer")
post106=np.zeros(1000)
for i in range(1000):
    post106[i] = beatbox.You.all_simulated_universes[-1-i].fn[106] 
    
n, bins, patches = plt.hist(post106, 20, normed=1, facecolor='green', alpha=0.75)

plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.axvline(beatbox.You.all_simulated_universes[0].fn[106])

#plt.axvline(We.fn[5])

plt.show()