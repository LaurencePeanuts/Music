#Get the noisless reconstruction
#Create an empty Universe
beatbox.You.all_reconstructed_universes= np.append(beatbox.You.all_reconstructed_universes, beatbox.Universe())
#Reconstruct the f_n form the noisless a_lm and put them in the empty universe
beatbox.You.solve_for_3D_potential(beatbox.You.all_simulated_universes[0].ay2ayreal_for_inference(beatbox.You.all_simulated_universes[0].ay), print_alpha = 0)
beatbox.You.all_reconstructed_universes[-1].fn=beatbox.You.reconstrunct_fn
beatbox.You.all_reconstructed_universes[-1].transform_3D_potential_into_alm( usedefault=1, fn=1)
# Make the phi map for that Universe
beatbox.You.all_reconstructed_universes[-1].rearrange_fn_from_vector_to_grid()
beatbox.You.all_reconstructed_universes[-1].evaluate_potential_given_fourier_coefficients()

print (beatbox.You.all_reconstructed_universes[-1].fn-beatbox.You.all_simulated_universes[0].fn)[:(len(beatbox.You.all_reconstructed_universes[-1].fn)/2)]
print (beatbox.You.all_reconstructed_universes[-1].fn-beatbox.You.all_simulated_universes[0].fn)[(len(beatbox.You.all_reconstructed_universes[-1].fn)/2):]

#Calculate the bias vectors
realbias = (beatbox.You.all_reconstructed_universes[-1].fn - beatbox.You.all_simulated_universes[0].fn)[:(len(beatbox.You.all_reconstructed_universes[-1].fn)/2)]
imagbias = (beatbox.You.all_reconstructed_universes[-1].fn  -beatbox.You.all_simulated_universes[0].fn)[(len(beatbox.You.all_reconstructed_universes[-1].fn)/2):]

#Print the bias stats
print np.mean(realbias), np.std(realbias), np.max(realbias), np.min(realbias)
print np.mean(imagbias), np.std(imagbias), np.max(imagbias), np.min(imagbias)