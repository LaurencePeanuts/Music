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

#print (beatbox.You.all_reconstructed_universes[-1].fn-beatbox.You.all_simulated_universes[0].fn)[:(len(beatbox.You.all_reconstructed_universes[-1].fn)/2)]
#print (beatbox.You.all_reconstructed_universes[-1].fn-beatbox.You.all_simulated_universes[0].fn)[(len(beatbox.You.all_reconstructed_universes[-1].fn)/2):]

#Calculate the bias vectors
realbias = (beatbox.You.all_reconstructed_universes[-1].fn - beatbox.You.all_simulated_universes[0].fn)[:(len(beatbox.You.all_reconstructed_universes[-1].fn)/2)]
imagbias = (beatbox.You.all_reconstructed_universes[-1].fn  -beatbox.You.all_simulated_universes[0].fn)[(len(beatbox.You.all_reconstructed_universes[-1].fn)/2):]

#Print the bias stats
print np.mean(realbias), np.std(realbias), np.max(realbias), np.min(realbias)
print np.mean(imagbias), np.std(imagbias), np.max(imagbias), np.min(imagbias)



#Save a histogram of the real & imaginary biases
n, bins, patches = plt.hist(realbias, 80, normed=0, facecolor='green', alpha=0.75)

n2, bins2, patches2 = plt.hist(imagbias, 10, normed=0, facecolor='yellow', alpha=0.75)


path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)

try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise    



plt.savefig(path_to_save+'/bias_histogram.png')




plt.show()