beatbox.You.generate_realizations_from_posterior(beatbox.You.all_reconstructed_universes[0].fn, number_of_realizations=1000)

#for k in range(10):
#    beatbox.You.all_simulated_universes[-1-k].show_CMB_T_map(from_perspective_of =  "observer")
   
#for k in range(10):
#    beatbox.You.all_simulated_universes[-1-k].show_CMB_T_map(from_perspective_of = s"observer")
post106=np.zeros(1000)
for i in range(1000):
    post106[i] = beatbox.You.all_simulated_universes[-1-i].fn[num]

   
n, bins, patches = plt.hist(post106, 20, normed=1, facecolor='green', alpha=0.75)

#plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.axvline(beatbox.You.all_simulated_universes[2].fn[num])
#plt.axvline(beatbox.You.all_reconstructed_universes[-1].fn[num])

#n2, bins2, patches2 = plt.hist(postn, 100, normed=1, facecolor='yellow', alpha=0.75)

#plt.axvline(We.fn[5])

t = np.linspace(-5,5,1000)
plt.plot(t, norm.pdf(t, loc=beatbox.You.all_reconstructed_universes[0].fn[num], scale=np.sqrt(beatbox.You.inv_A[num, num])), 'r-', lw=2, label='frozen pdf')
plt.plot(t, 1./(np.sqrt(beatbox.You.inv_A[num, num]*2.*np.pi))*np.exp(-0.5*(t-beatbox.You.all_reconstructed_universes[0].fn[num])**2./beatbox.You.inv_A[num, num]), 'k:', lw=2)

path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)

try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise    

plt.savefig(path_to_save+'/posterior_f11.png')



plt.show()