sum_fn = np.zeros(numreal)
for n in range(numreal):
    for i in range(len(beatbox.You.all_simulated_universes[-1].fn)):
        sum_fn[n] = sum_fn[n]+(beatbox.You.all_simulated_universes[-1].fn[i]-beatbox.You.all_reconstructed_universes[n].fn[i])**2
    
    sum_fn[n] = sum_fn[n]/len(beatbox.You.all_simulated_universes[-1].fn)

n, bins, patches = plt.hist(sum_fn, 50, normed=0, facecolor='green', alpha=0.75)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.xlabel('FoM')
# plt.ylabel('Probability')
plt.title(r'Histogram of values of the Figure of Merit for '+str(numreal)+' realizations of noise')


# #plt.axis([0.5, 4, 0, 2])
plt.grid(True)
path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)

try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise    



plt.savefig(path_to_save+'/figureofmerit.png')



plt.show()
print np.mean(sum_fn)