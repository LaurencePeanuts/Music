from scipy.stats import chi2
from scipy.special import erfinv
probabilities2=chi2.cdf(chi2vals, len(beatbox.You.all_simulated_universes[-1].fn)) 
sigmas_dev=erfinv(probabilities2)*np.sqrt(2)

import matplotlib.pyplot as plt
n, bins, patches = plt.hist(sigmas_dev,40, facecolor='green', alpha=0.75)
#plt.gca().set_xscale("log")
#n, bins, patches = plt.hist(pvals,bins=np.linspace(0.9, 1.0, 100), facecolor='green', alpha=0.75)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.xlabel(r'$\sigma$')
plt.ylabel('Probability')
plt.title(r'$\sigma$ for 1000 mocks with $l_{max}=8$, $n_{max}=2$')
#plt.xlim(0.95,1)

#x = np.linspace(norm.ppf(0.01, loc=1, scale=0.19), norm.ppf(0.99, loc=1,  scale=0.19), 100)
#x=np.linspace(0,1.6, 100)
#rv = norm()
#plt.plot(x, norm.pdf(x, loc=1.02, scale=0.14), 'k-', lw=2, label='frozen pdf')


#plt.axis([0, 10, 0, 1])
plt.grid(True)

plt.axvline(np.mean(sigmas_dev))
print np.mean(sigmas_dev)
print len(np.squeeze(np.where(sigmas_dev<1.)))*100./len(sigmas_dev)

path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)

try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise  

plt.savefig(path_to_save+'/sigma_histogram.png')

plt.show()

