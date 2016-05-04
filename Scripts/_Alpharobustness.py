#-----------------------------------------------

# Make a histogram with the alpha values for those numreal reconstructions

# First load the alphas saved
alphas = np.loadtxt("RobustnessAnalysis/alpha_lmax"+str(beatbox.Universe.truncated_lmax)+"lmin"+str(beatbox.Universe.truncated_lmin)+"nmax"+str(beatbox.Universe.truncated_nmax)+"nmin"+str(beatbox.Universe.truncated_nmin)+".txt")

# Make the histogram
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import multivariate_normal

n, bins, patches = plt.hist(alphas, 40, normed=1, facecolor='green', alpha=0.75)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# bug: the Tex interpretor doesn't work properly...
plt.xlabel(r'$\alpha$')
plt.ylabel('Probability')
plt.title(r'Histogram of $\alpha$ values')


# Find the best fit mu and sigma
sigma, mu = np.mgrid[0.001:0.01:.001, 0.9:1.2:.01]
posterior = np.zeros(sigma.shape)

for i in range(sigma.shape[0]):
    for j in range(sigma.shape[1]):
        posterior[i, j] = multivariate_normal.pdf(alphas, mu[i, j]*np.ones(alphas.shape), sigma[i, j]*np.diag(np.ones(alphas.shape)))

a=np.isinf(posterior)
if a.any():
    sigma, mu = np.mgrid[np.std(alphas)**2*0.8:np.std(alphas)**2*1.2:np.std(alphas)**2*0.4/11, 0.999*np.mean(alphas): 1.001*np.mean(alphas): np.mean(alphas)*0.002/11]
    posterior = np.zeros(sigma.shape)
    posterior1 = np.zeros(sigma.shape)
    posterior2 = np.zeros(sigma.shape)
    posterior3 = np.zeros(sigma.shape)
    posterior4 = np.zeros(sigma.shape)
    posterior5 = np.zeros(sigma.shape)
    alphas1=alphas[:np.floor(len(alphas)/5).astype(int)]
    alphas2=alphas[np.floor(len(alphas)/5).astype(int):np.floor(2*len(alphas)/5).astype(int)]
    alphas3=alphas[np.floor(2*len(alphas)/5).astype(int):np.floor(3*len(alphas)/5).astype(int)]
    alphas4=alphas[np.floor(3*len(alphas)/5).astype(int):np.floor(4*len(alphas)/5).astype(int)]
    alphas5=alphas[np.floor(4*len(alphas)/5).astype(int):]

    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            posterior1[i, j] = multivariate_normal.pdf(alphas1, mu[i, j]*np.ones(alphas1.shape), sigma[i, j]*np.diag(np.ones(alphas1.shape)))
    posterior1=posterior1/np.amax(posterior1)
        
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            posterior2[i, j] = multivariate_normal.pdf(alphas2, mu[i, j]*np.ones(alphas2.shape), sigma[i, j]*np.diag(np.ones(alphas2.shape)))
    posterior2=posterior2/np.amax(posterior2)
    
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            posterior3[i, j] = multivariate_normal.pdf(alphas3, mu[i, j]*np.ones(alphas3.shape), sigma[i, j]*np.diag(np.ones(alphas3.shape)))
    posterior3=posterior3/np.amax(posterior3)
    
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            posterior4[i, j] = multivariate_normal.pdf(alphas4, mu[i, j]*np.ones(alphas4.shape), sigma[i, j]*np.diag(np.ones(alphas4.shape)))
    posterior4=posterior4/np.amax(posterior4)
    
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            posterior5[i, j] = multivariate_normal.pdf(alphas5, mu[i, j]*np.ones(alphas5.shape), sigma[i, j]*np.diag(np.ones(alphas5.shape)))
    posterior5=posterior5/np.amax(posterior5)
        
    posterior=posterior1*posterior2*posterior3*posterior4*posterior5
        
max_ind = np.argmax(posterior)



x = np.linspace(0,1.6, 100)
rv = norm()
plt.plot(x, norm.pdf(x, loc=mu.reshape(-1)[max_ind], scale=np.sqrt(sigma.reshape(-1)[max_ind])), 'k-', lw=2, label='frozen pdf')


#plt.axis([40, 160, 0, 0.03])
plt.grid(True)


plt.savefig('RobustnessAnalysis/alphahist_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)+'.png')


plt.show()

#-------------------------------------------

# Calculate the 'recovery success' for the alpha parameter:
sum_alpha = 0
for i in range(len(alphas)):
    sum_alpha = sum_alpha+(alphas[i]-1)**2
    
sum_alpha = sum_alpha/len(alphas)


print 'For l_max = ' + str(beatbox.Multiverse.truncated_lmax) + ' and n_max = ' + str(beatbox.Multiverse.truncated_nmax) + ', the RMS error on the alpha recovery is ' + str(np.sqrt(sum_alpha)) + '.' 