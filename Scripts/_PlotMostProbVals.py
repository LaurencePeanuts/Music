#Get the noisless reconstruction
execfile ('_MostProbValsBias.py')

#--------------------------------------------------

# Now plot the most probable values

import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import multivariate_normal


num=11

postn=np.zeros(1000)
for i in range(1000):
    postn[i] = beatbox.You.all_reconstructed_universes[i].fn[num] 
    
n, bins, patches = plt.hist(postn, 20, normed=1, facecolor='green', alpha=0.75)

plt.xlabel('$f_n$ values for $n=$'+str(num))
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ most\ probable\ reconstructed\ values:}$')
plt.axis([-1.5, 0., 0, 5])
plt.grid(True)

plt.axvline(beatbox.You.all_simulated_universes[0].fn[num])
plt.axvline(beatbox.You.all_reconstructed_universes[-1].fn[num], color='y')

#N = (2.*pi)**(len(beatbox.You.all_simulated_universes[-1-i].fn)/2.) / (np.linalg.det(beatbox.You.A)**0.5) * np.exp(0.5*np.dot(beatbox.You.all_reconstructed_universes[0].fn.T, np.dot(beatbox.You.inv_A, beatbox.You.all_reconstructed_universes[0].fn)))
#t = np.linspace(-3.,3.,6./200.)
t = np.linspace(-4.5,2,1000)
#plt.plot(t, norm.pdf(t, loc=beatbox.You.all_reconstructed_universes[-1].fn[num], scale=1./np.sqrt(beatbox.You.A[num, num])), 'k-', lw=1, label='frozen pdf')

#plt.plot(t, np.linalg.det(beatbox.You.A)**0.5/(2.*np.pi)**(0.5*len(beatbox.You.all_reconstructed_universes[0].fn)) *np.exp(-0.5*(beatbox.You.all_reconstructed_universes[0].fn[num]-t)**2. * beatbox.You.A[num, num]), 'k-', lw=2, label='frozen pdf')


#xj=np.append(beatbox.You.all_reconstructed_universes[0].fn[:num], beatbox.You.all_reconstructed_universes[0].fn[num+1:])
#Ani=np.append(beatbox.You.A[num, :num], beatbox.You.A[num, 1+num:])
#Aij_1112 = np.append(beatbox.You.A[:num, :num], beatbox.You.A[:num, 1+num:], axis=1)
#Aij_2122 = np.append(beatbox.You.A[1+num:, :num], beatbox.You.A[1+num:, 1+num:], axis=1)
#Aij = np.append(Aij_1112, Aij_2122, axis=0)

#first = (t - 1./beatbox.You.A[num, num]*(beatbox.You.A[num, num]*beatbox.You.all_reconstructed_universes[0].fn[num]-np.dot(Ani,xj)))**2. * beatbox.You.A[num, num]
#second =  beatbox.You.all_reconstructed_universes[0].fn[num]**2 * beatbox.You.A[num, num]
#third = (np.dot(Ani, xj))**2 * 1./beatbox.You.A[num, num]
#fourth = np.dot(xj.T , np.dot(Aij,xj))

#plt.plot(t, np.linalg.det(beatbox.You.A)**0.5/(2.*np.pi)**(0.5*len(beatbox.You.all_reconstructed_universes[0].fn)) *np.exp(-0.5 * (first-second+third+fourth )), 'r-', lw=1)




plt.plot(t, norm.pdf(t, loc=beatbox.You.all_reconstructed_universes[-1].fn[num], scale=np.sqrt(beatbox.You.inv_A[num, num])), 'r-', lw=1, label='frozen pdf')

path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)

try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise    

plt.savefig(path_to_save +'/mostprobvalues_f11.png')




plt.show()