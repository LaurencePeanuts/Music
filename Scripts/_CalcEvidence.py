from numpy.linalg import det
ind = np.where(beatbox.Universe.kfilter>0)
PS_long = np.zeros(2*len(ind[1]))
PS_long[:len(ind[1])] = (beatbox.You.all_simulated_universes[0].Power_Spectrum[ind])
PS_long[len(ind[1]):] = (beatbox.You.all_simulated_universes[0].Power_Spectrum[ind])
PS = np.zeros(len(ind[1]))
PS[:len(ind[1])/2] = PS_long[:len(ind[1])/2]
PS[len(ind[1])/2:] = PS_long[len(ind[1]):3*len(ind[1])/2]

PiFact = 1./(2.*np.pi)**(len(MockUniverse[i].ay_real)/2.)
DetDenom = 1./(det(beatbox.You.C_yy)*det(np.diag(PS))*det(beatbox.You.A))**(1./2.)
EvidenceVector = np.zeros(numreal)

for i in range(numreal):

    inv_Cf=np.diag(1./PS)
    #print det(beatbox.You.C_yy), det(np.diag(PS)), 1./det(beatbox.You.A)
    #print 1./(det(beatbox.You.C_yy)*det(np.diag(PS))*det(beatbox.You.A))**(1./2.)
    
    #print np.exp(0.5* np.dot(MockUniverse[i].ay_real.T , np.dot ( np.dot( beatbox.You.inv_Cyy, np.dot( beatbox.You.R_real, np.dot(beatbox.You.inv_A, np.dot(beatbox.You.R_real.T , beatbox.You.inv_Cyy) ) ) ) - beatbox.You.inv_Cyy , MockUniverse[i].ay_real )  ))
    NumExp = np.exp(0.5* np.dot(MockUniverse[i].ay_real.T , np.dot ( np.dot( beatbox.You.inv_Cyy, np.dot( beatbox.You.R_real, np.dot(beatbox.You.inv_A, np.dot(beatbox.You.R_real.T , beatbox.You.inv_Cyy) ) ) ) - beatbox.You.inv_Cyy , MockUniverse[i].ay_real )  ))
    #print 1./(2.*np.pi)**(len(MockUniverse[i].ay_real)/2.)
    EvidenceVector[i] = PiFact * NumExp * DetDenom


n, bins, patches = plt.hist(EvidenceVector, np.logspace(-14, 1, 40), normed=0, facecolor='green', alpha=0.75)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.gca().set_xscale("log")

plt.xlabel(r'Evidence')
plt.ylabel('Probability')
plt.title(r'Histogram of Evidence values for ' + str(numreal) + ' mocks')
# x = np.linspace(-50,50, 100)
# rv = norm()
# plt.plot(x, norm.pdf(x, loc=mu.reshape(-1)[max_ind], scale=np.sqrt(sigma.reshape(-1)[max_ind])), 'k-', lw=2, label='frozen pdf')

# #plt.plot(x, norm.pdf(x, loc=1.027, scale=0.005), 'b-', lw=2, label='frozen pdf')


# #plt.axis([0.5, 4, 0, 2])
plt.grid(True)
path_to_save='RobustnessAnalysis/rob_plt_lmax'+str(beatbox.Universe.truncated_lmax)+'_lmin'+str(beatbox.Universe.truncated_lmin)+'_nmax'+str(beatbox.Universe.truncated_nmax)+'_nmin'+str(beatbox.Universe.truncated_nmin)


try: 
    os.makedirs(path_to_save)
except OSError:
    if not os.path.isdir(path_to_save):
        raise    



plt.savefig(path_to_save+'/Evidence.png')
print np.mean(np.log10(EvidenceVector))
print np.mean((EvidenceVector))

plt.show()