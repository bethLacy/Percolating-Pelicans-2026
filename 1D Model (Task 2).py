import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage #measurements deprecated, updated version to use ndimage instead

nsamp = 1000
L = 1000
p = 0.90

allarea = np.array([])

for i in range(nsamp):
    z = np.random.rand(L)
    m = z<p
    lw, num = ndimage.label(m)
    labelList = np.arange(1, num + 1)
    area = ndimage.sum(m, lw, labelList)
    allarea = np.append(allarea,area) 

n,sbins = np.histogram(allarea,bins=int(max(allarea)))
s = 0.5*(sbins[1:]+sbins[:-1])
nsp = n/(L*nsamp)
sxi = -1.0/np.log(p)
nsptheory = (1-p)**2*np.exp(-s/sxi)

plt.figure()

plt.plot(s, nsp, 'o', label='Simulation')
plt.plot(s, nsptheory, '-', label='Theory')

plt.xlabel('$s$')
plt.ylabel('$n(s,p)$')
plt.legend()

plt.show(block=True)