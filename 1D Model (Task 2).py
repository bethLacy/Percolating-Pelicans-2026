import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage #measurements deprecated, updated version to use ndimage instead


nsamp = 1000 #set number of simulations
L = 1000 #set number of sites per simulation

#CHANGE ME:
ps = [0.90, 0.92, 0.94, 0.96, 0.98] #set individual site occupancy probabilities to investigate

#initialise arrays to store results
avSizeSims = []
avSizeTheorys = []
charSizeSims = []
charSizeTheorys = []

for p in ps: #repeat the investigation for all the different p values we want to try

    #Cluster number density (for each p):

    #This section is taken from pages 25-26 of the book (with a few updates)
    allarea = np.array([]) #initialise array to store sizes of clusters

    for i in range(nsamp): #repeat for nsamp different simulations
        z = np.random.rand(L) #create array z of L random numbers 0-1
        m = z<p #create array m of whether each site is occupied (a site is occupied if z < p at that position)
        lw, num = ndimage.label(m) #num gives us the number of clusters identified
        labelList = np.arange(1, num + 1) 
        area = ndimage.sum(m, lw, labelList) #area gives us the size of each cluster in this simulation
        allarea = np.append(allarea,area) #append the areas in 'area' to 'allarea', so 'allarea' will contain the sizes of all clusters

    #At this point, we have run nsamp different simulations giving us the main output of allarea
    #allarea contains the size of all different clusters accross all the different simulations

    #this section uses a histogram to divide up the area sizes into bins, so we can plot the density against the interval midpoint
    n,sbins = np.histogram(allarea,bins=int(max(allarea)))
    s = 0.5*(sbins[1:]+sbins[:-1])
    nsp = n/(L*nsamp)

    #this section determines the theoretical cluster number density for each bin based on a typical cluster size
    sxi = -1.0/np.log(p)
    nsptheory = (1-p)**2*np.exp(-s/sxi)

    #Plot of the cluster number density (simulated and theoretical) at each cluster size interval
    plt.figure()
    plt.plot(s, nsp, 'o', label='Simulation')
    plt.plot(s, nsptheory, '-', label='Theory')
    plt.xlabel('$s$ (cluster size)')
    plt.ylabel('$n(s,p)$ (cluster number density)')
    plt.title(f'Cluster number density, p = {p}')
    plt.legend()
    plt.savefig(f'CND_p={p}.eps', format='eps')
    plt.show(block=True)

    #Set up for results output
    print("")
    print(f"p = {p}:")
    print("")


    #Average cluster size (for each p):

    #Calculate results
    avSizeSim = np.sum(allarea**2)/np.sum(allarea)
    avSizeTheory = (1+p)/(1-p)

    #display results
    print(f"Average cluster size (simulated) = {avSizeSim}")
    print(f"Average cluster size (theoretical) = {avSizeTheory}")
    print("")

    #save results for later
    avSizeSims.append(avSizeSim)
    avSizeTheorys.append(avSizeTheory)


    #Characteristic cluster size (for each p):

    #calculate results
    lognsp = np.log(nsp[nsp > 0]) #log non-zero cnds
    gradient, intercept = np.polyfit(s[nsp > 0], lognsp, 1) #fit straight line
    charSizeSim = -1/gradient #the characteristic cluster size is -1 over the gradient we found
    charSizeTheory = sxi #we already found the theoretical ccs earlier as sxi

    #display results
    print(f"Characteristic cluster size (simulated) = {charSizeSim}")
    print(f"Characteristic cluster size (theoretical) = {charSizeTheory}")
    print("")

    #save results for later
    charSizeSims.append(charSizeSim)
    charSizeTheorys.append(charSizeTheory)

    
#Scaling exponents (done only once using all ps results combined):

#Convert to numpy arrays for safety
avSizeSims = np.array(avSizeSims)
avSizeTheorys = np.array(avSizeTheorys)
charSizeSims = np.array(charSizeSims)
charSizeTheorys = np.array(charSizeTheorys)
psArray = np.array(ps)

#set deltaP
deltaP = 1-psArray

#Find sigma exponents (using characteristic cluster size)
charGradientSim, charInterceptSim = np.polyfit(np.log(deltaP), np.log(charSizeSims), 1)
sigmaSim = -1/charGradientSim
charGradientTheory, charInterceptTheory = np.polyfit(np.log(deltaP), np.log(charSizeTheorys), 1)
sigmaTheory = -1/charGradientTheory

#Find gamma exponents (using average cluster size)
avGradientSim, avInterceptSim = np.polyfit(np.log(deltaP), np.log(avSizeSims), 1)
gammaSim = -avGradientSim
avGradientTheory, avInterceptTheory = np.polyfit(np.log(deltaP), np.log(avSizeTheorys), 1)
gammaTheory = -avGradientTheory

#Display results
print("")
print("Scaling exponents:")
print(f"Sigma (simulated) = {sigmaSim}")
print(f"Sigma (theoretical) = {sigmaTheory}")
print(f"Gamma (simulated) = {gammaSim}")
print(f"Gamma (theoretical) = {gammaTheory}")


#Plot average cluster size results:
plt.figure()
plt.plot(ps, avSizeSims, 'r.--', label='Simulation')
plt.plot(ps, avSizeTheorys, 'b.-.', label='Theory')
plt.xlabel('$p$ (probability of individual site occupation)')
plt.ylabel('Average cluster size')
plt.legend()
plt.savefig(f'AverageClusterSize.eps', format='eps')
plt.show(block=True)

#Plot characteristic cluster size results:
plt.figure()
plt.plot(ps, charSizeSims, 'r.--', label='Simulation')
plt.plot(ps, charSizeTheorys, 'b.-.', label='Theory')
plt.xlabel('$p$ (probability of individual site occupation)')
plt.ylabel('Characteristic cluster size')
plt.legend()
plt.savefig(f'CharClusterSize.eps', format='eps')
plt.show(block=True)



#VARIABLES:

#nsamp: number of simulations
#L: number of sites in the simulation (length of the lattice)
#p: probability of a site being occupied

#allarea: array containing the sizes of every cluster found (across all simulations)
#z: array of L random numbers, used to compare to p to decide which sites are occupied
#m: array of booleans of whether or not each site is occupied
#lw: array that assigns a label to each site, giving the same to connected occupied sites, essentially identifying clusters
#num: number of clusters identified in current simulation
#labelList: contains labels of all occupied clusters
#area: array containing size of each cluster in current simulation

#n: array containing the number of clusters of each interval of cluster size across all simulations
#sbins: array containing interval boundaries used to determine n
#s: array containing the midpoints of each interval, used for plotting
#nsp: simulated cluster number density
#sxi: the typical theoretical size of clusters in the current setup
#nsptheory: the theoretical cluster number density

#avSizeSim: the simulated average cluster size
#avSizeTheory: the theoretical average cluster size
#avSizeSims: array of the simulated average cluster sizes (for all ps)
#avSizeTheorys: array of the theoretical average cluster sizes (for all ps)

#lognsp: logs of all the non-zero cnds
#intercept, gradient: constants of the best fit straight line of the lognsps

#charSizeSim: the simulated characteristic cluster size
#charSizeTheory: the theoretical characteristic cluster size
#charSizeSims: array of the simulated characteristic cluster sizes (for all ps)
#charSizeTheorys: array of the theoretical characteristic cluster sizes (for all ps)