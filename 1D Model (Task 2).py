import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage #measurements deprecated, updated version to use ndimage instead


nsamp = 1000 #set number of simulations
L = 1000 #set number of sites per simulation

#CHANGE ME:
ps = [0.5, 0.75, 0.9] #set individual site occupancy probabilities to investigate

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
    charSizeSim = 0 #TODO: replace placeholder value
    charSizeTheory = sxi

    #display results
    print(f"Characteristic cluster size (simulated) = {charSizeSim}")
    print(f"Characteristic cluster size (theoretical) = {charSizeTheory}")
    print("")

    #save results for later
    charSizeSims.append(charSizeSim)
    charSizeTheorys.append(charSizeTheory)

    
#Scaling exponents (done only once using all ps results combined):

#TODO: Determine the exponents corresponding to the characteristic and average cluster size



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

#charSizeSim: the simulated characteristic cluster size
#charSizeTheory: the theoretical characteristic cluster size
#charSizeSims: array of the simulated characteristic cluster sizes (for all ps)
#charSizeTheorys: array of the theoretical characteristic cluster sizes (for all ps)