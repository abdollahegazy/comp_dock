import numpy as np
import scipy.stats as stats
import matplotlib as plt
import random
import glob
import os


#rng = np.random.default_rng()
#print(stats.kstest(stats.uniform.rvs(size=100, random_state=rng), stats.norm.cdf))



#define chemid
bindingEnergies = np.load("../Data/bindingEnergies.npz", allow_pickle=True)
speciesList = bindingEnergies['dim1']
molNums = bindingEnergies['dim2']
results = bindingEnergies['data']



#randomList = [] #would be best to make this into an array

for i in range(0,500):
	n = random.random()
	randomList.append(n)
	print(i, n)


controlEnergies = np.empty(1,)
species = ['Arabidopsis','DouglasFir','Eucalyptus','Human']

def analyse_set(s,data):
for s in species:
	results = np.empty(shape = (len(species)), dtype = object)
	uniprot = np.empty(shape = (len(species)), dtype = object)

for i in range(len(uniprot)):
    uniprot[i] = []
    for d in glob.glob('../Simulations/%s/docking/*' % species[i]):
        uid = os.path.basename(d)
        uniprot[i].append(uid)
    results[i] = analyse_set(species[i],uniprot[i])    
    

	#define uid
#	for u in uid:
#		for c in chemid:
#			stats.kstest(stats.uniform.rvs(randomList,bindingEnergies)
			#save this into a set for use in graph?
