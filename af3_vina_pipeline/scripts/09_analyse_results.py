import pickle
import json
import shutil
import numpy as np
import os
import glob


species = ['Arabidopsis','DouglasFir','Eucalyptus','Human']

molList = os.listdir('../SmallMolecules/molfiles')
molNums = [mol.split('.')[0] for mol in molList]
print(molNums)

#def analyse_set(s,data,mol):
def analyse_set(s,cmpd,data):
    print('protein\t\tAlphafold\tPockets\conf\tEnergy')
    results = []
    for uid in data:
        s_alpha = '-'
        if os.path.exists('../Data/%s/alphafold_%s.pkl' % (s, uid)):
            with open('../Data/%s/alphafold_%s.pkl' % (s,uid), 'rb') as f:
                dataPkl = pickle.load(f)
            alpha_plddt=dataPkl['plddt']
            s_alpha = '%.2f' % (np.average(alpha_plddt)/100)

        # s_rosetta = '-'
        # if os.path.exists('../Data/%s/rosetta_%s.npz' % (s,uid)): #remove
            # data = np.load('../Data/%s/rosetta_%s.npz' % (s,uid)) #
            # rosetta_plddt=data['lddt'] ############################
            # s_rosetta = '%.2f' % (np.average(rosetta_plddt))#######
        

        files = glob.glob('../Simulations/%s/docking/%s/%s/out_c*_p*.log' % (s,uid,cmpd))
        energies = []
        for f in files:
            fIn = open(f, 'r')
            for line in fIn:
                word = line.split()
                if len(word)>1:
                    if word[0]=='1':
                        energies.append(float(word[1]))
                        break
            fIn.close()
        if len(energies)>0:
            results.append(np.min(energies))
            print('%-10s\t%s\t\t%.2f\t\t%.2f' % (uid, s_alpha, len(files)/12, np.min(energies)))
        else:
            print('%-10s\t%s\t\t%.2f\t\t  -' % (uid, s_alpha, len(files)/12))
            print(s, cmpd)
            results.append(np.inf)
    return np.array(results)

uniprot = np.empty(shape = (len(species),len(molNums)), dtype = np.object)
results = np.empty(shape = (len(species),len(molNums)), dtype = np.object)
#molecule = np.empty(shape = (len(species)), dtype = np.object)


for i in range(len(species)):
    for j in range(len(molNums)):
        uniprot[i][j] = []
        for d in glob.glob('../Simulations/%s/docking/*' % species[i]):
            uid = os.path.basename(d)
            uniprot[i][j].append(uid)
        results[i][j] = analyse_set(species[i],molNums[j],uniprot[i][j])
        #results[i] = analyse_set(species[i],uniprot[i],molNums)    
    

np.savez('../Data/bindingEnergies.npz', data=results,dim1=species,dim2=molNums)