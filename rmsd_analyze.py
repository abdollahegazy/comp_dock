import json
from tqdm import tqdm
import numpy as np

data = json.load(open("rmsd_conf3.json"))['rmsd']
min_rmsd = 9999999
max_rmsd = 0
min_protein = max_protein = ''

for species in data.keys():

    species_data = data[species]
    
    species_val = []
    
    for protein in tqdm(species_data,
        desc=f'iterating proteins for {species}',
        total=len(species_data)):

        ligand_data = species_data[protein]
        lig_vals = list(ligand_data.values())

        species_val.extend(lig_vals)

        species_data[protein]['mean'] = float(np.mean(lig_vals))
        species_data[protein]['median'] = float(np.median(lig_vals))

    species_data['mean']    = float(np.mean(species_val))
    species_data['median'] = float(np.median(species_val))


for species,species_data in data.items():
    for protein,protein_data in species_data.items():

        if protein in ('mean','median'):
            continue

        protein_avg = protein_data.get('mean')

        if protein_avg is None:
            continue

        if protein_avg < min_rmsd:
            min_rmsd = protein_avg
            min_protein = f'{species}/{protein}'

        if protein_avg > max_rmsd:
            max_rmsd = protein_avg
            max_protein = f'{species}/{protein}'


print(f"minimum avg rmsd protein: {min_protein}, with RMSD {min_rmsd:.2f} Å\n")
print(f"maximum avg rmsd protein: {max_protein}, with RMSD {max_rmsd:.2f} Å")
        # break

with open('rmsd_conf3_with_mean_median.json','w') as f:
    json.dump(data,f,indent=4)
        