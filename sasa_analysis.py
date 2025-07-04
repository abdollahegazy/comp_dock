import freesasa
freesasa.setVerbosity(freesasa.nowarnings)

import json
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

sasas  = json.load(open('sasa_conf3_md_proteins.json'))
rmsd = json.load(open('rmsd_conf3_with_mean_median.json'))
protein_rmsds = json.load(open('protein_rmsd_conf3.json'))['protein_rmsd']

avg_sasas = {species: sum(d.values()) / len(d) for species, d in sasas.items()}

print(avg_sasas)
# quit()
species_list = []
x = []
y = []
for species,protein_data in sasas.items():
    
    for protein_name,sasa in protein_data.items():
        try:
            sasa = sasas[species][protein_name]
            avg_rmsd = rmsd[species][protein_name]['median']
            protein_rmsd = protein_rmsds[species][protein_name]
            x.append(sasa)
            y.append(avg_rmsd)
            species_list.append(species)
        except:
            continue
        # print(protein_name,sasa,avg_rmsd)



slope, intercept = np.polyfit(x, y, 1)
y_pred = np.array(x) * slope + intercept
# Compute R^2 (coefficient of determination)
ss_res = np.sum((np.array(y) - y_pred) ** 2)
ss_tot = np.sum((np.array(y) - np.mean(y)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

species_names = sorted(list(set(species_list)))
colors = plt.cm.tab10(range(len(species_names)))
species_to_color = {s: c for s, c in zip(species_names, colors)}

# Scatter each species separately for legend
plt.figure()
for s in species_names:
    xs = [x_i for x_i, sp in zip(x, species_list) if sp == s]
    ys = [y_i for y_i, sp in zip(y, species_list) if sp == s]
    plt.scatter(xs, ys, label=s, color=species_to_color[s])


print(len(x))
x_min, x_max = min(x), max(x)
x_fit = np.linspace(x_min, x_max, 100)
y_fit = slope * x_fit + intercept
plt.plot(x_fit, y_fit, color='black', linestyle='--', label='Best fit')

plt.xlabel("Protein SASA")
plt.ylabel("Avg Ligand RMSD")
plt.title(f"Protein SASA vs. Ligand RMSD, colored by species (R = {r_squared**.5:.3f})")
plt.legend()
plt.savefig("protein_md_sasa_vs_ligand_rmsd_median.png")
plt.show()

# plt.xlabel("Protein SASA")
# plt.ylabel("Avg Ligand RMSD")
# plt.title(f"SASA vs. Ligand RMSD, colored by species (R² = {r_squared:.3f})")
# plt.legend()
# plt.savefig("test_species_colored.png")
# plt.show()

# xy_sorted = sorted(zip(x, y), key=lambda pair: pair[0])
# x_sorted, y_sorted = zip(*xy_sorted)

# plt.scatter(x_sorted, y_sorted)
# # plt.savefig("test")

# plt.scatter(x,y)
# plt.savefig("test")