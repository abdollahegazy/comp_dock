import os

faulty_dir = "ISSUES"
fixed_dir = "FIXED"

for species in os.listdir(faulty_dir):
    species_path = os.path.join(faulty_dir,species)
    fixed_species_path = fixed_dir + "/"+ species
    os.makedirs(fixed_species_path,exist_ok=True)
    for protein in os.listdir(species_path):
        faulty_protein_dir = os.path.join(species_path,protein)

        with open(fixed_species_path + "/" + protein,'w') as f:
            with open(faulty_protein_dir) as g:
                for line in g:
                    if "UNK" in line:
                        continue
                    f.write(line)
