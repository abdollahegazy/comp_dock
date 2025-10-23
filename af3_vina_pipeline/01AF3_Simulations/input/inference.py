import os
import json
from chemspipy import ChemSpider
from pprint import pprint
from tqdm import tqdm  
from copy import deepcopy

cs = ChemSpider('TzfV80skfR1yCAe7oY4y06Q2s5Kzlff4a1NTNVq2')

OUT_PATH = "/tank/abdolla/project/dockingAF3"
os.makedirs(OUT_PATH, exist_ok=True)

SPECIES = ["Arabidopsis", "DouglasFir", "Eucalyptus", "Human"]

LIGANDS = [6309, 780, 4444606, 6085, 2169, 23208, 4444418, 6223, 3698, 395716,
           30999, 21105998, 388386, 13876103, 917, 2424, 4450907, 7930, 10194105, 84989, 108426, 393012,
           4518347, 514]

LIGAND_SMILES = {}

print("Fetching SMILES strings...")
for ligand in tqdm(LIGANDS, desc="Ligands"):
    LIGAND_SMILES[ligand] = cs.get_details(ligand)['smiles']

print("=" * 10)
print("LIGAND INFORMATION")
pprint(LIGAND_SMILES)
print("=" * 10)

for s in tqdm(SPECIES, desc="Species"):
    print(f"\nWORKING ON SPECIES {s}")
    species_out = os.path.join(OUT_PATH, s)
    os.makedirs(species_out, exist_ok=True)

    msa_path = os.path.join(s, "MSA")
    msa_files = os.listdir(msa_path)

    for protein_MSA_json_file in tqdm(msa_files, desc=f"{s} proteins", leave=True):
        protein_MSA_json_path = os.path.join(msa_path, protein_MSA_json_file)
        protein_id = protein_MSA_json_file.split("_")[0]
        protein_out = os.path.join(species_out, protein_id)
        os.makedirs(protein_out, exist_ok=True)

        protein_af3_json = json.load(open(protein_MSA_json_path))

        for ligand in tqdm(LIGANDS, desc=f"Ligands for {protein_id}", leave=True):
            protein_af3_json_temp = deepcopy(protein_af3_json)
            ligand_out = os.path.join(protein_out, f"{ligand}.json")
            
            protein_af3_json_temp['name'] = f"{ligand}"

            protein_af3_json_temp['sequences'].append(
                {
                    "ligand": {
                        "id": ["B"],
                        "smiles": LIGAND_SMILES[ligand]
                    }
                }
            )

            with open(ligand_out, "w") as f:
                json.dump(protein_af3_json_temp, f, indent=2)
