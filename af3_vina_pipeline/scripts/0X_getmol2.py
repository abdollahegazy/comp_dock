from chemspipy import ChemSpider
import subprocess
import os

# ChemSpider API key
cs = ChemSpider('TzfV80skfR1yCAe7oY4y06Q2s5Kzlff4a1NTNVq2')

# Ensure output directory exists
output_dir = "../SmallMolecules/ligand_mol2"
os.makedirs(output_dir, exist_ok=True)

# Temp directory for MOL files
mol_dir = "../SmallMolecules/temp_mol"
os.makedirs(mol_dir, exist_ok=True)

# List of compound IDs
compound_ids = [
    6309, 780, 4444606, 6085, 2169, 23208, 4444418, 6223, 3698, 395716,
    30999, 21105998, 388386, 13876103, 917, 2424, 4450907, 7930, 10194105,
    84989, 108426, 393012, 4518347, 514
]

# Loop through each compound
for cid in compound_ids:
    mol2_file = os.path.join(output_dir, f"{cid}.mol2")
    mol_file = os.path.join(mol_dir, f"{cid}.mol")

    # Skip if already processed
    if os.path.exists(mol2_file):
        continue

    try:
        print(f"Processing CID {cid}...")

        # Get ChemSpider record and download MOL file
        record = cs.get_details(cid)
        mol_data = record['mol3D']

        with open(mol_file, 'w') as f:
            f.write(mol_data)

        # Convert MOL to MOL2 using Open Babel
        cmd = f'obabel "{mol_file}" -O "{mol2_file}"'
        subprocess.run(cmd, shell=True, check=True)

    except Exception as e:
        print(f"Failed to process CID {cid}: {e}")
