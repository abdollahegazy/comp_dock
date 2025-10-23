from pathlib import Path
from utils import run_command_in_conda, species, mol_list
import re

from rdkit import Chem
from meeko import MoleculePreparation
config = {'keep_nonpolar_hydrogens': False, 'hydrate': False, 'flexible_amides': False, 'macrocycle': False, 'min_ring_size': 7, 'max_ring_size': 33, 'rigidify_bonds_smarts': [], 'rigidify_bonds_indices': [], 'double_bond_penalty': 50, 'atom_type_smarts': {}, 'pH_value': '7.4', 'is_protein_sidechain': False, 'remove_index_map': False, 'remove_smiles': False}
from meeko import obutils
import os

import sys

species = [sys.argv[1]]

pattern = re.compile(r"^protein_conf(\d+)_pose(\d+)$")

print(mol_list)

for s in species:
    num_omit = 0
    pockets = Path(f"../pockets/{s}")

    for i, uid in enumerate(pockets.iterdir()):
        uid = uid.name
        i += 1
        print(f"Processing {s}/{uid}. Complete {i}/500")



        for i,ligand_id in enumerate(mol_list):
            print(f"Processing {ligand_id} for {s}/{uid}. {i+1}/{len(mol_list)} complete.")
            ligand_dir = Path(f"../docking/{s}/{uid}/{ligand_id}")

            if not ligand_dir.exists():
                print(f"Skipping {ligand_dir} (directory not found) for {s}/{uid}")
                continue

            for protein_file in ligand_dir.iterdir():
                # Only process files with `.py` extension
                if protein_file.suffix != '.py':
                    continue
                
                match = pattern.match(protein_file.stem)
                if not match:
                    continue
                
                conf_idx, pose = map(int, match.groups())
                ligand_pdb_path = f"../docking/{s}/{uid}/{ligand_id}/ligand_conf{conf_idx}_pose{pose}.pdb"
                out_pdbqt_path = f"../docking/{s}/{uid}/{ligand_id}/ligand_conf{conf_idx}_pose{pose}.pdbqt"
                pH = 7.4

                if Path(out_pdbqt_path).exists():
                    print(f"Skipping {out_pdbqt_path} (already exists)")
                    continue

                try:
                    frmt = os.path.splitext(ligand_pdb_path)[1][1:]
                    with open(ligand_pdb_path) as f:
                        input_string = f.read()
                    obmol_supplier = obutils.OBMolSupplier(input_string, frmt)
                    for mol in obmol_supplier:
                        preparator = MoleculePreparation.from_config(config)
                        preparator.prepare(mol)
                        ligand_prepared = preparator.write_pdbqt_string()
                except Exception as e:
                    print(f"Skipping {ligand_pdb_path} for {s}/{uid}  (RDKit load error)")
                    continue
                
                try:
                    run_command_in_conda(
                        command=f"mk_prepare_ligand.py -i {ligand_pdb_path} --pH {pH} -o {out_pdbqt_path}",
                        env="vina"
                    )
                except Exception as e:
                    print(f"Error in mk_prepare: species={s}, uid={uid}, conf={conf_idx}, pose={pose}")



