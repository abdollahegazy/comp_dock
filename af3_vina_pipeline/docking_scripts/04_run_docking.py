from pathlib import Path
from utils import run_command_in_conda,species,mol_list
import re

pattern = re.compile(r"^protein_conf(\d+)_pose(\d+)$")


import sys

species = [sys.argv[1]]


error_species = set()

for s in species:
    pockets = Path(f"../pockets/{s}")

    for i,uid in enumerate(pockets.iterdir()):
        uid = uid.name
        i+=1
        print(f"Processing {s}/{uid}. Complete {i}/500")


        for ligand_id in mol_list:
            ligand_dir = Path(f"../docking/{s}/{uid}/{ligand_id}")


            if not ligand_dir.exists():
                print(f"Skipping {ligand_dir} (directory not found) for {s}/{uid}")
                continue
            
            for protein_file in ligand_dir.iterdir():

                #bc technically itll iter over all file sbut
                #that will be doulbe counting
                if protein_file.suffix != '.py':
                    continue

                print(f"Docking for {s}/{uid}")
                
                match = pattern.match(protein_file.stem)
                conf,pose = map(int,match.groups())

                out_file = f"../docking/{s}/{uid}/{ligand_id}/out_conf{conf}_pose{pose}.log"

                # checker = f"/tank/abdolla/docking/{s}/{uid}/{ligand_id}/out_poses_c{conf}_p{pose}.pdbqt"
                
                # if Path(checker).exists():
                #     print(f"Skipping {checker} (already exists)")
                #     continue
                
                pdbqt_file = f"/tank/abdolla/docking/{s}/{uid}/{ligand_id}/ligand_conf{conf}_pose{pose}.pdbqt"

                if not Path(pdbqt_file).exists():
                    print(f"Skipping {out_file} (pdbqt file does not exist)")
                    continue


                try:
                    run_command_in_conda(
                        f"python3 {protein_file} > {out_file}",
                        env="vina"
                    )
                except:
                    error_species.add(f"{s}/{uid}/{ligand_id}/protein_conf{conf}_pocket{pose}")


print("THESE WERE THE ERRORS")
print(error_species)

