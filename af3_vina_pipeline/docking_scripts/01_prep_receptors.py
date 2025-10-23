from pathlib import Path
from utils import run_command_in_conda,species



def prepare_receptor(s,uid,conf):
    receptor_file = f"../docking/{s}/{uid}/protein_conf{conf}.pdbqt"
    input_file = f"../Data/{s}/{uid}/protein_conf{conf}.pdb"

    #this is bc i had to find the py2->3 ported version of ADFR
    adfr_prepare_receptor_file = "prepare_receptor4"
    
    run_command_in_conda(
    command=f"{adfr_prepare_receptor_file} -r {input_file} -o {receptor_file}",
    env="adfr"
    )

for s in species:
    pockets = Path(f"../pockets/{s}")

    for i,uid in enumerate(pockets.iterdir()):
        uid = uid.name
        i+=1
        print(f"Processing {s}/{uid}. Complete {i}/500")

        docking_dir = Path(f"../docking/{s}/{uid}")


        for conf in range(0,12):
            center_file = Path(f"../pockets/{s}/{uid}/protein_conf{conf}/centers.txt")
            pdbqt_file = Path(f"../docking/{s}/{uid}/protein_conf{conf}.pdbqt")

            if center_file and not pdbqt_file.exists():
                docking_dir.mkdir(parents=True,exist_ok=True)
                prepare_receptor(s,uid,conf)



