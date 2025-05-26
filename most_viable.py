"""
Script to see which conformations/posses have the most files to use
from duncans MD simulations
"""
from pprint import pprint
from pathlib import Path
from tqdm import tqdm


LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930, # 2169 omitted bc its always empty
]

SPECIES = ["Human", "DouglasFir", "Arabidopsis", "Eucalyptus"]


def count_missing_protein_confs():
    """
    result on duncan proteins:
    {'0': 313,
    '1': 255,
    '10': 235,
    '11': 228,
    '2': 250,
    '3': 250,
    '4': 238,
    '5': 251,
    '6': 237,
    '7': 246,
    '8': 227,
    '9': 226} # i didnt fully run the other script but it also had 39K out 44K (2000 proteins * 22 ligands) docked poses, so ill use conf 9 pose 0

    on my proteins had 0 issues.

    """

    missing_protein_conf = {}

    for conf in range(0,12):
        no_protein_conf = 0
        for species in SPECIES:
            PROTEIN_PATH = Path(f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking")
            # PROTEIN_PATH = Path(f"../docking/docking/{species}")
            for protein_path in tqdm(
                    PROTEIN_PATH.iterdir(),
                    total=len(list(PROTEIN_PATH.iterdir())),
                    desc=f'iterating {species} on conf {conf}'):
                    if protein_path.is_file():
                        continue

                    protein = protein_path / f"protein_c{conf}.pdbqt"

                    if not protein.exists():
                        no_protein_conf+=1
                        continue
        if conf not in missing_protein_conf:
            missing_protein_conf[conf] = no_protein_conf
    
    # pprint(missing_protein_conf)
    return missing_protein_conf
            

def count_missing_ligands():
     #tuples are (conf,pose)

    missing_ligands_conf_poses = {}

    for conf in range(9,10):
        for species in SPECIES:
            PROTEIN_PATH = Path(f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking")
            # PROTEIN_PATH = Path(f"../docking/docking/{species}")
            for protein_path in tqdm(
                    PROTEIN_PATH.iterdir(),
                    total=len(list(PROTEIN_PATH.iterdir())),
                    desc=f'iterating {species} on conf {conf}'):
                    if protein_path.is_file():
                        continue
                    
                    for ligand in LIGANDS:
                        i = p = 0


                        ligand_path = protein_path / f"{ligand}" 
                        ligand_c_p = ligand_path / f"out_poses_c{conf}_p{p}.pdbqt"
                        while ligand_c_p.exists():
                                
                            if (conf,p) not in missing_ligands_conf_poses:
                                missing_ligands_conf_poses[(conf,p)] = 0
                            missing_ligands_conf_poses[(conf,p)]+=1
                            p+=1
                            ligand_c_p = ligand_path / f"out_poses_c{conf}_p{p}.pdbqt"

        break

    print(missing_ligands_conf_poses)           


    #                 if not protein.exists():
    #                     no_protein_conf+=1
    #                     continue
    #     if conf not in missing_protein_conf:
    #         missing_protein_conf[conf] = no_protein_conf
    
    # # pprint(missing_protein_conf)
    # return missing_protein_conf


    
if __name__ == "__main__":


    # count_missing_protein_confs()
    count_missing_ligands()
#obabel -ipdbqt /dogwood/tank/Duncan/IMPACTS2022/Simulations/DouglasFir/docking/A0A0F7C9C2/917/out_poses_c0_p0.pdbqt -opdb -O ligand_pose.pdb -f 1 -l 1
