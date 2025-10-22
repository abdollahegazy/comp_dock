"""
Script to see which conformations/posses have the most files to use
from duncans MD simulations
"""
from collections import defaultdict
from pprint import pprint
from pathlib import Path
from tqdm import tqdm
import os,time

LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930,# 2169,514 # omitted bc its always empty ??
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
    '9': 226} 
    # i didnt fully run the other script but it also had 39K out 44K (2000 proteins * 22 ligands) docked poses, so ill use conf 9 pocket 0
    #nvm gonna use conf 3 pocket 0 bc lower pose better and that appears to be good

    on my proteins had 0 issues.

    """

    missing_protein_conf = {}

    for conf in range(3,4):
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
            

# def count_missing_ligands():
#      #tuples are (conf,pose)

#     missing_ligands_conf_poses = defaultdict(int)
#     # confs = [1,2,3,4]

#     st = time.time()
#     s = p = l = 0 
#     for conf in range(3,4):
#         for species in SPECIES:
#             PROTEIN_PATH = Path(f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking")

#             protein_dirs = [p for p in PROTEIN_PATH.iterdir() if p.is_dir()]
#             s+=1
#             for protein_path in tqdm(
#                     protein_dirs, desc=f'iterating {species} on conf {conf}'):
#                     p+=1

#                     for ligand in LIGANDS:
#                         ligand_path = protein_path / str(ligand)
                                         
#                         if not ligand_path.exists():
#                             continue
#                         l+=1

#                         pose_files = [f for f in os.listdir(ligand_path) if f.startswith(f"out_poses_c{conf}_p") and f.endswith(".pdbqt")]

                        # print(pose_files)
                    
                    # for ligand in LIGANDS:
                    #     p = 0


                    #     ligand_path = protein_path / f"{ligand}" 
                    #     ligand_c_p = ligand_path / f"out_poses_c{conf}_p{p}.pdbqt"
                    #     while ligand_c_p.exists():
                                
                    #         if (conf,p) not in missing_ligands_conf_poses:
                    #             missing_ligands_conf_poses[(conf,p)] = 0
                    #         missing_ligands_conf_poses[(conf,p)]+=1
                    #         p+=1
                    #         ligand_c_p = ligand_path / f"out_poses_c{conf}_p{p}.pdbqt"

        # et = time.time()
        # print(missing_ligands_conf_poses)           
        # print(f"Took {et-st:.2f} seconds to iterate over {s} species, {p} proteins, {l} ligands")


    #                 if not protein.exists():
    #                     no_protein_conf+=1
    #                     continue
    #     if conf not in missing_protein_conf:
    #         missing_protein_conf[conf] = no_protein_conf
    
    # # pprint(missing_protein_conf)
    # return missing_protein_conf


# import re 

# def count_missing_ligands_bulk():
#     """Do one filesystem walk instead of 44k operations"""
#     missing_ligands_conf_poses = defaultdict(int)
#     pose_pattern = re.compile(r'out_poses_c(\d+)_p(\d+)\.pdbqt')
    
#     total_processed = 0
    
#     for conf in range(3, 4):
#         for species in SPECIES:
#             base_path = f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking"
            
#             # Single os.walk instead of nested loops
#             print(f"Walking {species}...")
#             # dir_count = sum(1 for _, _, _ in os.walk(base_path))
#             for root, dirs, files in tqdm(os.walk(base_path), desc=f'{species}'):                
#                 # Check if we're in a ligand directory
#                 path_parts = Path(root).parts
#                 if len(path_parts) >= 3:
#                     ligand_name = path_parts[-1]

#                     try:
#                         int(ligand_name)
#                     except:
#                         continue
                    
#                     # print(int(ligand_name))
#                     if int(ligand_name) in LIGANDS:
                        
#                         # Process all pose files in this directory
#                         pose_nums = []
#                         for filename in files:
#                             match = pose_pattern.match(filename)
#                             if match and int(match.group(1)) == conf:
#                                 pose_nums.append(int(match.group(2)))
                        
#                         if pose_nums:
#                             pose_nums.sort()
#                             for i, pose_num in enumerate(pose_nums):
#                                 if pose_num == i:
#                                     missing_ligands_conf_poses[(conf, pose_num)] += 1
#                                 else:
#                                     break
                        
#                         total_processed += 1

#             break
    
#     print(f"Processed {total_processed} ligand directories")
#     return dict(missing_ligands_conf_poses)


#dont mind down here u can remove this if u want cuh u prob want the above 2 funcs

def count_valid_ligands_bulk():
    """Bulk version using os.walk"""
    import os
    
    ligand_counts = {str(k):0 for k in LIGANDS}
    target_file = "out_poses_c3_p0.pdbqt"
    
    for species in SPECIES:
        base_path = f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking"
        
        for root, dirs, files in tqdm(os.walk(base_path)):
            if target_file in files:
                # Extract ligand name from path
                ligand_name = Path(root).name
                if int(ligand_name) in LIGANDS:
                    ligand_counts[ligand_name] += 1

    # Print results
    print(f"\nValid c3_p0 files per ligand:")
    for ligand, count in sorted(ligand_counts.items()):
        print(f"{ligand}: {count}")
    
    return dict(ligand_counts)
    
if __name__ == "__main__":
    # count_missing_protein_confs()
    # count_missing_ligands()

    pprint(count_valid_ligands_bulk())