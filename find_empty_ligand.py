from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from pprint import pprint
from tqdm import tqdm
import os

SPECIES = ["Human", "DouglasFir", "Arabidopsis", "Eucalyptus"]

def count_ligands_for_species(species):
    MD_PATH = Path(f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking/")
    ligand_counts = Counter()

    try:
        for protein in tqdm(MD_PATH.iterdir(),desc=f'processing {species}', leave=False,total=len(list(MD_PATH.iterdir()))):
            for ligand in protein.iterdir():
                if ligand.is_file():
                    continue

                if ligand.name not in ligand_counts:
                    ligand_counts[ligand.name] = 0

                for f in ligand.iterdir():
                    if f.suffix == '.pdbqt':
                        ligand_counts[ligand.name] += 1
    except Exception as e:
        print(f"Error processing {species}: {e}")
    
    return ligand_counts


if __name__ == "__main__":
    all_counts = Counter()

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = {executor.submit(count_ligands_for_species, species): species for species in SPECIES}
        for future in tqdm(as_completed(futures), total=len(futures), desc="processing species"):
            counts = future.result()
            all_counts.update(counts)

    pprint(dict(all_counts))

# output gives
# {'10194105': 75158,
#  '108426': 75156,
#  '13876103': 75184,
#  '21105998': 75184,
#  '2169': 0,
#  '23208': 75186,
#  '2424': 75178,
#  '30999': 75164,
#  '3698': 75170,
#  '388386': 75184,
#  '393012': 75126,
#  '395716': 75176,
#  '4444418': 75162,
#  '4444606': 75177,
#  '4450907': 75158,
#  '4518347': 75130,
#  '6085': 75186,
#  '6223': 75148,
#  '6309': 75189,
#  '780': 75181,
#  '7930': 74969,
#  '84989': 75130,
#  '917': 75180}
#so 2169 is the bad ligand
