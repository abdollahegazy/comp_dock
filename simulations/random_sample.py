import random
from pathlib import Path
from collections import defaultdict
from pprint import pprint

random.seed(420)

AF_COMPLEX = Path("../MD_complexes_c3")
MD_COMPLEXES = Path("")
k = 20

def get_all_protein_paths():
    """
    just gets a list of every protein. important if u want to sample proteins first then ligands (decrease chance of just
    using multiple ligands from a single protein)
    """
    all_proteins_paths = []

    for species in AF_COMPLEX.iterdir():

        for protein in species.iterdir():
            all_proteins_paths.append(protein)

    return all_proteins_paths

def analyze_sample(sample):
    """
    analyzes sample by seeing breakdown of species and ligands
    TODO: also incorporate base count to ensure no bias
    """
    species_counts = defaultdict(int)
    protein_counts = defaultdict(int)
    ligand_counts = defaultdict(int)

    for observation in sample:
        name_parts = observation.parts
        species = name_parts[2]
        protein = name_parts[3]
        ligand = name_parts[4]

        species_counts[species] += 1
        protein_counts[protein] += 1
        ligand_counts[ligand] += 1
    
    analysis =   {
        "species": dict(species_counts),
        "proteins": dict(protein_counts),
        "ligands": dict(ligand_counts),
    }
    pprint(analysis)
def two_stage_sampling():
    sampled_complexes = []

    all_proteins_paths = get_all_protein_paths()

    k_random_protein_paths = random.sample(all_proteins_paths,k=k)

    for proteins in k_random_protein_paths:
        curr_protein_ligands = list(proteins.iterdir())
        sampled_complex = random.sample(curr_protein_ligands,k=1)
        sampled_complexes.append(sampled_complex[0])
    
    return sampled_complexes



if __name__ == "__main__":
    two_stage_sample = two_stage_sampling()
    analyze_sample(two_stage_sample)
    [print(p) for p in two_stage_sample]