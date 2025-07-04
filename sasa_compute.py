import freesasa
freesasa.setVerbosity(freesasa.nowarnings)

import hashlib
import json
import MDAnalysis as mda
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from pathlib import Path
from tqdm import tqdm
import numpy as np

DATA_PATH = Path("./temp")
AF3_PATH = Path("./AF_complexes")

SPECIES = ['Arabidopsis','Eucalyptus','Human','DouglasFir']
LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930,
]

sasas = {s: {} for s in SPECIES}

for species in DATA_PATH.iterdir():
    species_name = species.name
    proteins = list(species.iterdir())
    for protein in tqdm(proteins, desc=f'iterating proteins for {species_name}'):

        protein_id = protein.name
        
        pdb_md = protein / 'protein.pdb'

        if not pdb_md.exists(): continue  # skip if MD protein missing

        u_md = mda.Universe(pdb_md)
        resnames_md = u_md.residues.resnames
        u_md.residues.resnames = np.where(resnames_md=="HIS","HSE",resnames_md)

        prot_md = u_md.select_atoms("protein")

        sasa_md = SASAAnalysis(mda.Universe(pdb_md).select_atoms("protein")).run().results.total_area[0]



        af3_sasas = {}
        sasa_cache = {}
        af3_protein_base = AF3_PATH / species_name / protein_id

        c = 1

        for lig in LIGANDS:
            pdb_af3 = af3_protein_base / str(lig) / "complex.pdb"
            if not pdb_af3.exists(): continue

            u_af3 = mda.Universe(pdb_af3)
            resnames_af3 = u_af3.residues.resnames
            u_af3.residues.resnames = np.where(resnames_af3=="HIS","HSE",resnames_af3)

            prot_af3 = u_af3.select_atoms("protein")

            coords = prot_af3.positions.round(3).tobytes()
            h = hashlib.md5(coords).hexdigest()

            if h not in sasa_cache:
                c+=1
                sasa_cache[h] = SASAAnalysis(u_af3).run().results.total_area[0]
            af3_sasas[lig] = sasa_cache[h]

        if not af3_sasas:
            continue  # skip if no AF3 complexes found


        sasas[species_name][protein_id] = {
            "md_sasa":float(sasa_md),
            'af3_sasa':af3_sasas
        }


with open("sasa_md_and_af3_per_ligand.json", 'w') as f:
    json.dump(sasas, f, indent=4)


