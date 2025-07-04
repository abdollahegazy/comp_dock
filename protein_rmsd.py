import json
import sys
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
from datetime import datetime

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

CONFORMATION = sys.argv[1]  # e.g., 3

AF3_COMPLEXES = Path("./AF_complexes/")
MD_COMPLEXES = Path(f"./MD_complexes_c{CONFORMATION}/")

SPECIES = ["Arabidopsis", "DouglasFir", "Eucalyptus", "Human"]

def calculate_protein_rmsd(af3_complex, md_complex):
    u_md = mda.Universe(md_complex)
    u_af3 = mda.Universe(af3_complex)
    prot_md = u_md.select_atoms("protein and name CA")
    prot_af3 = u_af3.select_atoms("protein and name CA")

    if prot_md.n_atoms != prot_af3.n_atoms:
        raise ValueError(
            f"Protein atom count mismatch: {md_complex} ({prot_md.n_atoms}), {af3_complex} ({prot_af3.n_atoms})"
        )

    align.alignto(prot_md, prot_af3)
    r = rmsd(prot_md.positions, prot_af3.positions, center=False, superposition=False)
    return r

def run_protein_rmsd():
    detailed_prot_rmsd = defaultdict(dict)  # {species: {protein_name: rmsd}}

    for species in MD_COMPLEXES.iterdir():
        species_name = species.name
        complexes = list(species.iterdir())
        for prtn_path_md in tqdm(complexes, total=len(complexes), desc=f"protein RMSD for {species_name}"):
            prtn_path_af3 = AF3_COMPLEXES / species_name / prtn_path_md.name

            # Pick any ligand dir just to load the complex.pdb for the protein structure (choose the first found)
            ligand_dirs = sorted((prtn_path_md).iterdir())
            found = False
            for ligand_dir in ligand_dirs:
                md_pdb = ligand_dir / "complex.pdb"
                af3_pdb = prtn_path_af3 / ligand_dir.name / "complex.pdb"
                if md_pdb.exists() and af3_pdb.exists():
                    try:
                        r = calculate_protein_rmsd(af3_pdb, md_pdb)
                        detailed_prot_rmsd[species_name][prtn_path_md.name] = r
                        found = True
                        break
                    except Exception as e:
                        print(f"Error with {species_name}/{prtn_path_md.name}: {e}")
                        continue
            if not found:
                print(f"No valid complexes for {species_name}/{prtn_path_md.name}")

    return detailed_prot_rmsd

def to_dict(obj):
    if isinstance(obj, defaultdict):
        obj = dict(obj)
    if isinstance(obj, dict):
        return {k: to_dict(v) for k, v in obj.items()}
    else:
        return obj

if __name__ == "__main__":
    detailed_prot_rmsd = run_protein_rmsd()
    output = {
        "metadata": {
            "generated_on": datetime.now().isoformat(),
            "description": "Per-protein backbone RMSD (AF3 vs MD) by species",
            "conformation": CONFORMATION,
        },
        "protein_rmsd": to_dict(detailed_prot_rmsd)
    }
    outname = f"protein_rmsd_conf{CONFORMATION}.json"
    with open(outname, "w") as f:
        json.dump(output, f, indent=4)
    print(f"Wrote {outname}")
