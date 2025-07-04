import freesasa
freesasa.setVerbosity(freesasa.nowarnings)

import json
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from pathlib import Path
from tqdm import tqdm
import numpy as np

# Paths
DATA_PATH = Path("./temp")                    # your MD “protein.pdb” directory
AF3_PATH   = Path("./AF_complexes")           # AF3 complexes
MD_PATH    = Path("./MD_complexes_c3")        # MD complexes per ligand

# Parameters
SPECIES = ['Arabidopsis','Eucalyptus','Human','DouglasFir']
LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930,
]
POCKET_CUTOFF = 5.0  # Å

# Output structure
sasas = {s: {} for s in SPECIES}

for species_dir in DATA_PATH.iterdir():
    species = species_dir.name
    for protein_dir in tqdm(list(species_dir.iterdir()),
                            desc=f"{species} proteins"):
        pid = protein_dir.name
        pdb_md = protein_dir / "protein.pdb"
        if not pdb_md.exists():
            continue

        # 1) Global MD SASA (single protein.pdb)
        try:
            sasa_md = SASAAnalysis(
                mda.Universe(pdb_md).select_atoms("protein")
            ).run().results.total_area[0]
        except Exception as e:
            print(f"[MD global failed] {species}/{pid}: {e}")
            continue

        # 2) Per‐ligand AF3 SASA & pocket AF3 SASA
        af3_global_sasas = []
        af3_pocket_sasas = []
        af3_base = AF3_PATH / species / pid

        # 3) Per‐ligand MD pocket SASA
        md_pocket_sasas = []

        for lig in LIGANDS:
            af3_pdb = af3_base / str(lig) / "complex.pdb"
            md_pdb  = MD_PATH / species / pid / str(lig) / "complex.pdb"

            # — AF3 global
            if af3_pdb.exists():
                u_af3 = mda.Universe(af3_pdb)
                prot_af3 = u_af3.select_atoms("protein")
                try:
                    af3_global_sasas.append(
                        SASAAnalysis(prot_af3).run().results.total_area[0]
                    )
                except Exception as e:
                    print(f"[AF3 global failed] {species}/{pid}, lig {lig}: {e}")

                # — AF3 pocket
                lig_sel = u_af3.select_atoms("not protein and not resname HOH")
                if len(lig_sel) > 0:
                    dmat = distance_array(prot_af3.positions,
                                          lig_sel.positions)
                    close_idx = np.unique(np.where(dmat < POCKET_CUTOFF)[0])
                    if close_idx.size:
                        try:
                            af3_pocket_sasas.append(
                                SASAAnalysis(prot_af3[close_idx]).run().results.total_area[0]
                            )
                        except Exception as e:
                            print(f"[AF3 pocket failed] {species}/{pid}, lig {lig}: {e}")

            # — MD pocket
            if md_pdb.exists():
                u_md = mda.Universe(md_pdb)
                prot_md_complex = u_md.select_atoms("protein")
                lig_md = u_md.select_atoms("not protein and not resname HOH")
                if len(lig_md) > 0:
                    dmat_md = distance_array(prot_md_complex.positions,
                                             lig_md.positions)
                    close_idx_md = np.unique(np.where(dmat_md < POCKET_CUTOFF)[0])
                    if close_idx_md.size:
                        try:
                            md_pocket_sasas.append(
                                SASAAnalysis(prot_md_complex[close_idx_md])
                                    .run().results.total_area[0]
                            )
                        except Exception as e:
                            print(f"[MD pocket failed] {species}/{pid}, lig {lig}: {e}")

        # Skip proteins with no AF3 data
        if not af3_global_sasas:
            continue

        # Compute means
        sasa_af3_global   = float(np.mean(af3_global_sasas))
        sasa_af3_pocket   = float(np.mean(af3_pocket_sasas)) if af3_pocket_sasas else None
        sasa_md_pocket     = float(np.mean(md_pocket_sasas))   if md_pocket_sasas else None

        # Derive metrics
        delta_global      = sasa_md - sasa_af3_global
        abs_delta_global  = abs(delta_global)
        ratio_global      = sasa_md / sasa_af3_global
        avg_global        = 0.5 * (sasa_md + sasa_af3_global)

        pocket_delta      = (sasa_md_pocket - sasa_af3_pocket)     if (sasa_md_pocket and sasa_af3_pocket) else None
        pocket_ratio      = (sasa_md_pocket / sasa_af3_pocket)     if (sasa_md_pocket and sasa_af3_pocket) else None
        pocket_avg        = (0.5 * (sasa_md_pocket + sasa_af3_pocket)) if (sasa_md_pocket and sasa_af3_pocket) else None

        # Store
        sasas[species][pid] = {
            "md_global":      float(sasa_md),
            "af3_global":     sasa_af3_global,
            "delta_global":   float(delta_global),
            "abs_delta":      float(abs_delta_global),
            "ratio_global":   float(ratio_global),
            "avg_global":     float(avg_global),
            "md_pocket":      sasa_md_pocket,
            "af3_pocket":     sasa_af3_pocket,
            "delta_pocket":   pocket_delta,
            "ratio_pocket":   pocket_ratio,
            "avg_pocket":     pocket_avg
        }

# Write out
with open("sasa_conf3_extended_metrics.json", "w") as f:
    json.dump(sasas, f, indent=4)
