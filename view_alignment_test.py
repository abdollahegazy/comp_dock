from pymol import cmd
from MDAnalysis import Merge, Universe
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.PDB import PDBWriter
import os

# # # 1. Load the two structures
# # cmd.load("MD_complexes_c3/Arabidopsis/a0a1p8aql0/780/complex.pdb", "MD")
# # cmd.load("AF_complexes/Arabidopsis/a0a1p8aql0/780/complex.pdb", "AF3")

# # # 2. Superpose MD → AF3 on Cα atoms
# # cmd.align("MD and name CA", "AF3 and name CA")

# # # 3. Style & color
# # cmd.hide("everything", "all")
# # cmd.show("cartoon", "AF3")
# # cmd.show("cartoon", "MD")
# # cmd.color("blue", "AF3")
# # cmd.color("red",  "MD")

# # # 4. Zoom and render
# # cmd.zoom()  
# # cmd.png("aligned_complexes.png", width=1200, height=800, ray=1)



def write_viz_pdb(af3_pdb, md_pdb, out_pdb="viz_merged.pdb"):

    print("AF3 path:", af3_pdb, "exists?", os.path.isfile(af3_pdb))
    print("MD  path:", md_pdb,  "exists?", os.path.isfile(md_pdb))

    u_af3 = Universe(af3_pdb)
    u_md  = Universe(md_pdb)
    print(" u_af3 total atoms:", len(u_af3.atoms))
    print(" u_md  total atoms:", len(u_md.atoms))

    # 1) Select only heavy‐atom protein + ligand
    prot_af3 = u_af3.select_atoms("protein and not element H")
    lig_af3  = u_af3.select_atoms("segid B and not element H")
    prot_md  = u_md.select_atoms("protein and not element H")
    lig_md   = u_md.select_atoms("segid L and not element H")

    assert len(prot_af3), "No AF3 protein found!"
    assert len(prot_md),  "No MD protein found!"
    assert len(lig_af3) == len(lig_md), "Ligand atom counts differ!"

    # 2) Align MD to AF3 on CA atoms
    align.alignto(u_md, u_af3, select="protein and name CA")

    # 3) Merge (don't assign segids before)
    merged = Merge(prot_af3, prot_md, lig_af3, lig_md)

    # 4) Assign segids **AFTER** merging, by slicing
    n_af3   = prot_af3.n_atoms
    n_md    = prot_md.n_atoms
    n_ligaf = lig_af3.n_atoms
    n_ligmd = lig_md.n_atoms

    merged.atoms[:n_af3].segments.segids = ["AFP"]
    merged.atoms[n_af3:n_af3+n_md].segments.segids = ["MDP"]
    merged.atoms[n_af3+n_md:n_af3+n_md+n_ligaf].segments.segids = ["AFL"]
    merged.atoms[n_af3+n_md+n_ligaf:].segments.segids = ["MDL"]

    with PDBWriter(out_pdb) as w:
        w.write(merged.atoms)

    print(f"\nWrote {out_pdb}:")
    print(f"  AF3 protein  = {prot_af3.n_atoms} atoms")
    print(f"  MD  protein  = {prot_md.n_atoms} atoms")
    print(f"  AF3 ligand   = {lig_af3.n_atoms} atoms")
    print(f"  MD  ligand   = {lig_md.n_atoms} atoms")
    print(f"  TOTAL        = {len(merged.atoms)} atoms")



write_viz_pdb(
    "AF_complexes/Arabidopsis/a0a7g2e326/6309/complex.pdb",
    "MD_complexes_c3/Arabidopsis/a0a7g2e326/6309/complex.pdb",
    out_pdb="merged_for_vmd.pdb"
)
