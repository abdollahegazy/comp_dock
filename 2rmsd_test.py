import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict


AF3_COMPLEXES = Path("./AF_complexes/")
MD_COMPLEXES = Path("./MD_complexes_c3/")
SPECIES = ["Arabidopsis","DouglasFir","Eucalyptus","Human"]

LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930, # 2169 omitted bc its always empty on MD sims
]

def calculate_ligand_rmsd(af3_complex,md_complex):

    u_md = mda.Universe(md_complex)
    u_af3 = mda.Universe(af3_complex)

    prot_md = u_md.select_atoms("protein and name CA")  #Ca backbone because just using 'backbone' resulted in 1 extra atom 
    prot_af3 = u_af3.select_atoms("protein and name CA")
    #there was soem issue about Oxygen being gone in the MD i believe

    #proteins need to be same size
    if prot_md.n_atoms != prot_af3.n_atoms:
        raise ValueError(f"MD ligand has {prot_md.n_atoms} atoms, AF3 ligand has {prot_af3.n_atoms} atoms")
    
    #aligns MD onto AF3
    align.alignto(u_md, u_af3, select="protein and name CA") #align proteins
    #just realized i align the universes, should align the proteins

    lig_md = u_md.select_atoms("segid L and not element H") #basically the MD ligand is missing an extra H always?
    lig_af3 = u_af3.select_atoms("segid B and not element H")

    if lig_md.n_atoms != lig_af3.n_atoms:
        raise ValueError(f"MD ligand has {lig_md.n_atoms} atoms, AF3 ligand has {lig_af3.n_atoms} atoms")

    # compute ligand RMSD (no further superposition needed)
    r = rmsd(lig_md.positions,
            lig_af3.positions,
            center=False,      # already aligned
            superposition=False)
    
    return r

#manually inspect some outlying structures

#mean rmsd for a given protein across all ligands
#with protein's surface area/ solvet accessible surface area (MDAnalysis, sasa)
#logically follows as ligand has more area to chose from
#plot against the proteins sasa
#want to see if bigger proteins give more rmsd across ligands

#remaining rmsd info
#how reppresenetive are these results of a protein in the wild
#take sample of our proteins and run in MD to see if any changes occur (2ndary structure,binding sites?)

def run_rmsd():


    species_rmsd = defaultdict(lambda: [0, 0.0]) #dict in form of species: (#proteins,tot RMSD)
    ligand_rmsd  = defaultdict(lambda: [0, 0.0]) #dict in form of ligand_id: (#proteins,tot RMSD)


    for species in MD_COMPLEXES.iterdir():

        complexes = list(species.iterdir())
        species_name = species.name

        for prtn_path_md in tqdm(complexes,total=len(complexes),desc=f"running RMSD for {species}"):
            prtn_path_af3 = AF3_COMPLEXES / species.name / prtn_path_md.name

            for ligand in LIGANDS:
                cmplx_md = prtn_path_md / str(ligand) / "complex.pdb"

                cmplx_af3 = prtn_path_af3 / str(ligand) / "complex.pdb"
                if not (cmplx_md.exists() and cmplx_af3.exists()):
                    continue

                r = calculate_ligand_rmsd(cmplx_af3,cmplx_md)

                species_stats = species_rmsd[species_name]
                species_stats[0] += 1       # increment count
                species_stats[1] += r       # accumulate RMSD

                lig_stats = ligand_rmsd[ligand]
                lig_stats[0] += 1
                lig_stats[1] += r
                # break
            # break

    return species_rmsd,ligand_rmsd


if __name__ == "__main__":
    species_rmsd,ligand_rmsd = run_rmsd()
    species_rmsd = dict(species_rmsd)
    ligand_rmsd = dict(ligand_rmsd)

    print(species_rmsd)
    print()
    print(ligand_rmsd)


    for name, (n, total) in species_rmsd.items():
        avg = total / n if n else float('nan')
        print(f"{name:12s}  {n:3d} complexes  avg RMSD = {avg:.3f} Å")

    print()

    for lig, (n, total) in ligand_rmsd.items():
        avg = total / n if n else float('nan')
        print(f"{lig:8d}  {n:3d} hits      avg RMSD = {avg:.3f} Å")




    # md_set = {(atom.resid, atom.segid, atom.name) for atom in prot_md}
    # af3_set = {(atom.resid, atom.segid, atom.name) for atom in prot_af3}

    # 2) find the MD atoms that aren’t in AF3
    # extra = af3_set - md_set

    # print(f"Found {len(extra)} extra MD backbone atoms:")
    # for resid, segid, name in sorted(extra):
    #     print(f"  resid={resid:4d}, segid={segid!r}, atom name={name!r}")