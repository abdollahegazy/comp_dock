
"""
Script to merge Duncans protein confs and their ligands 
oriented to binding pocekts into just one docked
structure
"""
from pathlib import Path
from tqdm import tqdm
from pprint import pprint
from subprocess import run,DEVNULL
from collections import defaultdict

from Bio.PDB import PDBParser, PDBIO, Structure


LIGANDS = [
    10194105, 21105998, 2424, 388386, 4444418, 4518347, 6309, 84989,
    108426, 30999, 393012, 4444606, 6085, 780, 917,
    13876103, 23208, 3698, 395716, 4450907, 6223, 7930, # 2169 omitted bc its always empty
]

SPECIES = ["Human", "DouglasFir", "Arabidopsis", "Eucalyptus"]
OUT_DIR = Path("./complexes")
TEMP_DIR = Path("./temp")

def cache_obabel(conf=3,pose=0):
    
    for species in SPECIES:

        PROTEIN_PATH = Path(f"/dogwood/tank/Duncan/IMPACTS2022/Simulations/{species}/docking")

        Path.mkdir(TEMP_DIR/species,exist_ok=True)
        proteins = list(PROTEIN_PATH.iterdir())

        for protein_dir in tqdm(proteins,total=len(proteins),desc=f"processing ligands for {species}"):
            protein_path = protein_dir / f"protein_c{conf}.pdbqt"

            if not protein_path.exists():
                continue
            protein_id = protein_path.parts[-2].lower()

            Path.mkdir(TEMP_DIR/species/protein_id,exist_ok=True)

            protein_out_path = TEMP_DIR / species/ protein_id / "protein.pdb"

            run(["obabel", "-ipdbqt", str(protein_path), "-O", str(protein_out_path)], 
                check=True,
                stderr=DEVNULL,
                stdout=DEVNULL
                )

            for ligand in LIGANDS:
                ligand_path = protein_dir / f"{ligand}" / f"out_poses_c{conf}_p{pose}.pdbqt"

                if not ligand_path.exists():
                    continue
            
                Path.mkdir(TEMP_DIR/species/protein_id/str(ligand),exist_ok=True)
                ligand_out_path = TEMP_DIR / species / protein_id / str(ligand) / f"ligand_c{conf}_p{pose}.pdb"
                cleaned_ligand_out_path = TEMP_DIR / species / protein_id / str(ligand) / f"ligand_c{conf}_p{pose}_cleaned.pdb"
                run([ "obabel", "-ipdbqt", str(ligand_path), "-f", "1", "-l", "1", "-O", str(ligand_out_path) ], 
                    check=True,
                    stderr=DEVNULL,
                    stdout=DEVNULL
                    )
                clean_ligand_pdb(ligand_out_path,cleaned_ligand_out_path)


def clean_ligand_pdb(input_path, output_path):
    counters = defaultdict(int)
    with open(input_path, 'r') as inp, open(output_path, 'w') as out:
        current_residue = None
        for line in inp:
            if line.startswith(('ATOM  ', 'HETATM')):
                # Key residues by (chain, resseq, resname)
                chain = line[21]
                resseq = line[22:26]
                resname = line[17:20]
                residue_key = (chain, resseq, resname)
                # Reset counters when residue changes
                if residue_key != current_residue:
                    counters.clear()
                    current_residue = residue_key
                elem = line[76:78].strip() or line[12:16].strip()[0]
                counters[elem] += 1
                new_name = f"{elem}{counters[elem]}"
                # Write atom name in cols 12-16, right-justified
                newline = (
                    line[:12] +
                    new_name.rjust(4) +
                    line[16:]
                )
                # Optionally, change record name to HETATM
                newline = 'HETATM' + newline[6:]
                out.write(newline)
            else:
                out.write(line)


def merge(protein_path,ligand_path,out_path,parser):
    protein_structure = parser.get_structure("PROT", protein_path)
    ligand_structure = parser.get_structure("LIG", ligand_path)


    for model in protein_structure:
        for chain in model:
            chain.id = 'A'
    
    for model in ligand_structure:
        for chain in model:
            chain.id = 'L'
    
    combined = Structure.Structure("COMPLEX")
    model = protein_structure[0].copy()
    for chain in ligand_structure[0]:
        model.add(chain.copy())
    combined.add(model)

    io = PDBIO()
    io.set_structure(combined)
    io.save(str(out_path))
    
def merge_pdb_files(parser):
    for species in SPECIES:
        temp_path = TEMP_DIR / species
        Path.mkdir(OUT_DIR/species,exist_ok=True)
        n = len(list(temp_path.iterdir()))
        for protein_dir in tqdm(temp_path.iterdir(),total=n,desc=f'merging for {species}'):

            protein_path = protein_dir / "protein.pdb"
            protein_id = protein_path.parts[-2].lower()
            
            Path.mkdir(OUT_DIR/species/protein_id,exist_ok=True)

            for ligand in LIGANDS:
                ligand_path = protein_dir / str(ligand) / "ligand_c3_p0_cleaned.pdb"
                if not ligand_path.exists():
                    continue
                Path.mkdir(OUT_DIR/species/protein_id/str(ligand),exist_ok=True)
                complex_out_path = OUT_DIR/species/protein_id/str(ligand) / "complex.pdb"
                merge(protein_path,ligand_path,complex_out_path,parser)

if __name__ == "__main__":

    Path.mkdir(OUT_DIR,exist_ok=True)
    Path.mkdir(TEMP_DIR,exist_ok=True)

    # cache_obabel()
    
    parser = PDBParser(QUIET=True)
    merge_pdb_files(parser)


    

