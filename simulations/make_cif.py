import subprocess
from pathlib import Path

COMPLEX_DIR = Path("complexes")

def convert():

    for species in COMPLEX_DIR.iterdir():
        for protein in species.iterdir():
            for ligand in protein.iterdir():
                print(ligand)
                af_complex = ligand / "AF3_complex.pdb"
                md_complex = ligand / "MD_complex.pdb"

                # af_cif = ligand / "AF3_complex.mmcif"
                # md_cif = ligand / "MD_complex.mmcif"
                
                subprocess.run([
                    "pdb2cif", str(af_complex)] #,"-O",str(af_cif)]
                )
                subprocess.run([
                    "pdb2cif", str(md_complex)] #,"-O",str(md_cif)]
                )
                # quit()



if __name__ == "__main__":
    convert()