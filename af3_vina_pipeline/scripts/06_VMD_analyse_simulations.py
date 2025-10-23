import os
import glob
import re

from vmd import evaltcl 
from vmd import Molecule 

species = ["DouglasFir","Eucalyptus","Arabidopsis","Human"]

for s in species:
    base_dir = f'../{s}'
    for d in glob.glob(f'{base_dir}/*'):

        uid = os.path.basename(d)

        if (os.path.exists((os.path.join(base_dir, uid)))
            and 
            not os.path.exists(f'../Data/{s}/{uid}/protein_conf11.pdb')):

            mol = Molecule.Molecule()  
            mol.load(f'{base_dir}/{uid}/system.psf', 'psf')
            mol.load(f'{base_dir}/{uid}/system.pdb', 'pdb')


            for i in range(1, 12):
                coor_file = f"{base_dir}/{uid}/run{i:03d}.coor"
                mol.load(coor_file, 'namdbin')

            evaltcl('package require topotools')
            evaltcl(f'set prot [atomselect {int(mol)} protein]')
            evaltcl(f'topo -molid {int(mol)} -sel $prot guessatom element mass')

            path = f'../Data/{s}/{uid}'
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)

            evaltcl(f'$prot writepsf {path}/protein.psf')

            for f in range(0, mol.numFrames()):
                evaltcl(f'$prot frame {f}')
                evaltcl(f'$prot writepdb {path}/protein_conf{f}.pdb')
            evaltcl('$prot delete')
            mol.delete()

quit()





# # check for existing runXXX.coor files (ignore eq.coor and eq2.coor)
# existing_files = [f for f in os.listdir(d) 
#                 if f.startswith('run') 
#                 and f.endswith('.coor') 
#                 and re.match(r"run\d{3}\.coor", f)]

# expected_files = [f"run{i:03d}.coor" for i in range(1, 12)]
# missing_files = [f for f in expected_files if f not in existing_files]

# # check for extra run files (like run012.coor)
# extra_files = [f for f in existing_files if int(f[3:6]) > 11]

# if missing_files or extra_files:
# print(f"Skipping {uid} in {s} → Missing files: {missing_files} | Extra files: {extra_files}")
# continue  
