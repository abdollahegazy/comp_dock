import os
import subprocess
from pathlib import Path

from vmd import evaltcl # type: ignore
from vmd import Molecule # type: ignore

# Define paths and species
mol_list = [p.stem for p in Path('../SmallMolecules/pdbfiles').glob('*.pdb')]
species = ["Human"]

def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True,stdout=subprocess.DEVNULL)

def prepare_receptor(s, uid, conf):
    receptor_file = f"../docking/{s}/{uid}/protein_conf{conf}.pdbqt"
    input_file = f"../Data/{s}/{uid}/protein_conf{conf}.pdb"

    if not Path(receptor_file).exists():
        print(f"Preparing receptor for {s}/{uid}, conf {conf}")
        run_command(
            f"bash -c 'module purge && "
            f"module load Miniforge3 && "
            f"conda activate adfr && "
            f"python3 /mnt/ffs24/home/hegazyab/.conda/envs/adfr/lib/python3.8/site-packages/AutoDockTools/Utilities24/prepare_receptor4.py  -r {input_file} -o {receptor_file}'"
        )

def prepare_ligand(s, uid, csid, conf, counter, bx, by, bz):
    target_dir = Path(f"../docking/{s}/{uid}/{csid}")
    target_dir.mkdir(parents=True, exist_ok=True)

    Generate ligand TCL script
    tcl_file = target_dir / f"ligand_c{conf}_p{counter}.tcl"
    with tcl_file.open('w') as f:
        f.write(f"""mol new ../SmallMolecules/pdbfiles/{csid}.pdb
set a [atomselect top all]
$a moveby [vecsub {{{bx} {by} {bz}}} [measure center $a]]
$a writepdb ../docking/{s}/{uid}/{csid}/ligand_c{conf}_p{counter}.pdb
quit
""")

    # Generate docking Python script
    py_file = target_dir / f"protein_c{conf}_p{counter}.py"
    with py_file.open('w') as f:
        f.write(f"""
from vina import Vina
v = Vina(sf_name='vina')
v.set_receptor('../docking/{s}/{uid}/protein_conf{conf}.pdbqt')
v.set_ligand_from_file('../docking/{s}/{uid}/{csid}/ligand_c{conf}_p{counter}.pdbqt')
v.compute_vina_maps(center=[{bx}, {by}, {bz}], box_size=[20, 20, 20])

energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

energy_minimized = v.optimize()
print('Score after minimization: %.3f (kcal/mol)' % energy_minimized[0])

v.write_pose('../docking/{s}/{uid}/{csid}/out_minimize_c{conf}_p{counter}.pdbqt', overwrite=True)
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('../docking/{s}/{uid}/{csid}/out_poses_c{conf}_p{counter}.pdbqt', n_poses=5, overwrite=True)
""")
        

def vmd_prep(s, uid, csid, conf, counter,bx,by,bz):
    mol = Molecule.Molecule()

    mol.load(f'../SmallMolecules/pdbfiles/{csid}.pdb', 'pdb')

    evaltcl(f'set a [atomselect {int(mol)} all]')

    evaltcl(f'$a moveby [vecsub {{{bx} {by} {bz}}} [measure center $a]]')
    
    output_file = f'../docking/{s}/{uid}/{csid}/ligand_c{conf}_p{counter}.pdb'

    evaltcl(f'$a writepdb {output_file}')

    evaltcl('$a delete')

    mol.delete()



def run_docking(s, uid, csid, conf, counter,bx,by,bz):
    target_dir = Path(f"../docking/{s}/{uid}/{csid}")
    log_file = target_dir / f"out_c{conf}_p{counter}.log"
    
    if not log_file.exists():
        vmd_prep(s,uid,csid,conf,counter,bx,by,bz)



        # Prepare ligand
        run_command(
            f"bash -c 'module purge && "
            f"module load Miniforge3 && "
            f"conda activate vina && "
            f"mk_prepare_ligand.py -i {target_dir}/ligand_c{conf}_p{counter}.pdb "
            f"--pH 7.4 -o {target_dir}/ligand_c{conf}_p{counter}.pdbqt' "
        )
        
        # Run docking
        run_command(
            f"bash -c 'module purge && "
            f"module load Miniforge3 && "
            f"conda activate vina && "
            f"python3 {target_dir}/protein_c{conf}_p{counter}.py > {log_file}' "
        )

def main():
    for s in species:
        pockets = Path(f"../pockets/{s}")

        if not pockets.exists():
            continue
        
        for uid in pockets.iterdir():
            uid = uid.name
            target_dir = Path(f"../docking/{s}/{uid}")
            target_dir.mkdir(parents=True, exist_ok=True)
            
            for conf in range(12):
                centers_file = Path(f"../pockets/{s}/{uid}/protein_conf{conf}/centers.txt")
                
                if centers_file.exists():
                    prepare_receptor(s, uid, conf)

                    
                    for csid in mol_list:
                        counter = 0
                        with centers_file.open() as f:
                            for line in f:
                                bx, by, bz = map(float, line.split())
                                log_file = Path(f"../docking/{s}/{uid}/{csid}/out_c{conf}_p{counter}.log")
                                
                                if not log_file.exists():
                                    prepare_ligand(s, uid, csid, conf, counter, bx, by, bz)
                                    run_docking(s, uid, csid, conf, counter,bx,by,bz)
                                counter += 1

if __name__ == "__main__":
    main()
