from pathlib import Path
from utils import species,mol_list
from vmd import evaltcl # type: ignore
from vmd import Molecule # type: ignore



def vmd_ligand_prep(s,uid,ligand_id,conf,counter,bx,by,bz):
    mol = Molecule.Molecule()

    mol.load(f"../SmallMolecules/pdbfiles/{ligand_id}.pdb",'pdb')

    evaltcl(f'set a [atomselect {int(mol)} all]')

    evaltcl(f'$a moveby [vecsub {{{bx} {by} {bz}}} [measure center $a]]')
    
    output_file = f'../docking/{s}/{uid}/{ligand_id}/ligand_conf{conf}_pose{counter}.pdb'

    evaltcl(f'$a writepdb {output_file}')

    evaltcl('$a delete')

    mol.delete()


def generate_docking_script(s, uid, ligand_id, conf, counter, bx, by, bz):

    vina_script_template =f"""
from vina import Vina
v = Vina(sf_name='vina')
v.set_receptor('/tank/abdolla/docking/{s}/{uid}/protein_conf{conf}.pdbqt')
v.set_ligand_from_file('/tank/abdolla/docking/{s}/{uid}/{ligand_id}/ligand_conf{conf}_pose{counter}.pdbqt')
v.compute_vina_maps(center=[{bx}, {by}, {bz}], box_size=[20, 20, 20])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])

v.write_pose('/tank/abdolla/docking/{s}/{uid}/{ligand_id}/out_minimize_c{conf}_p{counter}.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('/tank/abdolla/docking/{s}/{uid}/{ligand_id}/out_poses_c{conf}_p{counter}.pdbqt', n_poses=5, overwrite=True)
"""
    
    return vina_script_template



for s in species:
    pockets = Path(f"../pockets/{s}")

    for i,uid in enumerate(pockets.iterdir()):
        uid = uid.name
        i+=1
        print(f"Processing {s}/{uid}. Complete {i}/500")
        for conf in range(12):
            center_file = Path(f"../pockets/{s}/{uid}/protein_conf{conf}/centers.txt")
            
            if not center_file.exists():
                continue

            for ligand_id in mol_list:

                ligand_id_dir = Path(f"../docking/{s}/{uid}/{ligand_id}")
                # ligand_id_dir.mkdir(parents=True,exist_ok=True)           

                
                with center_file.open() as f:
                    for pose_counter,line in enumerate(f):
                        bx,by,bz = map(float,line.split())

                        docking_script = generate_docking_script(s, uid, ligand_id, conf, pose_counter, bx, by, bz)

                        py_file = ligand_id_dir / f"protein_conf{conf}_pose{pose_counter}.py"
                        py_file.write_text(docking_script)


                        # vmd_ligand_prep(s,uid,ligand_id,conf,pose_counter,bx,by,bz)

