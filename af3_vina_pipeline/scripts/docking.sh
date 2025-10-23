module load Miniforge3
eval "$(conda shell.bash hook)"


molList=($(basename -s .pdb ../SmallMolecules/pdbfiles/*.pdb))

species=("Arabidopsis" "DouglasFir" "Eucalyptus" "Human")

conda activate adfr

for s in $species
do
    for uid in $(ls ../pockets/$s)
    do
        mkdir -p ../docking/$s/$uid
        for conf in {0..11}
        do
            file= ../pockets/$s/$uid/protein_conf$conf/centers.txt
            if test -f $file; then
                if [ ! -f ../docking/$s/$uid/protein_conf$conf.pdbqt]; then
                    prepare_receptor -r ../Data/$s/$uid/protein_conf$conf.pdb \
                    -o ../docking/$s/$uid/protein_c$conf.pdbqt
                    conda deactivate
                fi
                # We need to loop through the list of chemspider compounds
                for csid in $molList
                do
                    mkdir -p ../docking/$s/$uid/$csid
                    counter=0 
                    
                    to_do=""
                    while IFS=" *" read -r bx by bz
                    do
                        if [ ! -f ../docking/$s/$uid/$csid/out_c${conf}_p${counter}.log ]; then
                        echo "$bx $by $bz"
cat > ../docking/$s/$uid/$csid/protein_c${conf}_p${counter}.py <<END
from vina import Vina
v = Vina(sf_name='vina')
v.set_receptor('../Simulations/$s/docking/$uid/protein_c$conf.pdbqt')

v.set_ligand_from_file('../Simulations/$s/docking/$uid/$csid/ligand_c${conf}_p${counter}.pdbqt')
v.compute_vina_maps(center=[$bx, $by, $bz], box_size=[20, 20, 20])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('../Simulations/$s/docking/$uid/$csid/out_minimize_c${conf}_p${counter}.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('../Simulations/$s/docking/$uid/$csid/out_poses_c${conf}_p${counter}.pdbqt', n_poses=5, overwrite=True)
END
cat > ../Simulations/$s/docking/$uid/$csid/ligand_c${conf}_p${counter}.tcl <<END
mol new ../SmallMolecules/pdbfiles/${csid}.pdb
set a [atomselect top all]
\$a moveby [vecsub {$bx $by $bz} [measure center \$a]]
\$a writepdb ../Simulations/$s/docking/$uid/$csid/ligand_c${conf}_p${counter}.pdb
quit
END
                            to_do="${to_do} c${conf}_p${counter}"
                        fi
                        counter=$(expr $counter + 1)
                    done < "$file"
        
                    for f in ${to_do}
                    do
                        vmd -dispdev text -e ../Simulations/$s/docking/$uid/$csid/ligand_$f.tcl
                        conda activate vina
                        mk_prepare_ligand.py -i ../Simulations/$s/docking/$uid/$csid/ligand_$f.pdb \
                                          --pH 7.4 \
                                          -o ../Simulations/$s/docking/$uid/$csid/ligand_$f.pdbqt
                        #can I just add an if statement here? like the one below?
                        #if [ ! -f ../Simulations/$s/docking/$uid/$csid/out_$f.log ]; then
                        python3 ../Simulations/$s/docking/$uid/$csid/protein_$f.py > ../Simulations/$s/docking/$uid/$csid/out_$f.log
                        conda deactivate
                    done
                done
            fi
        done
    done
done
