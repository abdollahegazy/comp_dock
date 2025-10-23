
species="Arabidopsis DouglasFir Eucalyptus Human"

for s in $species
do
    for uid in $(ls ../Simulations/$s/pockets/)
    do
        mkdir -p ../Simulations/$s/docking/$uid
        mkdir -p ../Simulations/$s/docking/$uid/6309
        mv ../Simulations/$s/docking/$uid/ligand* ../Simulations/$s/docking/$uid/6309/
        mv ../Simulations/$s/docking/$uid/out* ../Simulations/$s/docking/$uid/6309/
        mv ../Simulations/$s/docking/$uid/*py ../Simulations/$s/docking/$uid/6309/
    done
done