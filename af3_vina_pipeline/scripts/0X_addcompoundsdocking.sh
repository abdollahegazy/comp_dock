
species="Arabidopsis DouglasFir Eucalyptus Human"

molNums=()

for f in ../SmallMolecules/molfiles/*.mol;
do
    #echo ${f%.mol}
    #echo('s/[^0-9]*//g')
    filename=$(basename -- "$f")
    numeric=$(echo "$filename" | sed 's/[^0-9]*//g')
    #echo "$numeric"
    molNums+=("$numeric")
done

echo "IDs: ${molNums[@]}"

for s in $species
do
    for uid in $(ls ../Simulations/$s/pockets/)
    do
        for small in "${molNums[@]}"
        do
            if [ ! -f ../Simulations/$s/docking/$uid/$small ] 
            then
                echo "mkdir -p ../Simulations/$s/docking/$uid/$small"
            else
                echo "directory already exists"
            fi
        done
    done
done
