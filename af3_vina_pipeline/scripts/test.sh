#!/bin/bash

molList=""

for d in $(ls ../SmallMolecules/pdbfiles/*)
    do
        id=${d##../SmallMolecules/pdbfiles/}
        echo $id
        uid=${id%.pdb}
        echo $uid
        molList="${molList} $uid"
done
echo $molList