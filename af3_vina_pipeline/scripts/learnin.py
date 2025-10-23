#this file probably shouldn't exist, I'm just learning how to write directory entries to a list
import os

molList = os.listdir('../SmallMolecules/molfiles')
molNums = [x.split('.')[0] for x in molList]

print(molNums)


#for d in $(ls ../SmallMolecules/pdbfiles/*)
#    do
#        id=${d##../SmallMolecules/pdbfiles/}
#        echo $id
#        uid=${id%.pdb}
#        echo $uid
#        molList="${molList} $uid"
#done