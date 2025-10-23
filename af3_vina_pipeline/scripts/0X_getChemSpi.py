from chemspipy import ChemSpider
import subprocess
import os
cs = ChemSpider('TzfV80skfR1yCAe7oY4y06Q2s5Kzlff4a1NTNVq2')

#c = cs.get_compound(780) #auxin
listofcompounds = [6309,780,4444606,6085,2169,23208,4444418,6223,3698,395716,
30999,21105998,388386,13876103,917,2424,4450907,7930,10194105,84989,108426,393012,
4518347,514]

for el in listofcompounds:
	if not os.path.exists("../SmallMolecules/pdbfiles/%d.pdb" % el):
		details = cs.get_details(el)
		fout = open("../SmallMolecules/spreadsheets/%d.csv", "w")
		fout.write("%d, %s, %f, %s\n" % (details['id'], details['smiles'], details['molecularWeight'], details['commonName']))
		fout.close()
		mol3D = details['mol3D']

		molFile = open("../SmallMolecules/molfiles/%d.mol" % el, "w") #this should specificy a less obtrusive location. Perhaps small compounds belong in ../Data/SmallCompounds?
		n = molFile.write(mol3D)
		molFile.close()
		subprocess.call("obabel ../SmallMolecules/molfiles/%d.mol -opdb -O../SmallMolecules/pdbfiles/%d.pdb" % (el,el), shell=True)
#Compound.mol_2d
#Compound.mol_3d
#return mol file for given compound


