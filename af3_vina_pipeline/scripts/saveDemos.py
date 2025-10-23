import os
import shutil

species = ['Arabidopsis','DouglasFir','Eucalyptus','Human']
subdirectories = ['alphafold']

for s in species:
	for sd in subdirectories:
		srcPath=f"../Simulations/{s}/{sd}/"
		print(srcPath)
		targetPath=f"/tank/Duncan/piplineGit/borenpipeline/demo/testSimulation/{s}/{sd}/"
		
		dirList = os.listdir(srcPath)
		dirEntry = [dir.split('.')[0] for dir in dirList]
		
		for i in range(99):
			print (dirEntry[i])
			print(i,srcPath+dirEntry[i])	
			srcDir = os.path.join(srcPath, dirEntry[i])
			targetDir = os.path.join(targetPath, dirEntry[i])
			shutil.copytree(srcDir, targetDir)
			print(f"Copied {dirEntry[i]} to {targetPath}")
		
		#os.system('cp -r {srcPath}+{dirEntry[i]},{targetPath}')
	
