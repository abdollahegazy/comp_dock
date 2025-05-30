from pymol import cmd

# 1. Load the two structures
cmd.load("MD_complexes_c3/Arabidopsis/a0a1p8aql0/780/complex.pdb", "MD")
cmd.load("AF_complexes/Arabidopsis/a0a1p8aql0/780/complex.pdb", "AF3")

# 2. Superpose MD → AF3 on Cα atoms
cmd.align("MD and name CA", "AF3 and name CA")

# 3. Style & color
cmd.hide("everything", "all")
cmd.show("cartoon", "AF3")
cmd.show("cartoon", "MD")
cmd.color("blue", "AF3")
cmd.color("red",  "MD")

# 4. Zoom and render
cmd.zoom()  
cmd.png("aligned_complexes.png", width=1200, height=800, ray=1)
