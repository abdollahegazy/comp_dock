display projection Orthographic
display depthcue off
axes location Off
display shadows on
display ambientocclusion on
color Display Background white
color Name C silver

set path ../Simulations/Arabidopsis/pockets/F4IRF5/protein_conf1

#protein selection

set mol [mol ../Simulations/Arabidopsis/equilibration/F4IRF5/restrain2.pdb]
mol modselect 0 $mol protein
mol modcolor 0 $mol Type
mol modstyle 0 $mol NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 0 $mol AOChalky
mol smoothrep $mol 0 5

for { set conf 1 } { $rep <= 3 } { incr rep } {
	set mol [mol $path/pocket$conf.pdb]

	
	#pocket selection
	
	mol addrep $mol
	mol modselect 0 $mol all
	mol modcolor 0 $mol Name
	mol modstyle 0 $mol VDW 1.000000 12.000000
	mol modmaterial 0 $mol AOChalky
	mol smoothrep $mol 0 5

	
}





