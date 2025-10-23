package require psfgen
package require solvate
package require autoionize
package require pbctools

set species {"Arabidopsis" "DouglasFir" "Eucalyptus" "Human"}

foreach species $species {
    set i 0
    foreach d [glob ../Simulations/af_pdb_processed/$species/*.pdb] {
        # puts $d
        set uid [file rootname [file tail $d]]
        
        if { ! [file exists ../Simulations/$species/equilibration/$uid] || 
            ! [file exists ../Simulations/$species/equilibration/$uid/$uid.psf] ||
            ! [file exists ../Simulations/$species/equilibration/$uid/$uid.pdb]} {
            resetpsf
            topology toppar/top_all36_prot.rtf
            pdbalias residue HIS HSE
            pdbalias atom ILE CD1 CD
            segment P {pdb ../Simulations/af_pdb_processed/$species/$uid.pdb}
            coordpdb ../Simulations/af_pdb_processed/$species/$uid.pdb P
            guesscoord
            file mkdir ../Simulations/$species/equilibration/$uid
            writepsf ../Simulations/$species/equilibration/$uid/$uid.psf
            writepdb ../Simulations/$species/equilibration/$uid/$uid.pdb
            
            solvate ../Simulations/$species/equilibration/$uid/$uid.psf ../Simulations/$species/equilibration/$uid/$uid.pdb -o ../Simulations/$species/equilibration/$uid/solvate -t 10
            autoionize -psf ../Simulations/$species/equilibration/$uid/solvate.psf -pdb ../Simulations/$species/equilibration/$uid/solvate.pdb -sc 0.15 -o ../Simulations/$species/equilibration/$uid/ionize
            
            mol new ../Simulations/af_pdb_processed/$species/$uid.pdb
            set prot [atomselect top "protein and name CA"]
            set molid [mol new ../Simulations/$species/equilibration/$uid/ionize.psf]
            mol addfile ../Simulations/$species/equilibration/$uid/ionize.pdb $molid
            animate write pdb ../Simulations/$species/equilibration/$uid/system.pdb
            animate write psf ../Simulations/$species/equilibration/$uid/system.psf
            pbc writexst ../Simulations/$species/equilibration/$uid/system.xsc
            set betas [$prot get beta]
            set resids [$prot get resid]
            for {set i 0} {$i < [llength $resids]} {incr i} {
                set prot2 [atomselect top "protein and resid [lindex $resids $i]"]
                $prot2 set beta [lindex $betas $i]
                $prot2 delete
            }
            $prot delete
            set all [atomselect top "all"]
            set sel [atomselect top "protein and name N C O CA CB"]
            set sel2 [atomselect top "protein and beta > 0.7 and name N C O CA CB"]
            $all set beta 0 
            $sel set beta 1
            animate write pdb ../Simulations/$species/equilibration/$uid/restrain.pdb
            $all set beta 0
            $sel2 set beta 1
            animate write pdb ../Simulations/$species/equilibration/$uid/restrain2.pdb
            $sel2 delete
            $sel delete
            $all delete
            
            foreach m [molinfo list] {mol delete $m}
            
            file copy -force equil/eq.sh ../Simulations/$species/equilibration/$uid/
            file copy -force equil/eq2.sh ../Simulations/$species/equilibration/$uid/
            file copy -force equil/run.sh ../Simulations/$species/equilibration/$uid/
            file copy -force equil/eq.namd ../Simulations/$species/equilibration/$uid/
            file copy -force equil/eq2.namd ../Simulations/$species/equilibration/$uid/
            file copy -force equil/run.namd ../Simulations/$species/equilibration/$uid/
            file copy -force equil/toppar ../Simulations/$species/equilibration/$uid/
            exec sed -i "s/DUMMY_NAME/$uid/g" ../Simulations/$species/equilibration/$uid/eq.sh
            exec sed -i "s/DUMMY_NAME/$uid/g" ../Simulations/$species/equilibration/$uid/eq2.sh
            exec sed -i "s/DUMMY_NAME/$uid/g" ../Simulations/$species/equilibration/$uid/run.sh
        }
    }
}
quit
