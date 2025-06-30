#this is newer than v3 btw

# Load required packages
package require psfgen
package require solvate
package require autoionize
package require pbctools

proc solvate_and_ionize {psf_file pdb_file output_prefix} {
    # Solvate the system
    solvate $psf_file $pdb_file -o "${output_prefix}_solvated" -minmax {{-60 -60 -60} {60 60 60}}
    puts "System solvated: ${output_prefix}_solvated"

    # Ionize the system
    autoionize -psf "${output_prefix}_solvated.psf" -pdb "${output_prefix}_solvated.pdb" -sc 0.15 -o "${output_prefix}_system"
    puts "System ionized: ${output_prefix}_system"
}

# Load topology files - this may need to change depending on the system, but I think it should be fine
topology top_all36_prot.rtf
topology toppar/top_all36_cgenff.rtf
topology toppar/2PG_deprot.str
topology /dogwood/home/dboren/toppar_c36_jul21/toppar_water_ions.str
# Alias for residues not recognized by CHARMM
#This is unlikely to be correct. You need to 
pdbalias residue HIS HSD ;# Alias HIS to HSD (neutral form with delta protonation)
pdbalias atom ILE CD1 CD ;# Example alias for atoms, adjust if required


#this loop was originally intended to loop over 7 directories with different proteins - you'll have to change it a bit

foreach dir [list "cnpgp_pga_mg" "rsg2npgp_pga_mg" "rsg2spgp_pga_mg" "rspgp_pga_mg" "sepgp_pga_mg" "tapgp_pga_mg" "tfpgp_pga_mg"] {
    puts "DIR IS $dir"
    mol new $dir/${dir}_model.cif
    #This cif file must be turned into a pdb.
    set selA [atomselect top "protein"] #this can probably just be "protein"
    $selA set segname PA
    $selA writepdb chainA.pdb
    segment PA {
        pdb chainA.pdb
    }
    coordpdb chainA.pdb PA
    guesscoord
    writepsf $dir.psf 
    writepdb $dir.pdb

    solvate_and_ionize $dir.psf $dir.pdb $dir
    pbc writexst $dir.xsc
    resetpsf
}




