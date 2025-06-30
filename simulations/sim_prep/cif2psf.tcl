# Load required packages
package require psfgen
package require solvate
package require autoionize
package require pbctools

proc solvate_and_ionize {psf_file pdb_file output_prefix} {
    solvate $psf_file $pdb_file -o "${output_prefix}_solvated" -minmax {{-60 -60 -60} {60 60 60}}
    puts "System solvated: ${output_prefix}_solvated"

    autoionize -psf "${output_prefix}_solvated.psf" -pdb "${output_prefix}_solvated.pdb" -sc 0.15 -o "${output_prefix}_system"
    puts "System ionized: ${output_prefix}_system"
}

# Load topology files
# topology top_all36_prot.rtf
# topology toppar/top_all36_cgenff.rtf
# topology toppar/2PG_deprot.str
# topology /dogwood/home/dboren/toppar_c36_jul21/toppar_water_ions.str

# Residue and atom aliasing
# pdbalias residue HIS HSD
# pdbalias atom ILE CD1 CD

# Loop through all complexes recursively
set complex_dirs [glob -nocomplain -type d complexes/*/*/*]

foreach dir $complex_dirs {
    puts "\nProcessing directory: $dir"

    # foreach tag {AF3 MD} {
    #     set pdb_path "${dir}/${tag}_complex.pdb"

    #     if {[file exists $pdb_path]} {
    #         set output_prefix "${dir}/${tag}"
    #         set psf_file "${output_prefix}.psf"
    #         set pdb_file "${output_prefix}.pdb"

    #         # Load and write initial PSF/PDB
    #         resetpsf
    #         segment PA {
    #             pdb $pdb_path
    #         }
    #         coordpdb $pdb_path PA
    #         guesscoord
    #         writepsf $psf_file
    #         writepdb $pdb_file

    #         # Solvate and ionize
    #         solvate_and_ionize $psf_file $pdb_file $output_prefix
    #         pbc writexst "${output_prefix}.xsc"
    #     } else {
    #         puts "WARNING: Missing $pdb_path"
    #     }
    # }
}
