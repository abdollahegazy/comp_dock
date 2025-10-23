# eval "$(conda shell.bash hook)"
# conda activate deepsurf

# species="Eucalyptus Arabidopsis DouglasFir"

# for s in $species; do
#     for d in $(ls -d ../Data/$s/*/); do

#         uid=$(basename "$d")
#         mkdir -p ../pockets/$s/$uid
        
#         for i in {0..11}; do
        
#             if [ ! -f "../pockets/$s/$uid/protein_conf${i}/centers.txt" ]; then
#                 echo "Running prediction for $s - $uid - protein_conf${i}"
#                 export PYTHONWARNINGS="ignore::FutureWarning"
#                 python3 /tank/abdolla/DeepSurf/predict.py -p "../Data/$s/$uid/protein_conf${i}.pdb" -mp /tank/abdolla/DeepSurf/models -o "../pockets/$s/$uid"
#             fi
#         done
#     done
# done

# conda deactivate



eval "$(conda shell.bash hook)"
conda activate deepsurf

species="Eucalyptus Arabidopsis DouglasFir"

for s in $species; do
    for d in $(ls -d ../Data/$s/*/); do

        uid=$(basename "$d")
        mkdir -p ../pockets/$s/$uid
        
        for i in {0..11}; do
            target_dir="../pockets/$s/$uid/protein_conf${i}"
            
            # Check if centers.txt is missing AND the directory is empty
            if [ ! -f "$target_dir/centers.txt" ] && [ ! -d "$target_dir" ] || [ -z "$(ls -A "$target_dir" 2>/dev/null)" ]; then
                echo "Running prediction for $s - $uid - protein_conf${i}"
                mkdir -p "$target_dir"
                export PYTHONWARNINGS="ignore::FutureWarning"
                python3 /tank/abdolla/DeepSurf/predict.py -p "../Data/$s/$uid/protein_conf${i}.pdb" -mp /tank/abdolla/DeepSurf/models -o "$target_dir"
            else
                echo "Skipping $s - $uid - protein_conf${i} (already processed or not empty)"
            fi
        done
    done
done

conda deactivate


# nohup bash 07_binding_pockets.sh > binding_pockets.out 2>&1 &
# tail -f binding_pockets.out 