# Conformations Directory

Contains 12 PDB conformation files extracted from MD simulations, along with associated PSF files.

## Contents

- **12 PDB conformations:** Extracted from NAMD trajectory data
- **PSF files:** Protein structure files with topology information
- **Additional tweaks:** Post-processing modifications

## Generation

Conformations are created by:
```
06_VMD_analyse_simulations.py
```

This script extracts frames from the MD trajectory and converts them to individual PDB files.