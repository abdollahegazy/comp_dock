# NAMD Simulation

Takes solvated and ionized proteins (from `05_VMD_prepare_simulations.tcl`) and runs NAMD simulations.

## Process

1. **Equilibration 1:** `eq1.namd`
2. **Equilibration 2:** `eq2.namd`  
3. **Production run:** `run.namd` (12ns)
4. **Output:** 12 `.coor` files → 12 PDB conformations

## Note

`a0a7g2e326` in Arabidopsis required 12+ hours for all equilibration simulations.