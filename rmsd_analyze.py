import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from pathlib import Path
from tqdm import tqdm
import pandas as pd

# Provided data, this is from conformation 9, just using for now.
ligand_data = {
    10194105: [1750, 34104.09647302055],
    21105998: [1751, 34047.16816866223], 
    2424: [1750, 34404.03048557364], 
    388386: [1750, 34081.10839882885], 
    4444418: [1750, 37028.99382134324], 
    4518347: [1748, 35725.6029172345],
    6309: [1751, 35554.56173628888], 
    84989: [1751, 33143.3151413672],
    108426: [1751, 33626.71413923767], 
    30999: [1751, 33561.45698269773], 
    393012: [1748, 35406.59849790922], 
    4444606: [1751, 35549.276919959106], 
    6085: [1751, 36063.06989944418], 
    780: [1751,35387.64860777813], 
    917: [1750, 34169.34985776826], 
    13876103: [1751, 33722.7221962813], 
    23208: [1751, 40213.39501269877], 
    3698: [1751, 33799.44586526322],
    395716: [1751, 34082.67336967359], 
    4450907: [1751, 34693.98421066502], 
    6223: [1749, 37050.18554313896], 
    7930: [1745, 34958.83133974555]} 

species_data = {'Arabidopsis': [10313, 235257.842973596], 
                'Human': [9236, 173483.84425913342], 
                'Eucalyptus': [10470, 212914.3264588643], 
                'DouglasFir': [8484, 148718.21589298834]}   

# Create DataFrame for ligands
df_lig = pd.DataFrame.from_dict(
    ligand_data, orient='index', columns=['count', 'total_rmsd']
)
df_lig.index.name = 'ligand_id'
df_lig['avg_rmsd'] = df_lig['total_rmsd'] / df_lig['count']
df_lig = df_lig.sort_values('avg_rmsd')

# Create DataFrame for species
df_sp = pd.DataFrame.from_dict(
    species_data, orient='index', columns=['count', 'total_rmsd']
)
df_sp.index.name = 'species'
df_sp['avg_rmsd'] = df_sp['total_rmsd'] / df_sp['count']
df_sp = df_sp.sort_values('avg_rmsd')
print(df_sp)
print()
print(df_lig)

total_count_lig = sum(c for c, _ in ligand_data.values())
total_rmsd_lig  = sum(r for _, r in ligand_data.values())
global_avg_lig  = total_rmsd_lig / total_count_lig

total_count_sp = sum(c for c, _ in species_data.values())
total_rmsd_sp  = sum(r for _, r in species_data.values())
global_avg_sp  = total_rmsd_sp / total_count_sp

# Prepare output DataFrame
df_global = pd.DataFrame([{
    'category': 'ligands',
    'total_count': total_count_lig,
    'total_rmsd': total_rmsd_lig,
    'global_avg_rmsd': global_avg_lig
}, {
    'category': 'species',
    'total_count': total_count_sp,
    'total_rmsd': total_rmsd_sp,
    'global_avg_rmsd': global_avg_sp
}])
print(df_global)
# # Display to user
# import ace_tools as tools; tools.display_dataframe_to_user(name="Ligand RMSD Summary", dataframe=df_lig)
# tools.display_dataframe_to_user(name="Species RMSD Summary", dataframe=df_sp)
