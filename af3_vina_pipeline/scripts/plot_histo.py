import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

atoms= pd.DataFrame(pd.read_csv('/dogwood/tank/Duncan/IMPACTS2022/scripts/aggreagte.csv'))
heads=np.array(atoms.columns.values)
for i in range(len(heads)):
	heads[i]= heads[i].removesuffix('_pct')
	if heads[i]=='DouglasFir':
		heads[i]='Douglas Fir'

hydro=atoms[atoms['Status']=='Hydrophobic Side Chains']
hydro=np.array(hydro.iloc[0:1,1:])
negative=atoms[atoms['Status']=='Negative Charged Side Chains']
negative=np.array(negative.iloc[0:1,1:])
polar=atoms[atoms['Status']=='Polar Uncharged Chains']
polar=np.array(polar.iloc[0:1,1:])
positive=atoms[atoms['Status']=='Positive Charged Side Chains']
positive=np.array(positive.iloc[0:1,1:])
special=atoms[atoms['Status']=='Special Cases']
special=np.array(special.iloc[0:1,1:])
aro=atoms[atoms['Status']=='Aromatic']
aro=np.array(aro.iloc[0:1,1:])
plt.bar(heads[1:5],hydro[0])
plt.bar(heads[1:5],negative[0],bottom=hydro[0],color='r')
plt.bar(heads[1:5],polar[0],bottom=hydro[0]+negative[0],color='y')
plt.bar(heads[1:5],positive[0],bottom=hydro[0]+negative[0]+polar[0],color='g')
plt.bar(heads[1:5],special[0],bottom=hydro[0]+negative[0]+polar[0]+positive[0],color='purple')
plt.bar(heads[1:5],aro[0],bottom=hydro[0]+negative[0]+polar[0]+positive[0]+special[0],color='orange')
plt.legend(['Hydrophobic','Negatively Charged','Polar Uncharged','Positively Charged','Special Cases','Aromatic'],loc="lower right",fontsize='8')
plt.savefig('/dogwood/tank/Duncan/IMPACTS2022/scripts/amino_acids_content.png',dpi=1200)
plt.ylabel("Proportion")
plt.show()
