import pandas as pd

df = pd.read_csv('bioactivity_preprocessed_data.csv')

# LIPINSKI descriptors
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
        
    baseData = np.arange(1,1)
    i=0
    for mol in moldata:
        
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])
        
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1
        
    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)
    
    return descriptors

df_lipinski = lipinski(df.canonical_smiles)

df_combined = pd.concat([df, df_lipinski], axis=1)

print(df_combined)

def pIC50(input):
    pIC50 = []
    
    for i in input['standard_value']:
        molar = i*(10**-9)
        pIC50.append(-np.log10(molar))
        
    input['pIC50'] = pIC50
    x = input.drop('standard_value', axis = 1)
    
    return x

df_final = pIC50(df_combined)

print(df_final.columns)
print(df_final.pIC50.describe())

import matplotlib.pyplot as plt

ic50_values = df_final['pIC50']
# plt.hist(ic50_values)
# plt.show()

# Exploratory Data Analysis (Chemical Space Analysis) via Lipinski descriptors

import seaborn as sns
sns.set(style = 'ticks')

plt.figure(figsize=(5.5, 5.5))
sns.countplot(x = 'bioactivity_class', data=df_final, edgecolor='black')
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')
plt.savefig('plot_bioactivity_class.pdf')

# Scatterplot of MW versus logP

plt.figure(figsize=(7.5, 7.5))

sns.scatterplot(x='MW', y = 'LogP', data=df_final, hue='bioactivity_class',
                size='pIC50', edgecolor='black', alpha=0.7)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.savefig('plot_MW_vs_LogP.pdf')


plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_final)

plt.xlabel('Bioactivity class', fontsize = 14, fontweight = 'bold')
plt.ylabel('pIC50 value', fontsize = 14, fontweight = 'bold')

plt.savefig('plot_ic50.pdf')