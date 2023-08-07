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
plt.hist(ic50_values)
plt.show()