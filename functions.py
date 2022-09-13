# importing necessary modules

import pandas as pd
import numpy as np

# importing from rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

# creating pandas dataframe for solubility data

sol_data = pd.read_csv('delaney-sol.csv')
mol_list = [Chem.MolFromSmiles(element) for element in sol_data.SMILES]

# function for determining logP MolWt numRB descriptors from sol_data

def solve_descriptors(smiles, verbose=False):
    '''
    
    '''
    i = 0
    desc_data=np.arange(1,1)
    labels = ['logP', 'MolWt', 'NumOfRB']
    moldata=[Chem.MolFromSmiles(elem) for elem in smiles]
    
    for mol in moldata:
        
        logP= Descriptors.MolLogP(mol)
        molWt=Descriptors.MolWt(mol)
        numRB=Descriptors.NumRotatableBonds(mol)
        
        row = np.array([logP, molWt, numRB])
        
        if i == 0:
            desc_data=row
        else:
            desc_data=np.vstack([desc_data, row])
        i+=1
        
    descriptors = pd.DataFrame(data=desc_data,columns=labels)
    
    return descriptors
    
    

# function for determining aromatic proportion of each mol in sol_data

def aromatic_proportion(smiles):
    
    moldata=[Chem.MolFromSmiles(elem) for elem in smiles]
    
    # function for determining number of aromatic atoms in a molecule
    def aromatic_atoms(molecule):
        bool_lst = [molecule.GetAtomWithIdx(i).GetIsAromatic() for i in range(molecule.GetNumAtoms())]
        count = []
        for bool in bool_lst:
            if bool:
                count.append(1)
        result = sum(count)
        return result
    
    desc_data = [aromatic_atoms(mol)/Descriptors.HeavyAtomCount(mol) for mol in moldata]
    AP_descriptors = pd.DataFrame(desc_data, columns=['AP'])
    
    return AP_descriptors


# creating X matrix

descriptors = solve_descriptors(sol_data.SMILES)
descriptors_ap = aromatic_proportion(sol_data.SMILES)

df = pd.concat([descriptors, descriptors_ap], axis=1)
print(df)

##run df.to_csv('delaney-descriptors.csv', index=False) to produce csv file of descriptors






