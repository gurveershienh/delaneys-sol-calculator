################################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################
# importing necessary modules
from tkinter import *
import pandas as pd
import numpy as np
# importing from rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
################################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################


# function for determining aromatic proportion of each mol in sol_data

def aromatic_proportion(molecule):
    
    
    # function for determining number of aromatic atoms in a molecule
    bool_lst = [molecule.GetAtomWithIdx(i).GetIsAromatic() for i in range(molecule.GetNumAtoms())]
    num_aa = []
    for bool in bool_lst:
        if bool:
            num_aa.append(1)
    result = sum(num_aa) / Descriptors.HeavyAtomCount(molecule)
    return result

def solve_descriptors(smiles, verbose=False):
    i=0
    smiles = smiles.replace(' ', '')
    smiles = smiles.split(',')
    desc_data=np.arange(1,1)
    labels = ['logP', 'MolWt', 'NumOfRB', 'AP']
    moldata=[]
    
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)
    
    for mol in moldata:
        
        logP= Descriptors.MolLogP(mol)
        molWt=Descriptors.MolWt(mol)
        numRB=Descriptors.NumRotatableBonds(mol)
        AP = aromatic_proportion(mol)
        
        row = np.array([logP, molWt, numRB, AP])
        if i == 0:
            desc_data=row
        else:
            desc_data=np.vstack([desc_data, row])
        i+=1
    
    print(desc_data,labels)
    descriptors = pd.DataFrame(data=desc_data,columns=labels)
    

print(pd.DataFrame(data=np.array([ 2.1965, 72.151,   2,      0    ]), columns=['logP', 'MolWt', 'NumOfRB', 'AP']))
root = Tk()
root.title('Delaney\'s Solubility Calculator')
root.resizable(width=False,height=False)

canvas = Canvas(root, height=600, width=500)
canvas.pack()

frame = Frame(root, bg='#cdf7dd')
frame.place(relheight=1, relwidth=1)

smiles = Entry(frame)
smiles.place(x=75, y =50, width=350, height=200)

calculate = Button(frame, text='Calculate', padx=10, pady=5, command=lambda:solve_descriptors( smiles.get() ) )
calculate.place(x=200,y=275,width=100, height=50)

lower_frame = Frame(root)
lower_frame.place(x=50,y=350, width=400, height= 225)

output = Label(lower_frame)
output.place(relx=0,rely=0)


root.mainloop()
