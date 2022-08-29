#!/usr/bin/env python3

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt


smi = 'CC(=CC(C)(C)CCCCCCCC(=O)O)C1CCC2C(=CC=C3CC(O)CC(O)C3)CCCC21C'
molecule = Chem.MolFromSmiles(smi)
fig = Draw.MolToMPL(molecule)

plt.title('Molecule')

fig.savefig('mol.jpeg', bbox_inches='tight')
 
