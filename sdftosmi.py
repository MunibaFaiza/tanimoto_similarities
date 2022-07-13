#!/usr/bin/env python3

import sys
from rdkit import Chem

def converter(file_name):
    sdf_file = Chem.SDMolSupplier(file_name)
    
    out_file = open('smiles.txt', "w")
    for mol in sdf_file:
        if mol is not None:             # avoiding compounds that cannot be loaded.
            smi = Chem.MolToSmiles(mol)
            out_file.write(f"{smi}\n")
    out_file.close()
if __name__ == "__main__":
    converter(sys.argv[1])
