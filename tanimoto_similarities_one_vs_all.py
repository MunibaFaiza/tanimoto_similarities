#!/usr/bin/env python3

import time
import random
import sys
from pathlib import Path
import seaborn as sns

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import SimilarityMaps

# show full results
np.set_printoptions(threshold=sys.maxsize)


# Defining query smile
query = Chem.MolFromSmiles('CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C')
fp1 = Chem.RDKFingerprint(query)

# Reading sdf of other compounds
suppl = Chem.SDMolSupplier('lig.sdf')
for mol in suppl:
    with open('sims_vs_all.txt', 'a') as outfile:
        fp2 = Chem.RDKFingerprint(mol)
        sim = DataStructs.TanimotoSimilarity(fp1,fp2)
        print(sim, file=outfile)