#!/usr/bin/env python


import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors,Crippen

import numpy as np
import pandas as pd


import os
import sys

from rdkit import Chem
import sys
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
import numpy as np
import argparse, sys, pickle, math, rdkit, matplotlib
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import SimilarityMaps
from numpy.random import randn
import warnings


from rdkit.Chem import Draw
from rdkit.Chem.Fraggle import FraggleSim
from rdkit.Chem.Draw import SimilarityMaps


from rdkit import rdBase, Chem, DataStructs
import pandas as pd
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdFingerprintGenerator

import pandas as pd
from rdkit.Chem import Draw, Lipinski, Crippen, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors as rdmd



file = open('AP_nohash_test.dat', 'w')


data_1 = open(sys.argv[1], 'r')
lines_1 = data_1.readlines()
data_1.close()

data_2 = open(sys.argv[2], 'r')
lines_2 = data_2.readlines()
data_2.close()

counter=0
counter_idx=0

for line in lines_1:
    
    q = Chem.MolFromSmiles(line)
    #fp1 = Pairs.GetAtomPairFingerprint(q)
    counter+=1
    print(counter)
    fp1 = rdmd.GetHashedAtomPairFingerprintAsBitVect(q,1024)
    #fp1 = AllChem.GetMorganFingerprintAsBitVect(q, 3, 4096)
    for line in lines_2:
        counter_idx+=1
        #print(counter_idx)
        m = Chem.MolFromSmiles(line)
        #fp2 = Pairs.GetAtomPairFingerprint(m)
        fp2 = rdmd.GetHashedAtomPairFingerprintAsBitVect(m,1024)
        #fp2 = AllChem.GetMorganFingerprintAsBitVect(q, 3, 4096)
        #print('counter',counter)
        #counter+=1
        #sim,match = FraggleSim.GetFraggleSimilarity(q, m, tverskyThresh=0.8)
        Tversky = DataStructs.TverskySimilarity(fp1, fp2, 0.95, 0.01)
        #print('Tversky',Tversky)
        #print('Fraggle',sim)
        if Tversky >= 0.7:

            print(line)
            #file.write(line)

file.close()
            





