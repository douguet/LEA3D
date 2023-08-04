#!/usr/bin/python3.7

#execute ./<>.py sdf_file nbconf sdf_outputname  

import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

# sys.argv[0] is the name of the program itself
molecule=sys.argv[1]

#######################################
# MAIN program

#hydrogens as they are:
#m = Chem.MolFromMolFile(molecule)
m = Chem.SDMolSupplier(molecule,removeHs=False)
#remove+add H
#m = Chem.SDMolSupplier(molecule)

for mi in m:
    logp = Descriptors.MolLogP(mi)
    print ("MolLogP %4.2f" % logp)

#######################################
