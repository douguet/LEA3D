#!/usr/bin/python3.7

#execute ./<>.py sdf_file nbconf sdf_outputname  

import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
#import numpy as np
#from subprocess import Popen,PIPE

# sys.argv[0] is the name of the program itself
molecule=sys.argv[1]
nbconf=int(sys.argv[2])
output=sys.argv[3]

#######################################
# MAIN program

#m = Chem.MolFromMolFile(molecule)
m = Chem.SDMolSupplier(molecule,removeHs=False)
#remove H
#m = Chem.SDMolSupplier(molecule)
w = Chem.SDWriter(output)

for mi in m:
    #AllChem.EmbedMultipleConfs(ibuH,clearConfs=True,numConfs=100,pruneRmsThresh=1)
    #AllChem.UFFOptimizeMolecule(m,confId=i)
    #AllChem.MMFFOptimizeMolecule(m2)

    #mi = Chem.AddHs(mi)
    
    #to read stereochemistry from 3D coordinates:
    #AssignAtomChiralTagsFromStructure(RDKit::ROMol {lvalue} mol, int confId=-1, bool replaceExistingTags=True)
    Chem.AssignAtomChiralTagsFromStructure(mi)
    #? not tested :
    #Chem.AssignStereochemistry(mi, cleanIt=False, force=False, flagPossibleStereoCenters=True)

    cids = AllChem.EmbedMultipleConfs(mi, numConfs=nbconf, randomSeed=1234,pruneRmsThresh=1)
    #print(len(cids))	
    for mj in cids:
        AllChem.MMFFOptimizeMolecule(mi,confId=mj)
        w.write(mi,confId=mj)

	#rmslist = []
	#AllChem.AlignMolConformers(cids, RMSlist=rmslist)
	#print RMSlist[1]	

w.flush()
w.close()
#print ("Done (see %s)" % output)

#######################################
