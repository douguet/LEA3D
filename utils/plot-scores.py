#!/usr/bin/python3.7

import matplotlib.pyplot as plt
import numpy as np
import os
import re
import sys

ARGV=[]
try:
    ARGV.append(sys.argv[1])
except:
    print('usage: plot-scores.py fitmoy.dat ')
    quit()

# sys.argv[0] is the name of the program itself
filefit=sys.argv[1]

filemol=open(filefit,'r')
get=filemol.read().split('\n')
filemol.close()

tabLignes=[]
compt=0
gen=0
while(compt < len(get)):
    testspace=[]
    testspace.append(re.split('', get[compt]))
    if(testspace[0][1]!=''):
        tabLignes.append(re.split('\s+', get[compt].strip()))
        gen=gen+1
    compt=compt+1

#nb-gen mean-score lowest-score largest-score
arr_gen=np.empty(shape=[gen,1], dtype='int')
arr_mean=np.empty(shape=[gen,1], dtype='float64')
arr_min=np.empty(shape=[gen,1], dtype='float64')
arr_max=np.empty(shape=[gen,1], dtype='float64')
compt2=0
while (compt2 < gen):
    arr_gen[compt2,0]=float(tabLignes[compt2][0])
    arr_mean[compt2,0]=float(tabLignes[compt2][1])
    arr_max[compt2,0]=float(tabLignes[compt2][2])
    arr_min[compt2,0]=float(tabLignes[compt2][3])
    compt2=compt2+1

# plot
fig, ax = plt.subplots()

plt.title("LEA3D results")
ax.plot(arr_gen, arr_mean, linewidth=2.0, label='Mean score')
ax.plot(arr_gen, arr_max, linewidth=2.0, label='Maximun score')
ax.plot(arr_gen, arr_min,linewidth=2.0, label='Minimun score')
ax.legend()

plt.xlabel("Generation number")
plt.ylabel("Fitness score in percentage")

plt.show()


