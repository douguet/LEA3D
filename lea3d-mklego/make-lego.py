#!/usr/bin/python3.7

import os, sys
import re
import os, shutil
import os.path
from sys import platform

ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: .py <sdf file to read and fragment>\n')
    quit()

# sys.argv[0] is the name of the program itself
filesdf=sys.argv[1]

flagcleanfiles=1;

# sys.argv[0] is the name of the program itself
leaexe=sys.argv[0]
leaexe=re.sub('\/make-lego\.py','',leaexe)

whichexe='linux'
if(whichexe in platform):
    # linux
    leaexec=leaexe + "/lea3d-CLASS_FGT.pl"
    leaexen=leaexe + "/nosel.pl"
    leaexem=leaexe + "/lea3d-MAKE_FGTS.pl"
    leaexef=leaexe + "/filterx.py"
    leaexed=leaexe + "/filter2dsdf.pl"
elif platform == "darwin":
    leaexec=leaexe + "/lea3d-CLASS_FGT.pl"
    leaexen=leaexe + "/nosel.pl"
    leaexem=leaexe + "/lea3d-MAKE_FGTS.pl"
    leaexef=leaexe + "/filterx.py"
    leaexed=leaexe + "/filter2dsdf.pl"
else:
    #windows with conda
    leaexec="perl " + leaexe + "\lea3d-CLASS_FGT.pl"
    leaexen="perl " + leaexe + "\\" + "nosel.pl"
    leaexem="perl " + leaexe + "\lea3d-MAKE_FGTS.pl"
    leaexef="python " + leaexe + "\\" + "filterx.py"
    leaexed="perl " + leaexe + "\\" + "filter2dsdf.pl"

#print("perl %s %s %s %s " % (leaexec,leaexen,leaexem,leaexef))

#######################################
def isexist(file_path):
    return os.path.exists(file_path)

def isexistnotempty(file_path):
    if(os.path.exists(file_path)):
        if(os.path.getsize(file_path)==0):
            fempty=0
        else:
            fempty=1
    else:
        fempty=0
    return fempty

#######################################
#def (filesdf) output=number of molecules
def nbsdf(fsdf):
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    whichend="$$$$"
    nbmol=0
    compt=0
    while(compt < len(getstr)):
        if(whichend in getstr[compt]):
            nbmol=nbmol+1
        compt=compt+1
    if(nbmol==0): #if .mol format
        nbmol=1
    return nbmol

#######################################

#######################################
# MAIN program

nbinput=nbsdf(filesdf)

#clean by removing 2D structures
cmd = '%s %s clean3D.sdf' % (leaexed,filesdf)
os.system(cmd)
if(isexistnotempty("clean3D.sdf")):
    nbclean3d=0
    nbclean3d=nbsdf('clean3D.sdf')

    #clean by dissociating salts
    cmd = '%s %s' % (leaexen,'clean3D.sdf')
    os.system(cmd)
    nbclean=0
    if(isexistnotempty("one.sdf")):
        if(isexistnotempty("dissociated.sdf")):
            #cat both in clean.sdf
            solfile=open('one.sdf', 'r')
            stran=solfile.readlines()
            solfile.close()
            mfile=open("clean.sdf", 'w')
            for f in stran:
                mfile.write(f)
            mfile.close()
            solfile=open('dissociated.sdf', 'r')
            stran=solfile.readlines()
            solfile.close()
            mfile=open("clean.sdf", 'a')
            for f in stran:
                mfile.write(f)
            mfile.close()
            if(isexist("two.sdf")):
                nbclean=nbsdf('two.sdf')
        else:
            shutil.copyfile("one.sdf","clean.sdf")

        #fragment
        cmd = '%s clean.sdf' % (leaexem)
        os.system(cmd)
        nbfgt=nbsdf('make_fgts.sdf')

        #remove fragments without X- dummy atoms
        cmd = '%s make_fgts.sdf 1' % (leaexef)
        os.system(cmd)
        nbfgtclean=nbsdf('filterx.sdf')

        #remove duplicates
        cmd = '%s filterx.sdf key X unique.sdf 0' % (leaexec)
        os.system(cmd)
        shutil.copyfile("unique.sdf","lea3d-legos.sdf")

        nboutput=nbsdf("lea3d-legos.sdf")
        print("\nInput %s molecules (%s with 3D coordinates and %s dissociated because are salt)- %s fragments (%s having X dummy atoms) - remove duplicates - Output in lea3d-legos.sdf %s fragments ready to be used in LEA3D" % (nbinput,nbclean3d,nbclean,nbfgt,nbfgtclean,nboutput))
    else:
        print("Error: empty output file after removing salts")
else:
    print("Error: empty output file after removing 2D structures")

if(flagcleanfiles==1):
    if(isexist("clean3D.sdf")):
        os.remove("clean3D.sdf")
    if(isexist("one.sdf")):
        os.remove("one.sdf")
    if(isexist("dissociated.sdf")):
        os.remove("dissociated.sdf")
    if(isexist("two.sdf")):
        os.remove("two.sdf")
    if(isexist("acyclic.sdf")):
        os.remove("acyclic.sdf")
    if(isexist("fused_rings.sdf")):
        os.remove("fused_rings.sdf")
    if(isexist("linker.sdf")):
        os.remove("linker.sdf")
    if(isexist("ring.sdf")):
        os.remove("ring.sdf")
    if(isexist("substituent.sdf")):
        os.remove("substituent.sdf")
    if(isexist("special.sdf")):
        os.remove("special.sdf")
    if(isexist("clean.sdf")):
        os.remove("clean.sdf")
    if(isexist("make_fgts.sdf")):
        os.remove("make_fgts.sdf")
    if(isexist("filterx.sdf")):
        os.remove("filterx.sdf")
    if(isexist("unique.sdf")):
        os.remove("unique.sdf")
    if(isexist("same_unique.sdf")):
        os.remove("same_unique.sdf")
    if(isexist("key")):
        os.remove("key")

