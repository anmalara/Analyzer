from ROOT import *
import numpy
import root_numpy
import time
import os
starttime=time.time()


def jvars(rootfile, first, last):
    jvars=["jet_pt","jet_eta", "jet_phi","jet_mass","jet_flavour", "jet_CSVv2", "jet_CMVA", "jet_TCHE", "jet_SSVHE"]
    f=TFile(rootfile)
    tree=f.Get("analyzer1/tree")
    t2=root_numpy.tree2array(tree, branches=jvars, selection="(jet_pt>30)&&(abs(jet_eta)<2.4)", start=first, stop=last)
    t2=root_numpy.rec2array(t2)
    numpy.save("jvarsNew_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy", t2)
    print t2.shape
    print time.time()-starttime
    f.Close()
    os.system("mv "+"jvarsNew_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy"+" /gpfs/ddn/users/lgiannini/NN/DataMiniAODNewValidation")
    

import sys
rootfile=sys.argv[1]
start=int(sys.argv[2])
stop=int(sys.argv[3])
jvars(rootfile, start, stop)

