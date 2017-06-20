from ROOT import *
import numpy
import root_numpy
import time
import os
starttime=time.time()


def jvars(rootfile, first, last):
    jvars=["jet_pt","jet_eta", "jet_phi","jet_mass","jet_flavour"]
    f=TFile(rootfile)
    tree=f.Get("analyzer1/tree")
    t2=root_numpy.tree2array(tree, branches=jvars, selection="(jet_pt>30)&&(abs(jet_eta)<2.4)", start=first, stop=last)
    t2=root_numpy.rec2array(t2)
    numpy.save("jvars_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy", t2)
    print t2.shape
    print time.time()-starttime
    f.Close()
    #os.system("mv "+"jvars_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy"+" /gpfs/ddn/users/lgiannini/NN/IntermediateTargets")
    
def svars(rootfile, first, last):
    stringa=["seed_pt","seed_eta","seed_phi","seed_mass","seed_dz","seed_dxy",
             "seed_3D_ip","seed_3D_sip","seed_2D_ip","seed_2D_sip","seed_3D_signedIp","seed_3D_signedSip","seed_2D_signedIp","seed_2D_signedSip",
             "seed_chi2reduced","seed_nPixelHits","seed_nHits","seed_jetAxisDistance","seed_jetAxisDlength"    ]
    f=TFile(rootfile)
    tree=f.Get("analyzer1/tree")    
    t2=root_numpy.tree2array(tree, branches=stringa, selection="(jet_pt>30)&&(abs(jet_eta)<2.4)", start=first, stop=last)
    ll=len(t2)
    t2=root_numpy.rec2array(t2)
    print t2.shape
    t2=t2.reshape((10,len(stringa),ll))
    t2=t2.swapaxes(0, 2)
#    t2=numpy.reshape(t2, (len(t2), len(stringa), 10))
#    print t2.shape
    numpy.save("svars"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy", t2)
    print time.time()-starttime
    print t2.shape

    f.Close()
    #os.system("mv "+"svars"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy"+" /gpfs/ddn/users/lgiannini/NN/IntermediateTargets")

def tvars(rootfile, first, last):
    stringa=["seed_pt","seed_eta","seed_phi","seed_mass","seed_dz","seed_dxy",
             "seed_3D_ip","seed_3D_sip","seed_2D_ip","seed_2D_sip","seed_3D_signedIp","seed_3D_signedSip","seed_2D_signedIp","seed_2D_signedSip",
             "seed_chi2reduced","seed_nPixelHits","seed_nHits","seed_jetAxisDistance","seed_jetAxisDlength",
             
             #more truth vars
             "seed_MC_pt",
            "seed_MC_eta",
            "seed_MC_phi",
            "seed_MC_mass",
            "seed_MC_dz",
            "seed_MC_dxy",
            "seed_MC_MomPdgId",
            "seed_MC_MomFlavour",
            "seed_MC_BChain",
            "seed_MC_DChain", #no B in this case
            "seed_MC_vx",
            "seed_MC_vy",
            "seed_MC_vz",
            "seed_MC_pvd"             
             
             ]
    
    
    stringa2=["nearTracks_pt","nearTracks_eta","nearTracks_phi","nearTracks_dz","nearTracks_dxy","nearTracks_mass","nearTracks_3D_ip","nearTracks_3D_sip",
          "nearTracks_2D_ip","nearTracks_2D_sip","nearTracks_PCAdist","nearTracks_PCAdsig","nearTracks_PCAonSeed_x","nearTracks_PCAonSeed_y","nearTracks_PCAonSeed_z",
          "nearTracks_PCAonSeed_xerr","nearTracks_PCAonSeed_yerr","nearTracks_PCAonSeed_zerr","nearTracks_PCAonTrack_x","nearTracks_PCAonTrack_y","nearTracks_PCAonTrack_z",
          "nearTracks_PCAonTrack_xerr","nearTracks_PCAonTrack_yerr","nearTracks_PCAonTrack_zerr","nearTracks_dotprodTrack","nearTracks_dotprodSeed","nearTracks_dotprodTrackSeed2D",
          "nearTracks_dotprodTrackSeed3D","nearTracks_dotprodTrackSeedVectors2D","nearTracks_dotprodTrackSeedVectors3D","nearTracks_PCAonSeed_pvd","nearTracks_PCAonTrack_pvd",
          "nearTracks_PCAjetAxis_dist","nearTracks_PCAjetMomenta_dotprod","nearTracks_PCAjetDirs_DEta","nearTracks_PCAjetDirs_DPhi",
          
          #more MC vars
          
            "nearTracks_MC_pt",
            "nearTracks_MC_eta",
            "nearTracks_MC_phi",
            "nearTracks_MC_dz",
            "nearTracks_MC_dxy",

            "nearTracks_MC_MomPdgId",
            "nearTracks_MC_MomFlavour",
            "nearTracks_MC_BChain",
            "nearTracks_MC_DChain", #no B in this case

            "nearTracks_MC_Track_vx",
            "nearTracks_MC_Track_vy",
            "nearTracks_MC_Track_vz",

            "nearTracks_MC_fromSeedVtx",
            "nearTracks_MC_fromSeedChain",
            "nearTracks_MC_pvd",
            #"nearTracks_MC_fromSeedVtx*(nearTracks_MC_pvd>0)"

          
          
          ]

    f=TFile(rootfile)
#    tree=f.Get("analyzer1/tree")
    tree=root_numpy.tree2array(f.Get('analyzer1/tree'),branches=stringa2, selection="(jet_pt>30)&&(abs(jet_eta)<2.4)", start=first, stop=last)
    print "loaded"
    tree2=root_numpy.rec2array(tree)
    print "s"
    print tree2.shape
    print round(time.time()-starttime,2), "reshape"
    tree3=tree2.reshape((200,51,len(tree)))
    tree3=tree3.reshape((10,51*20,len(tree)))
    print tree3.shape
    tree3=tree3.swapaxes(0, 2)
    t2=root_numpy.tree2array(f.Get('analyzer1/tree'), branches=stringa, selection="(jet_pt>30)&&(abs(jet_eta)<2.4)", start=first, stop=last)
    t2=root_numpy.rec2array(t2)
    print t2.shape
    t2=t2.reshape((10,len(stringa),len(tree)))
    t2=t2.swapaxes(0, 2)
    tree5=numpy.concatenate((t2, tree3), axis=1)
    print tree5.shape
    numpy.save("tvars_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy", tree5)
    print time.time()-starttime
    f.Close() 
    os.system("mv "+"tvars_"+str(first)+"_"+str(last)+"_"+rootfile.split(".")[0]+".npy"+" /gpfs/ddn/users/lgiannini/NN/DataRECOValidation")


import sys
rootfile=sys.argv[1]
start=int(sys.argv[2])
stop=int(sys.argv[3])
#jvars(rootfile, start, stop)
#svars(rootfile, start, stop)
tvars(rootfile, start, stop)
