import sys

files=[]

f=open("/scratch/lgiannini/test_class/CMSSW_8_1_0_pre16/src/seedAnalyzer/seedAnalyzer/test/All_miniAOD.txt", "r")
for line in f:
    files.append('root://xrootd-cms.infn.it/'+line.replace('\n', ''))


print files[int(sys.argv[2]):(int(sys.argv[3]))]
print len(files)


import FWCore.ParameterSet.Config as cms
process = cms.Process("Pippo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'


#process.options = cms.untracked.PSet()

#process.options.numberOfThreads = cms.untracked.uint32(1)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        files[int(sys.argv[2]):(int(sys.argv[3]))]
#    "file:021021EC-E0FD-E611-BBAD-0025905B85BA.root",

)#,
    #skipEvents=cms.untracked.uint32(36590)
)

#for MC jet flavour
#process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
#process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
#process.ak4JetFlavourInfos.jets = cms.InputTag("ak4PFJetsCHS")
#process.ak4JetFlavourInfos.hadronFlavourHasPriority = cms.bool(True)
#process.flavourSeq = cms.Sequence(
    #process.selectedHadronsAndPartons *
    #process.ak4JetFlavourInfos
#)


process.analyzer1 = cms.EDAnalyzer('AnalyzerSignedIP_MINIAOD_wArr',
 #jetMCSrc = cms.InputTag("ak4JetFlavourInfos"),
#    bits = cms.InputTag("TriggerResults","","HLT"),
#    prescales = cms.InputTag("patTrigger"),
#    objects = cms.InputTag("selectedPatTrigger"),
#    generatorInfo = cms.InputTag("generator"),
)

outfile='miniaodFiles/NewMiniArrays'+sys.argv[2]+'.root'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outfile) )


process.p = cms.Path(process.analyzer1)


