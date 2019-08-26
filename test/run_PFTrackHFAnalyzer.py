#------------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import FWCore.ParameterSet.VarParsing as VarParsing

#------------------------------------------------------------------------------------
# Declare the process and input variables
#------------------------------------------------------------------------------------
#process = cms.Process('NOISE',eras.Run2_50ns)#for 50ns 13 TeV data
#process = cms.Process('NOISE',eras.Run2_25ns)#for 25ns 13 TeV data
options = VarParsing.VarParsing ('analysis')
process = cms.Process("Trees",eras.Phase2) 

##
## Setup command line options
##
options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")

##
## Default
##
options.maxEvents = -1 # means all events
#options.skipEvents = 0 # default is 0.

##
## get and parse the command line arguments
##
options.parseArguments()

#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# TTbar sample
#
# RECO
options.inputFiles = [
#    'file:/afs/cern.ch/user/i/iashvili/step3_SinglePiPt25Eta2p8_4p2.root',
#
# noTrackLink samples 
#
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt10Eta3p45_3p54_n1000_noTrackLink_step3.root',
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiMinusPt10Eta3p45_3p54_n1000_noTrackLink_step3.root',     
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta2p8_4p2_n500_noTrackLink_step3.root',
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta3p45_3p54_n1000_noTrackLink_step3.root',
#
# TrackLink samples 
#
    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta3p45_3p54_n800_step3.root',
#
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt10Eta3p45_3p54_n1000_FakeTracks.root',
#
]

options.outputFile = 'hist_PFTrackHFAnalyzer.root'
#
#
#
print("maxEvents: ", options.maxEvents)
print("inputFiles: ", options.inputFiles)
print("outputFile: ", options.outputFile)

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(options.skipEvents) # default is 0.
)
readFiles.extend(options.inputFiles);

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

#------------------------------------------------------------------------------------
# import of standard configurations
#------------------------------------------------------------------------------------
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

#------------------------------------------------------------------------------------
# Set up our analyzer
#------------------------------------------------------------------------------------

#process.load("PFAnalysis.PFAnalyzers.PFTrackHFAnalyzer_cfi")
eventN = [1,3,4]
process.pfTrackHFAnalyzer = cms.EDAnalyzer("PFTrackHFAnalyzer",
                                           source_genpars = cms.untracked.InputTag('genParticles', ''),
                                           source_calopars = cms.untracked.InputTag('mix', 'MergedCaloTruth'),
                                           source_vertices = cms.untracked.InputTag('offlinePrimaryVertices', ''),
                                           source_pfcands = cms.untracked.InputTag('particleFlow', ''),
                                           source_pfclustersHF = cms.untracked.InputTag('particleFlowClusterHF', ''),
                                           source_pfrechitsHF = cms.untracked.InputTag('particleFlowRecHitHF', ''),
                                           source_pftracks = cms.untracked.InputTag('pfTrack', ''),
                                           source_tracks = cms.untracked.InputTag('generalTracks', ''),
                                           source_hfrechits = cms.untracked.InputTag('hfreco', ''),
                                           debug = cms.untracked.bool(True),
                                           debugRecHit = cms.untracked.bool(False),
                                           EventsToScan = cms.vint32(eventN),
                                           ptlow  = cms.double(0.),
                                           pthigh = cms.double(1000.),    
                                           etalow  = cms.double(0.),    
                                           etahigh = cms.double(10.),

)

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#------------------------------------------------------------------------------------
# Sequence definition
#------------------------------------------------------------------------------------
process.ana_step = cms.Sequence(
    process.pfTrackHFAnalyzer
)

#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.ana_step
)
