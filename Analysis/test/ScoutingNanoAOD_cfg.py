import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC',
    True,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterMuons',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to contain at least two muons'
)

params.register(
    'reducedInfo',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to store just the reduced information'
)

params.register(
    'trigProcess',
    'HLT',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData', 
    '112X_mcRun3_2021_realistic_v16', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagMC',
    'auto:phase1_2021_realistic',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'xsec',
    0.001,
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)

params.register(
    'inputFile', 
    'HLT2022_HLT_OR_GENSIMDIGIRAW.root', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Name of the input root file'
)

params.register(
    'output', 
    'scoutingNTuple.root', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Name of the output root file'
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
                                #params.inputFile#'file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220127_DYToLLM50Run3Summer21_asahasra/220127_135957/0000/HLT2022_HLT_1.root'
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_1.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_2.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_3.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_4.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_5.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_6.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_7.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_8.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_9.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_10.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_11.root',
                                'file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_12.root',
                            ])
                        )

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC :
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(params.output)
                               )
#from DarkPhotonAnalysis.DimuonAnalysis2018.TriggerPaths_cfi import getL1Conf
L1Info = ['L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu_12_5','L1_DoubleMu_15_7','L1_TripleMu_5_3_3','L1_TripleMu_5_5_3','L1_QuadMu0','L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4','L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18','L1_DoubleMu4_SQ_OS_dR_Max1p2','L1_SingleMu22','L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4','L1_DoubleMu4p5_SQ_OS_dR_Max1p2','L1_DoubleMu4p5_SQ_OS','L1_DoubleMu0er1p5_SQ_dR_Max1p4','L1_DoubleMu0er2p0_SQ_dR_Max1p4','L1_DoubleMu0_SQ']
# Make tree
process.mmtree = cms.EDAnalyzer('EGammaOnly_ScoutingNanoAOD',
                                triggerresults   = cms.InputTag("TriggerResults", "", params.trigProcess),
                                doL1 = cms.bool(False),
                                triggerConfiguration = cms.PSet(
                                    hltResults            = cms.InputTag('TriggerResults','','HLT'),
                                    l1tResults            = cms.InputTag(''),
                                    daqPartitions         = cms.uint32(1),
                                    l1tIgnoreMaskAndPrescale = cms.bool(False),
                                    throw                 = cms.bool(False)
                                ),
                                ReadPrescalesFromFile = cms.bool( False ),
                                AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1Seeds           = cms.vstring(L1Info),
                                beamspot        = cms.InputTag("hltOnlineBeamSpot"),
                                electrons        = cms.InputTag("hltScoutingEgammaPacker"),
                                photons          = cms.InputTag("hltScoutingEgammaPacker"),
                                gens = cms.InputTag("genParticles"),
                                hltEgammaExtra = cms.InputTag("hltEgammaHLTExtra"),
                                
                            )

process.p = cms.Path(                  process.mmtree)
