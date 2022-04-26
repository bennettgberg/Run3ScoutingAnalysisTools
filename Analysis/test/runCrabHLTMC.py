from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'Run3Scouting_DYToLL_Run3Winter21_crabRunSkim11'
config.General.workArea = 'crab_Run3Scouting_DYToLL_Run3Winter21'
config.General.transferOutputs = True
#config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ScoutingNanoAOD_cfg.py'

#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-RECO'
config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220201_DYToLLM50Run3Summer21_asahasra/'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.totalUnits = 20000
config.Data.publication = False
#config.Data.outputDatasetTag = 'Scouting_DYToLL_Run3Winter21_crabRunSkim9_asahasra'
config.Data.outputDatasetTag = 'Scouting_DYToLL_Run3Winter21_crabRunSkim11_bgreenbe'
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_BE_IIHE'
