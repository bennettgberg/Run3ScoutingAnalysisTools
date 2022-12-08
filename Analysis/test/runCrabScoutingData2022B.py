from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3Scouting2022B_crabRunSkim220823'
config.General.workArea = 'crabRunSkim220823_Run3Scouting2022B'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ScoutingNanoAOD_forCrab_cfg.py'
config.JobType.inputFiles = ['Cert_Collisions2022_355100_357101_Golden.json']

config.Data.inputDataset = '/ScoutingPFRun3/Run2022B-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.LumiMask = 'Cert_Collisions2022_355100_357101_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = 'Scouting2022B_crabRunSkim220823_asahasra'
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_BE_IIHE'
