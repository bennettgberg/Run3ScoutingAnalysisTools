// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/TrackReco/interface/fillCovariance.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

using namespace std;


class EGammaOnly_ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet&);
  ~EGammaOnly_ScoutingNanoAOD();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;
  
  const edm::EDGetTokenT<reco::BeamSpot>             	beamSpotToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  	electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >  	photonsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> > gensToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

        
	
  bool doL1;       
  triggerExpression::Data triggerCache_;
      
	
  // Generator-level information
  // Flags for the different types of triggers used in the analysis
  // For now we are interested in events passing either the single or double lepton triggers
  unsigned char                trig;
       
  edm::InputTag                algInputTag_;       
  edm::InputTag                extInputTag_;       
  edm::EDGetToken              algToken_;
  //l1t::L1TGlobalUtil          *l1GtUtils_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
       
  Float16_t beamspot_x;
  Float16_t beamspot_y;
  Float16_t beamspot_z;
        
  // Gen Level Lepton, Neutrino and DM particle
  UInt_t n_gen;
  vector<Int_t> genpart_pdg;
  vector<Float16_t> genpart_pt;
  vector<Float16_t> genpart_eta;
  vector<Float16_t> genpart_phi;
  vector<Float16_t> genpart_m;
  vector<Float16_t> genpart_vx;
  vector<Float16_t> genpart_vy;
  vector<Float16_t> genpart_vz;
  vector<Int_t> genpart_nmoms;
  vector<Int_t> genpart_mompdg;
  vector<Bool_t> genpart_isPromptFS;
  vector<Bool_t> genpart_isPromptDec;
  vector<Bool_t> genpart_isDirectPromptTauDecayProdFS;

  UInt_t n_genpartmomZ;
  vector<Int_t> genpartmomZ_pdg;
  vector<Float16_t> genpartmomZ_pt;
  vector<Float16_t> genpartmomZ_eta;
  vector<Float16_t> genpartmomZ_phi;
  vector<Float16_t> genpartmomZ_m;
  vector<Float16_t> genpartmomZ_vx;
  vector<Float16_t> genpartmomZ_vy;
  vector<Float16_t> genpartmomZ_vz;
  vector<Bool_t> genpartmomZ_isPromptFS;
  vector<Bool_t> genpartmomZ_isPromptDec;
  vector<Bool_t> genpartmomZ_isDirectPromptTauDecayProdFS;

  //Electron
  const static int 	max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> 	    Electron_pt;
  vector<Float16_t>        Electron_eta;
  vector<Float16_t>        Electron_phi;
  vector<Float16_t>	    Electron_m;
  vector<Float16_t>        Electron_d0;
  vector<Float16_t>	    Electron_dz;
  vector<Float16_t>	    Electron_detain;
  vector<Float16_t>	    Electron_dphiin;
  vector<Float16_t>	    Electron_sigmaietaieta;
  vector<Float16_t>	    Electron_hoe;
  vector<Float16_t>	    Electron_ooemoop;
  vector<Int_t>	    Electron_missinghits;
  vector<Int_t>        Electron_charge;
  vector<Float16_t>        Electron_ecaliso;
  vector<Float16_t>	    Electron_hcaliso;
  vector<Float16_t>        Electron_tkiso;
  vector<Float16_t>        Electron_r9;
  vector<Float16_t>        Electron_smin;
  vector<Float16_t>        Electron_smaj;
  vector<UInt_t>        Electron_seedid;
  vector<vector<Float16_t>> Electron_energymatrix;
  vector<vector<Float16_t>> Electron_timingmatrix;

  //trying to get Electron supercluster energy (should be the same as E in H/E)
  vector<Float16_t>   Electron_energy;
  //Photon
  const static int 	max_pho = 1000;
  UInt_t n_pho;
  vector<Float16_t> 	    	Photon_pt;
  vector<Float16_t>        	Photon_eta;
  vector<Float16_t>        	Photon_phi;
  vector<Float16_t>	    	Photon_m;
  vector<Float16_t>	    	Photon_sigmaietaieta;
  vector<Float16_t>	    	Photon_hoe;
  vector<Float16_t>        	Photon_ecaliso;
  vector<Float16_t>	    	Photon_hcaliso;
  vector<Float16_t>	    	Photon_trkiso;
  vector<Float16_t>	    	Photon_r9;
  vector<Float16_t>	    	Photon_smin;
  vector<Float16_t>	    	Photon_smaj;
  vector<Float16_t>	    	Photon_seedid;
  vector<vector<Float16_t>> Photon_energymatrix;
  vector<vector<Float16_t>> Photon_timingmatrix;

  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;

};

EGammaOnly_ScoutingNanoAOD::EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
  beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken           (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))), 
  gensToken                (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("gens"))),
  doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
  usesResource("TFileService");
  if (doL1) {
   algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
   extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
   algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
   l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    /* l1GtUtils_ = new l1t::L1TGlobalUtil(iConfig,consumesCollector());*/	
   l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(
    iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree"       , "tree");

  // Event weights
    
  tree->Branch("lumSec"		, &lumSec			 , "lumSec/i" );
  tree->Branch("run"			, &run				 , "run/i" );
    
  // Triggers
  tree->Branch("trig"                 , &trig                          , "trig/b");
  tree->Branch("l1Result"		, "std::vector<bool>"             ,&l1Result_	, 32000, 0);		
  // Pileup info
  //tree->Branch("nvtx"                 , &nvtx                          , "nvtx/i"       );

  tree->Branch("beamspot_x", &beamspot_x, "beamspot_x/f");
  tree->Branch("beamspot_y", &beamspot_y, "beamspot_y/f");
  tree->Branch("beamspot_z", &beamspot_z, "beamspot_z/f");

  // Gen level particles
  tree->Branch("n_genpart", &n_gen, "n_genpart/i");
  tree->Branch("genpart_pdg", &genpart_pdg);
  tree->Branch("genpart_pt", &genpart_pt);
  tree->Branch("genpart_eta", &genpart_eta);
  tree->Branch("genpart_phi", &genpart_phi);
  tree->Branch("genpart_m", &genpart_m);
  tree->Branch("genpart_vx", &genpart_vx);
  tree->Branch("genpart_vy", &genpart_vy);
  tree->Branch("genpart_vz", &genpart_vz);
  tree->Branch("genpart_nmoms", &genpart_nmoms);
  tree->Branch("genpart_mompdg", &genpart_mompdg);
  tree->Branch("genpart_isPromptFS", &genpart_isPromptFS);
  tree->Branch("genpart_isPromptDec", &genpart_isPromptDec);
  tree->Branch("genpart_isDirectPromptTauDecayProdFS", &genpart_isDirectPromptTauDecayProdFS);

  tree->Branch("n_genpartmomZ", &n_genpartmomZ, "n_genpartmomZ/i");
  tree->Branch("genpartmomZ_pdg", &genpartmomZ_pdg);
  tree->Branch("genpartmomZ_pt", &genpartmomZ_pt);
  tree->Branch("genpartmomZ_eta", &genpartmomZ_eta);
  tree->Branch("genpartmomZ_phi", &genpartmomZ_phi);
  tree->Branch("genpartmomZ_m", &genpartmomZ_m);
  tree->Branch("genpartmomZ_vx", &genpartmomZ_vx);
  tree->Branch("genpartmomZ_vy", &genpartmomZ_vy);
  tree->Branch("genpartmomZ_vz", &genpartmomZ_vz);
  tree->Branch("genpartmomZ_isPromptFS", &genpartmomZ_isPromptFS);
  tree->Branch("genpartmomZ_isPromptDec", &genpartmomZ_isPromptDec);
  tree->Branch("genpartmomZ_isDirectPromptTauDecayProdFS", &genpartmomZ_isDirectPromptTauDecayProdFS);

  //Electrons
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  tree->Branch("Electron_m", &Electron_m);
  tree->Branch("Electron_d0", &Electron_d0);
  tree->Branch("Electron_dz", &Electron_dz);
  tree->Branch("Electron_detain", &Electron_detain);
  tree->Branch("Electron_dphiin", &Electron_dphiin);
  tree->Branch("Electron_sigmaietaieta", &Electron_sigmaietaieta);
  tree->Branch("Electron_hoe", &Electron_hoe);
  tree->Branch("Electron_ooemoop", &Electron_ooemoop);
  tree->Branch("Electron_missinghits", &Electron_missinghits);
  tree->Branch("Electron_charge", &Electron_charge);
  tree->Branch("Electron_ecaliso", &Electron_ecaliso);
  tree->Branch("Electron_hcaliso", &Electron_hcaliso);
  tree->Branch("Electron_tkiso", &Electron_tkiso);
  tree->Branch("Electron_r9", &Electron_r9);
  tree->Branch("Electron_smin", &Electron_smaj);
  tree->Branch("Electron_smaj", &Electron_smin);
  tree->Branch("Electron_seedid", &Electron_seedid);
  tree->Branch("Electron_energymatrix", &Electron_energymatrix);
  tree->Branch("Electron_timingmatrix", &Electron_timingmatrix);
  
  //Photons
  tree->Branch("n_pho", &n_pho, "n_pho/i");
  tree->Branch("Photon_pt", &Photon_pt);
  tree->Branch("Photon_eta", &Photon_eta);
  tree->Branch("Photon_phi", &Photon_phi);	
  tree->Branch("Photon_m", &Photon_m);
  tree->Branch("Photon_sigmaietaieta", &Photon_sigmaietaieta);
  tree->Branch("Photon_hoe", &Photon_hoe);
  tree->Branch("Photon_ecaliso", &Photon_ecaliso);
  tree->Branch("Photon_hcaliso", &Photon_hcaliso);
  tree->Branch("Photon_trkiso", &Photon_trkiso);
  tree->Branch("Photon_r9", &Photon_r9);
  tree->Branch("Photon_smin", &Photon_smin);
  tree->Branch("Photon_smaj", &Photon_smaj);
  tree->Branch("Photon_seedid", &Photon_seedid);
  tree->Branch("Photon_energymatrix", &Photon_energymatrix);
  tree->Branch("Photon_timingmatrix", &Photon_timingmatrix);
}


EGammaOnly_ScoutingNanoAOD::~EGammaOnly_ScoutingNanoAOD() {
}

void EGammaOnly_ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;
      
  //cout << "about to try to get the handles." << endl;
  // Handles to the EDM content
  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);
    
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);

  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  bool eleValid = electronsH.isValid();
  bool beamValid = beamSpotH.isValid();
  //if(!eleValid) {
  //  cout << "Error: electronsH invalid!!!" << endl;
  //  if(beamValid) cout << "beamSpotH is valid tho" << endl;
  //  else cout << "beamSpotH is invalid too." << endl;
  //}  
  //else {
  //  cout << "electronsH is valid." << endl;
  //}

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  bool phoValid = photonsH.isValid();

  Handle<vector<reco::GenParticle> >gensH;
  iEvent.getByToken(gensToken, gensH);

  //Handle<Float_t>gensT0H;
  //iEvent.getByToken(gensT0Token, gensT0H);

  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();
  //cout<<"Run number: "<<run<<". Lumi: "<<lumSec<<endl;
  //cout << "handles and junk gotten." << endl; 
  // Which triggers fired
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
    if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   1; // DST_L1DoubleMu_CaloScouting_PFScouting
    if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   2; // DST_DoubleMu3_Mass10_CaloScouting_PFScouting
    if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   4; // DST_ZeroBias_CaloScouting_PFScouting
    if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   8; // DST_L1HTT_CaloScouting_PFScouting
    if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  16; // DST_CaloJet40_CaloScouting_PFScouting
    if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  32; // DST_HT250_CaloScouting
    if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  64; // DST_HT410_PFScouting
    if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig += 128; // DST_HT450_PFScouting
  }

  beamspot_x = beamSpotH->x0();
  beamspot_y = beamSpotH->y0();
  beamspot_z = beamSpotH->z0();
  
//  /////////////////// for genmatching /////////////////////////////////////////////
  n_gen=0;
  n_genpartmomZ = 0;
  for (auto gen_iter = gensH->begin(); gen_iter != gensH->end(); ++gen_iter) {
    if((std::abs(gen_iter->pdgId())==11 || std::abs(gen_iter->pdgId())==13 || std::abs(gen_iter->pdgId())==15)) {
      genpart_pdg.push_back(gen_iter->pdgId());
      genpart_pt.push_back(gen_iter->pt());
      genpart_eta.push_back(gen_iter->eta());
      genpart_phi.push_back(gen_iter->phi());
      genpart_m.push_back(gen_iter->mass());
      genpart_vx.push_back(gen_iter->vx());
      genpart_vy.push_back(gen_iter->vy());
      genpart_vz.push_back(gen_iter->vz());
      genpart_nmoms.push_back(gen_iter->numberOfMothers());
      genpart_mompdg.push_back(gen_iter->mother(0)->pdgId());
      genpart_isPromptFS.push_back(gen_iter->isPromptFinalState());
      genpart_isPromptDec.push_back(gen_iter->isPromptDecayed());
      genpart_isDirectPromptTauDecayProdFS.push_back(gen_iter->isDirectPromptTauDecayProductFinalState());
      n_gen++;
    }
    //if((std::abs(gen_iter->pdgId())==11 || std::abs(gen_iter->pdgId())==13 || std::abs(gen_iter->pdgId())==15) && gen_iter->isLastCopy() ){ 
    //std::cout<<n_gen<<"\t"<<gen_iter->pdgId()<<"\t"<<gen_iter->status()<<"\t"<<gen_iter->pt()<<std::endl;
    //std::cout<<gen_iter->pdgId()<<"\t"<<gen_iter->status()<<"\t"<<gen_iter->pt()<<"\t"<<gen_iter->numberOfDaughters()<<"\t"<<gen_iter->numberOfMothers()<<"\t"<<gen_iter->fromHardProcessBeforeFSR()<<"\t"<<gen_iter->fromHardProcessDecayed()<<"\t"<<gen_iter->fromHardProcessFinalState()<<std::endl;
    //}
    if( ((std::abs(gen_iter->pdgId())==11 || std::abs(gen_iter->pdgId())==13 || std::abs(gen_iter->pdgId())==15)) && 
	(std::abs(gen_iter->mother(0)->pdgId())==23 || std::abs(gen_iter->mother(0)->pdgId())==22)) {
      genpartmomZ_pdg.push_back(gen_iter->pdgId());
      genpartmomZ_pt.push_back(gen_iter->pt());
      genpartmomZ_eta.push_back(gen_iter->eta());
      genpartmomZ_phi.push_back(gen_iter->phi());
      genpartmomZ_m.push_back(gen_iter->mass());
      genpartmomZ_vx.push_back(gen_iter->vx());
      genpartmomZ_vy.push_back(gen_iter->vy());
      genpartmomZ_vz.push_back(gen_iter->vz());      
      genpartmomZ_isPromptFS.push_back(gen_iter->isPromptFinalState());
      genpartmomZ_isPromptDec.push_back(gen_iter->isPromptDecayed());
      genpartmomZ_isDirectPromptTauDecayProdFS.push_back(gen_iter->isDirectPromptTauDecayProductFinalState());
      n_genpartmomZ++;
    }
  }
 ////////////////////// for gen matching (I think) ////////////////////////////////////////// 
  n_ele = 0;
  //cout << "now about to start the electrons loop." << endl;
  //for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
  if(eleValid)
  for (auto &ele : *electronsH) {
    auto *electrons_iter = &ele;
    //cout << "started the loop." << endl;
    Electron_pt.push_back(electrons_iter->pt());
    Electron_eta.push_back(electrons_iter->eta());
    Electron_phi.push_back(electrons_iter->phi());	
    Electron_m.push_back(electrons_iter->m());
    Electron_d0.push_back(electrons_iter->d0());
    Electron_dz.push_back(electrons_iter->dz());
    Electron_detain.push_back(electrons_iter->dEtaIn());
    Electron_dphiin.push_back(electrons_iter->dPhiIn());
    Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
    Electron_hoe.push_back(electrons_iter->hOverE());	
    Electron_ooemoop.push_back(electrons_iter->ooEMOop());
    Electron_missinghits.push_back(electrons_iter->missingHits());
    Electron_charge.push_back(electrons_iter->charge());
    Electron_ecaliso.push_back(electrons_iter->ecalIso());
    Electron_hcaliso.push_back(electrons_iter->hcalIso());
    Electron_tkiso.push_back(electrons_iter->trackIso());
    Electron_r9.push_back(electrons_iter->r9());
    Electron_smin.push_back(electrons_iter->sMin());
    Electron_smaj.push_back(electrons_iter->sMaj());
    Electron_seedid.push_back(electrons_iter->seedId());
    Electron_energymatrix.push_back(electrons_iter->energyMatrix());
    Electron_timingmatrix.push_back(electrons_iter->timingMatrix());
    n_ele++;
    //cout << "finished an iteration of the electrons loop." << endl;
  }
  //cout << "electrons loop finished." << endl; 
  n_pho = 0;
  
  if(phoValid)
  for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
    Photon_pt.push_back(photons_iter->pt());
    Photon_eta.push_back(photons_iter->eta());
    Photon_phi.push_back(photons_iter->phi());
    Photon_m.push_back(photons_iter->m());
    Photon_sigmaietaieta.push_back(photons_iter->sigmaIetaIeta());
    Photon_hoe.push_back(photons_iter->hOverE());
    Photon_ecaliso.push_back(photons_iter->ecalIso());
    Photon_hcaliso.push_back(photons_iter->hcalIso());
    Photon_trkiso.push_back(photons_iter->trkIso());
    Photon_r9.push_back(photons_iter->r9());
    Photon_smin.push_back(photons_iter->sMin());
    Photon_smaj.push_back(photons_iter->sMaj());
    Photon_seedid.push_back(photons_iter->seedId());
    Photon_energymatrix.push_back(photons_iter->energyMatrix());
    Photon_timingmatrix.push_back(photons_iter->timingMatrix());
    
    n_pho++;
  }
  
 if (doL1) {
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    	for( int r = 0; r<280; r++){
	string name ("empty");
	bool algoName_ = false;
	algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
	cout << "getAlgNameFromBit = " << algoName_  << endl;
	cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
	}
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;	
			
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      l1Result_.push_back( l1htbit );
      }
 }

 

 tree->Fill();	
 clearVars();

}

void EGammaOnly_ScoutingNanoAOD::clearVars(){
  genpart_pdg.clear();
  genpart_pt.clear();
  genpart_eta.clear();
  genpart_phi.clear();
  genpart_m.clear();
  genpart_vx.clear();
  genpart_vy.clear();
  genpart_vz.clear();
  genpart_nmoms.clear();
  genpart_mompdg.clear();
  genpart_isPromptFS.clear();
  genpart_isPromptDec.clear();
  genpart_isDirectPromptTauDecayProdFS.clear();
  genpartmomZ_pdg.clear();
  genpartmomZ_pt.clear();
  genpartmomZ_eta.clear();
  genpartmomZ_phi.clear();
  genpartmomZ_m.clear();
  genpartmomZ_vx.clear();
  genpartmomZ_vy.clear();
  genpartmomZ_vz.clear();
  genpartmomZ_isPromptFS.clear();
  genpartmomZ_isPromptDec.clear();
  genpartmomZ_isDirectPromptTauDecayProdFS.clear();
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_d0.clear();
  Electron_dz.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_hoe.clear();
  Electron_ooemoop.clear();
  Electron_missinghits.clear();
  Electron_charge.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_tkiso.clear();
  Electron_r9.clear();
  Electron_smin.clear();
  Electron_smaj.clear();
  Electron_seedid.clear();
  Electron_energymatrix.clear();
  Electron_timingmatrix.clear();
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_m.clear();
  Photon_sigmaietaieta.clear();
  Photon_hoe.clear();
  Photon_ecaliso.clear();
  Photon_hcaliso.clear();
  Photon_trkiso.clear();
  Photon_r9.clear();
  Photon_smin.clear();
  Photon_smaj.clear();
  Photon_seedid.clear();
  Photon_energymatrix.clear();
  Photon_timingmatrix.clear();
}

void EGammaOnly_ScoutingNanoAOD::beginJob() {
  
}

void EGammaOnly_ScoutingNanoAOD::endJob() {
}

void EGammaOnly_ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // HLT paths

  triggerPathsVector.push_back("DST_DoubleMu1_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_Mass10_PFScouting_v*");
  triggerPathsVector.push_back("DST_L1HTT_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_CaloJet40_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT250_CaloScouting_v*");
  triggerPathsVector.push_back("DST_HT410_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT450_PFScouting_v*");

  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }
}

void EGammaOnly_ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(EGammaOnly_ScoutingNanoAOD);
