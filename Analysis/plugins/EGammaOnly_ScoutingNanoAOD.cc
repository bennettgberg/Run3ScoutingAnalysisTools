// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>
#include <boost/phoenix.hpp>

// ROOT includes
#include "TLorentzVector.h"
#include "TPRegexp.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"
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
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

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

  // Flag for simualtion processing
  Bool_t isMC;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> > gensToken;

  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> > pVtxToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> > muonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> > electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> > photonsToken;
  const edm::EDGetTokenT<double> rhoToken;

  // For L1 input
  edm::InputTag algInputTag_;
  edm::InputTag extInputTag_;
  edm::EDGetToken algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  
  // L1 trigger branches
  bool doL1;
  std::vector<std::string> l1Seeds_;
  std::vector<bool> l1Result_;

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
  vector<Bool_t> genpart_fromHardProcessBeforeFSR;
  vector<Bool_t> genpart_fromHardProcessDecayed;
  vector<Bool_t> genpart_fromHardProcessFS;
  vector<Bool_t> genpart_isHardProcess;
  vector<Bool_t> genpart_isLastCopy;
  vector<Bool_t> genpart_isLastCopyBeforeFSR;
  vector<Bool_t> genpart_isPromptFS;
  vector<Bool_t> genpart_isPromptDec;
  vector<Bool_t> genpart_isDirectPromptTauDecayProdFS;

  // Primary Vertex
  const static int max_pv = 1000;
  UInt_t n_pVtx;
  vector<Float16_t> pVtx_x;
  vector<Float16_t> pVtx_y;
  vector<Float16_t> pVtx_z;
  vector<Float16_t> pVtx_xError;
  vector<Float16_t> pVtx_yError;
  vector<Float16_t> pVtx_zError;
  vector<Int_t> pVtx_trksize;
  vector<Float16_t> pVtx_chi2;
  vector<Int_t> pVtx_ndof;
  vector<Bool_t> pVtx_isvalidvtx;
        
  // Muon
  const static int max_mu = 1000;
  UInt_t n_mu;
  vector<Float16_t> Muon_pt;
  vector<Float16_t> Muon_eta;
  vector<Float16_t> Muon_phi;
  vector<Float16_t> Muon_m;
  vector<unsigned int> Muon_type;
  vector<Float16_t> Muon_charge;
  vector<Float16_t> Muon_normchi2;
  vector<Float16_t> Muon_ecaliso;
  vector<Float16_t> Muon_hcaliso;
  vector<Float16_t> Muon_trkiso;
  vector<int> Muon_nvalidpixelhits;
  vector<int> Muon_nvalidstriphits;
  vector<Float16_t> Muon_trkchi2;
  vector<Float16_t> Muon_ndof;
  vector<Float16_t> Muon_dxy;
  vector<Float16_t> Muon_dz;
  vector<Float16_t> Muon_trkqoverp;
  vector<Float16_t> Muon_trklambda;
  vector<Float16_t> Muon_trkpt;
  vector<Float16_t> Muon_trkphi;
  vector<Float16_t> Muon_trketa;
  vector<Float16_t> Muon_dxyError;
  vector<Float16_t> Muon_dzError;
  vector<Float16_t> Muon_trkqoverperror;
  vector<Float16_t> Muon_trklambdaerror;
  vector<Float16_t> Muon_trkphierror;
  vector<Float16_t> Muon_trkdsz;
  vector<Float16_t> Muon_trkdszerror;
  vector<Float16_t> Muon_vx;
  vector<Float16_t> Muon_vy;
  vector<Float16_t> Muon_vz;

  //Electron
  const static int max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> Electron_pt;
  vector<Float16_t> Electron_eta;
  vector<Float16_t> Electron_phi;
  //vector<Float16_t> Electron_m;
  vector<Float16_t> Electron_d0;
  vector<Float16_t> Electron_dz;
  vector<Float16_t> Electron_detain;
  vector<Float16_t> Electron_dphiin;
  vector<Float16_t> Electron_sigmaietaieta;
  vector<Float16_t> Electron_hoe;
  vector<Float16_t> Electron_ooemoop;
  vector<Int_t>	Electron_missinghits;
  vector<Int_t> Electron_charge;
  vector<Float16_t> Electron_ecaliso;
  vector<Float16_t> Electron_hcaliso;
  vector<Float16_t> Electron_tkiso;
  vector<Float16_t> Electron_r9;
  vector<Float16_t> Electron_smin;
  vector<Float16_t> Electron_smaj;
  vector<UInt_t> Electron_seedid;
  vector<vector<Float16_t>> Electron_energymatrix;
  vector<vector<uint32_t>> Electron_detids;
  vector<vector<Float16_t>> Electron_timingmatrix;
  vector<bool> Electron_rechitzerosuppression;

  //trying to get Electron supercluster energy (should be the same as E in H/E)
  vector<Float16_t>   Electron_energy;
  //Photon
  const static int max_pho = 1000;
  UInt_t n_pho;
  vector<Float16_t> Photon_pt;
  vector<Float16_t> Photon_eta;
  vector<Float16_t> Photon_phi;
  vector<Float16_t> Photon_m;
  vector<Float16_t> Photon_sigmaietaieta;
  vector<Float16_t> Photon_hoe;
  vector<Float16_t> Photon_ecaliso;
  vector<Float16_t> Photon_hcaliso;
  vector<Float16_t> Photon_trkiso;
  vector<Float16_t> Photon_r9;
  vector<Float16_t> Photon_smin;
  vector<Float16_t> Photon_smaj;
  vector<Float16_t> Photon_seedid;
  vector<vector<Float16_t>> Photon_energymatrix;
  vector<vector<uint32_t>> Photon_detids;
  vector<vector<Float16_t>> Photon_timingmatrix;
  vector<bool> Photon_rechitzerosuppression;

  // Rho
  UInt_t n_rhoval;
  vector<Float16_t> rho;

  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int event;
  int lumSec;

};

EGammaOnly_ScoutingNanoAOD::EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  isMC(iConfig.existsAs<bool>("isMC")?iConfig.getParameter<bool>("isMC"):true), 
  gensToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gens"))),
  pVtxToken(consumes<std::vector<Run3ScoutingVertex> >(iConfig.getParameter<edm::InputTag>("primaryVtx"))), 
  muonsToken(consumes<std::vector<Run3ScoutingMuon> >(iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken(consumes<std::vector<Run3ScoutingElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken(consumes<std::vector<Run3ScoutingPhoton> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  doL1(iConfig.existsAs<bool>("doL1")?iConfig.getParameter<bool>("doL1"):false) {

  usesResource("TFileService");

  // If doL1, access necessary heads for L1 algorithm
  if(doL1) {
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }
  
  // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  // Event details  
  tree->Branch("lumSec", &lumSec, "lumSec/i" );
  tree->Branch("run", &run, "run/i" );
  tree->Branch("event", &event, "event/i" );
    
  // L1 info
  tree->Branch("l1Result", "std::vector<bool>", &l1Result_, 32000, 0);

  // Gen level particles
  if(isMC) {
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
    tree->Branch("genpart_fromHardProcessBeforeFSR", &genpart_fromHardProcessBeforeFSR);
    tree->Branch("genpart_fromHardProcessDecayed", &genpart_fromHardProcessDecayed);
    tree->Branch("genpart_fromHardProcessFS", &genpart_fromHardProcessFS);
    tree->Branch("genpart_isHardProcess", &genpart_isHardProcess);
    tree->Branch("genpart_isLastCopy", &genpart_isLastCopy);
    tree->Branch("genpart_isLastCopyBeforeFSR", &genpart_isLastCopyBeforeFSR);
    tree->Branch("genpart_isPromptFS", &genpart_isPromptFS);
    tree->Branch("genpart_isPromptDec", &genpart_isPromptDec);
    tree->Branch("genpart_isDirectPromptTauDecayProdFS", &genpart_isDirectPromptTauDecayProdFS);
  }

  // Primary vertex info
  tree->Branch("n_pVtx", &n_pVtx, "n_pVtx/i");
  tree->Branch("pVtx_x", &pVtx_x);
  tree->Branch("pVtx_y", &pVtx_y);
  tree->Branch("pVtx_z", &pVtx_z);
  tree->Branch("pVtx_xError", &pVtx_xError);
  tree->Branch("pVtx_yError", &pVtx_yError);
  tree->Branch("pVtx_zError", &pVtx_zError);
  tree->Branch("pVtx_trksize", &pVtx_trksize);
  tree->Branch("pVtx_chi2", &pVtx_chi2);
  tree->Branch("pVtx_ndof", &pVtx_ndof);
  tree->Branch("pVtx_isvalidvtx", &pVtx_isvalidvtx);
  
  // Muons
  tree->Branch("n_mu", &n_mu, "n_mu/i");
  tree->Branch("Muon_pt", &Muon_pt);
  tree->Branch("Muon_eta", &Muon_eta);
  tree->Branch("Muon_phi", &Muon_phi);
  tree->Branch("Muon_m", &Muon_m);
  tree->Branch("Muon_type", &Muon_type);
  tree->Branch("Muon_charge", &Muon_charge);
  tree->Branch("Muon_normchi2", &Muon_normchi2);
  tree->Branch("Muon_ecaliso", &Muon_ecaliso);
  tree->Branch("Muon_hcaliso", &Muon_hcaliso);
  tree->Branch("Muon_trkiso", &Muon_trkiso);
  tree->Branch("Muon_validpixelhits", &Muon_nvalidpixelhits);
  tree->Branch("Muon_nvalidstriphits", &Muon_nvalidstriphits);
  tree->Branch("Muon_trkchi2", &Muon_trkchi2);
  tree->Branch("Muon_ndof", &Muon_ndof);
  tree->Branch("Muon_dxy", &Muon_dxy);
  tree->Branch("Muon_dz", &Muon_dz);
  tree->Branch("Muon_trkqoverp", &Muon_trkqoverp);
  tree->Branch("Muon_trklambda", &Muon_trklambda);
  tree->Branch("Muon_trkpt", &Muon_trkpt);
  tree->Branch("Muon_trkphi",  &Muon_trkphi);
  tree->Branch("Muon_trketa", &Muon_trketa);
  tree->Branch("Muon_dxyError", &Muon_dxyError);
  tree->Branch("Muon_dzError", &Muon_dzError);
  tree->Branch("Muon_trkqoverperror", &Muon_trkqoverperror);
  tree->Branch("Muon_trklambdaerror", &Muon_trklambdaerror);
  tree->Branch("Muon_trkphierror", &Muon_trkphierror);
  tree->Branch("Muon_trkdsz", &Muon_trkdsz);
  tree->Branch("Muon_trkdszerror", &Muon_trkdszerror);
  tree->Branch("Muon_vx", &Muon_vx);
  tree->Branch("Muon_vy", &Muon_vy);
  tree->Branch("Muon_vz", &Muon_vz);

  // Electrons
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  //tree->Branch("Electron_m", &Electron_m);
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
  tree->Branch("Electron_energymatrix", &Electron_energymatrix);
  tree->Branch("Electron_detids", &Electron_detids);
  tree->Branch("Electron_timingmatrix", &Electron_timingmatrix);
  tree->Branch("Electron_rechitzerosuppression", &Electron_rechitzerosuppression);
  
  // Photons
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
  tree->Branch("Photon_detids", &Photon_detids);
  tree->Branch("Photon_timingmatrix", &Photon_timingmatrix);
  tree->Branch("Photon_rechitzerosuppression", &Photon_rechitzerosuppression);

  // Rho
  tree->Branch("n_rhoval", &n_rhoval, "n_rhoval/i");
  tree->Branch("rho", &rho);
}

EGammaOnly_ScoutingNanoAOD::~EGammaOnly_ScoutingNanoAOD() {
}

void EGammaOnly_ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
<<<<<<< HEAD
  using namespace fastjet;
  using namespace fastjet::contrib;

  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);
  bool beamValid = beamSpotH.isValid();
=======
      
  Handle<vector<reco::GenParticle> >gensH;
  if(isMC) {
    iEvent.getByToken(gensToken, gensH);
  }

  Handle<vector<Run3ScoutingVertex> > pVtxH;
  iEvent.getByToken(pVtxToken, pVtxH);
  bool pVtxValid = pVtxH.isValid();

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  bool muoValid = muonsH.isValid();
>>>>>>> 088daaf93fd01d3290d41beb3b086d8f45127aaa

  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  bool eleValid = electronsH.isValid();

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  bool phoValid = photonsH.isValid();

  Handle<double > rhoH;
  iEvent.getByToken(rhoToken, rhoH);
  bool rhoValid = rhoH.isValid();

  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  ///////////////////// for genmatching /////////////////////////////////////////////

  n_gen=0;

  if(isMC) {
    if(gensH.isValid()) {
      
      for (auto gen_iter = gensH->begin(); gen_iter != gensH->end(); ++gen_iter) {
	
	int nmoms = 0;
	if((std::abs(gen_iter->pdgId())==11 || std::abs(gen_iter->pdgId())==13 || std::abs(gen_iter->pdgId())==15)) {
	  
	  genpart_pdg.push_back(gen_iter->pdgId());
	  genpart_pt.push_back(gen_iter->pt());
	  genpart_eta.push_back(gen_iter->eta());
	  genpart_phi.push_back(gen_iter->phi());
	  genpart_m.push_back(gen_iter->mass());
	  genpart_vx.push_back(gen_iter->vx());
	  genpart_vy.push_back(gen_iter->vy());
	  genpart_vz.push_back(gen_iter->vz());
	  
	  nmoms = gen_iter->numberOfMothers();
	  genpart_nmoms.push_back(nmoms);
	  if(gen_iter->numberOfMothers() > 0) genpart_mompdg.push_back(gen_iter->mother(0)->pdgId());
	  
	  genpart_fromHardProcessBeforeFSR.push_back(gen_iter->fromHardProcessBeforeFSR());
	  genpart_fromHardProcessDecayed.push_back(gen_iter->fromHardProcessDecayed());
	  genpart_fromHardProcessFS.push_back(gen_iter->fromHardProcessFinalState());
	  genpart_isHardProcess.push_back(gen_iter->isHardProcess());
	  genpart_isLastCopy.push_back(gen_iter->isLastCopy());
	  genpart_isLastCopyBeforeFSR.push_back(gen_iter->isLastCopyBeforeFSR());
	  genpart_isPromptFS.push_back(gen_iter->isPromptFinalState());
	  genpart_isPromptDec.push_back(gen_iter->isPromptDecayed());
	  genpart_isDirectPromptTauDecayProdFS.push_back(gen_iter->isDirectPromptTauDecayProductFinalState());
	  
	  n_gen++;
	} //end if gen lepton
	
      } //end for genpart loop
      
    } // end of gensH.isValid()
  } // end of isMC

  // Primary Vertex
  n_pVtx = 0;
  if(pVtxValid) {
    for (auto &pVtx : *pVtxH) {
      auto *pVtx_iter = &pVtx;
      pVtx_x.push_back(pVtx_iter->x());
      pVtx_y.push_back(pVtx_iter->y());
      pVtx_z.push_back(pVtx_iter->z());
      pVtx_xError.push_back(pVtx_iter->xError());
      pVtx_yError.push_back(pVtx_iter->yError());
      pVtx_zError.push_back(pVtx_iter->zError());
      pVtx_trksize.push_back(pVtx_iter->tracksSize());
      pVtx_chi2.push_back(pVtx_iter->chi2());
      pVtx_ndof.push_back(pVtx_iter->ndof());
      pVtx_isvalidvtx.push_back(pVtx_iter->isValidVtx());
      n_pVtx++;
    } // end primary vertex loop
  } // end pVtxValid condition

  // Muons
  n_mu = 0;
  if(muoValid) {
    for (auto &muo : *muonsH) {
      auto *muons_iter = &muo;
      Muon_pt.push_back(muons_iter->pt());
      Muon_eta.push_back(muons_iter->eta());
      Muon_phi.push_back(muons_iter->phi());
      Muon_m.push_back(muons_iter->m());
      Muon_type.push_back(muons_iter->type());
      Muon_charge.push_back(muons_iter->charge());
      Muon_normchi2.push_back(muons_iter->normalizedChi2());
      Muon_ecaliso.push_back(muons_iter->ecalIso());
      Muon_hcaliso.push_back(muons_iter->hcalIso());
      Muon_trkiso.push_back(muons_iter->trackIso());
      Muon_nvalidpixelhits.push_back(muons_iter->nValidPixelHits());
      Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
      Muon_trkchi2.push_back(muons_iter->trk_chi2());
      Muon_ndof.push_back(muons_iter->trk_ndof());
      Muon_dxy.push_back(muons_iter->trk_dxy());
      Muon_dz.push_back(muons_iter->trk_dz());
      Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
      Muon_trklambda.push_back(muons_iter->trk_lambda());
      Muon_trkpt.push_back(muons_iter->trk_pt());
      Muon_trkphi.push_back(muons_iter->trk_phi());
      Muon_trketa.push_back(muons_iter->trk_eta());
      Muon_dxyError.push_back(muons_iter->trk_dxyError());
      Muon_dzError.push_back(muons_iter->trk_dzError());
      Muon_trkqoverperror.push_back(muons_iter->trk_qoverpError());
      Muon_trklambdaerror.push_back(muons_iter->trk_lambdaError());
      Muon_trkphierror.push_back(muons_iter->trk_phiError());
      Muon_trkdsz.push_back(muons_iter->trk_dsz());
      Muon_trkdszerror.push_back(muons_iter->trk_dszError());
      Muon_vx.push_back(muons_iter->trk_vx());
      Muon_vy.push_back(muons_iter->trk_vy());
      Muon_vz.push_back(muons_iter->trk_vz());

      n_mu++;
    } // end muon loop
  } // end muoValid condition

  // Electrons
  n_ele = 0;
  if(eleValid) {
    for (auto &ele : *electronsH) {
      auto *electrons_iter = &ele;
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      //Electron_m.push_back(electrons_iter->m());
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
      Electron_detids.push_back(electrons_iter->detIds());
      Electron_timingmatrix.push_back(electrons_iter->timingMatrix());
      //Electron_rechitzerosuppression.push_back(electrons_iter->rechitZeroSuppression());
      n_ele++;

    } // end electron loop
  } // end eleValid condition

  // Muons
  n_pho = 0;  
  if(phoValid) {

    for (auto &pho : *photonsH) {
      auto *photons_iter = &pho;
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
      Photon_detids.push_back(photons_iter->detIds());
      Photon_timingmatrix.push_back(photons_iter->timingMatrix());
      //Photon_rechitzerosuppression.push_back(photons_iter->rechitZeroSuppression());
      n_pho++;

    } // end photon loop
  } // end phoValid condition

  // Rho
  n_rhoval = 0;
  if(rhoValid) {
    rho.push_back(*rhoH);
    n_rhoval++;
  }

  // Fill L1 seeds
  if(doL1) {
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    /*
    for( int r = 0; r<512; r++){
      string name("empty");
      bool algoName_ = false;
      algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
      l1GtUtils_->getAlgNameFromBit(r,name);
      cout << "getAlgNameFromBit = " << algoName_  << endl;
      cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
    }
    */
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
<<<<<<< HEAD
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
  genpart_fromHardProcessBeforeFSR.clear();
  genpart_isPromptFS.clear();
  genpart_isPromptDec.clear();
  genpart_isDirectPromptTauDecayProdFS.clear();
=======
  if(isMC) {
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
    genpart_fromHardProcessBeforeFSR.clear();
    genpart_fromHardProcessDecayed.clear();
    genpart_fromHardProcessFS.clear();
    genpart_isHardProcess.clear();
    genpart_isLastCopy.clear();
    genpart_isLastCopyBeforeFSR.clear();
    genpart_isPromptFS.clear();
    genpart_isPromptDec.clear();
    genpart_isDirectPromptTauDecayProdFS.clear();
  }
  pVtx_x.clear();
  pVtx_y.clear();
  pVtx_z.clear();
  pVtx_xError.clear();
  pVtx_yError.clear();
  pVtx_zError.clear();
  pVtx_trksize.clear();
  pVtx_chi2.clear();
  pVtx_ndof.clear();
  pVtx_isvalidvtx.clear();
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_type.clear();
  Muon_charge.clear();
  Muon_normchi2.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkchi2.clear();
  Muon_ndof.clear();
  Muon_dxy.clear();
  Muon_dz.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_dxyError.clear();
  Muon_dzError.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkphierror.clear();
  Muon_trkdsz.clear();
  Muon_trkdszerror.clear();
  Muon_vx.clear();
  Muon_vy.clear();
  Muon_vz.clear();
>>>>>>> 088daaf93fd01d3290d41beb3b086d8f45127aaa
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  //Electron_m.clear();
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
  Electron_detids.clear();
  Electron_timingmatrix.clear();
  Electron_rechitzerosuppression.clear();
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
  Photon_detids.clear();
  Photon_timingmatrix.clear();
  Photon_rechitzerosuppression.clear();
  rho.clear();
}

void EGammaOnly_ScoutingNanoAOD::beginJob() {
  
}

void EGammaOnly_ScoutingNanoAOD::endJob() {
}

void EGammaOnly_ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
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
