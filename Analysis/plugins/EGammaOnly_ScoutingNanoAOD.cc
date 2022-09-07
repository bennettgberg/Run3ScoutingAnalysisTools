// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include "TLorentzVector.h"
#include "TPRegexp.h"
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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

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

  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> > electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> > photonsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> > gensToken;

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
  vector<Bool_t> genpart_fromHardProcessBeforeFSR;
  vector<Bool_t> genpart_fromHardProcessDecayed;
  vector<Bool_t> genpart_fromHardProcessFS;
  vector<Bool_t> genpart_isHardProcess;
  vector<Bool_t> genpart_isLastCopy;
  vector<Bool_t> genpart_isLastCopyBeforeFSR;
  vector<Bool_t> genpart_isPromptFS;
  vector<Bool_t> genpart_isPromptDec;
  vector<Bool_t> genpart_isDirectPromptTauDecayProdFS;

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

  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;

};

EGammaOnly_ScoutingNanoAOD::EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  electronsToken(consumes<std::vector<Run3ScoutingElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken(consumes<std::vector<Run3ScoutingPhoton> >(iConfig.getParameter<edm::InputTag>("photons"))), 
  gensToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gens"))) {

  usesResource("TFileService");

  // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree"       , "tree");

  // Event details  
  tree->Branch("lumSec"		, &lumSec			 , "lumSec/i" );
  tree->Branch("run"			, &run				 , "run/i" );
    
  // Beamspot info
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
  tree->Branch("genpart_fromHardProcessBeforeFSR", &genpart_fromHardProcessBeforeFSR);
  tree->Branch("genpart_fromHardProcessDecayed", &genpart_fromHardProcessDecayed);
  tree->Branch("genpart_fromHardProcessFS", &genpart_fromHardProcessFS);
  tree->Branch("genpart_isHardProcess", &genpart_isHardProcess);
  tree->Branch("genpart_isLastCopy", &genpart_isLastCopy);
  tree->Branch("genpart_isLastCopyBeforeFSR", &genpart_isLastCopyBeforeFSR);
  tree->Branch("genpart_isPromptFS", &genpart_isPromptFS);
  tree->Branch("genpart_isPromptDec", &genpart_isPromptDec);
  tree->Branch("genpart_isDirectPromptTauDecayProdFS", &genpart_isDirectPromptTauDecayProdFS);
  
  //Electrons
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
  tree->Branch("Photon_detids", &Photon_detids);
  tree->Branch("Photon_timingmatrix", &Photon_timingmatrix);
  tree->Branch("Photon_rechitzerosuppression", &Photon_rechitzerosuppression);
}


EGammaOnly_ScoutingNanoAOD::~EGammaOnly_ScoutingNanoAOD() {
}

void EGammaOnly_ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;

  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);
  bool beamValid = beamSpotH.isValid();

  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  bool eleValid = electronsH.isValid();

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  bool phoValid = photonsH.isValid();

  Handle<vector<reco::GenParticle> >gensH;
  iEvent.getByToken(gensToken, gensH);

  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  // beamspot
  if(beamValid) {
    beamspot_x = beamSpotH->x0();
    beamspot_y = beamSpotH->y0();
    beamspot_z = beamSpotH->z0();
  }
  else {
    beamspot_x = 1e20;
    beamspot_y = 1e20;
    beamspot_z = 1e20;
  }
  
  ///////////////////// for genmatching /////////////////////////////////////////////

  n_gen=0;

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
  genpart_fromHardProcessBeforeFSR.clear();
  genpart_isPromptFS.clear();
  genpart_isPromptDec.clear();
  genpart_isDirectPromptTauDecayProdFS.clear();
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
