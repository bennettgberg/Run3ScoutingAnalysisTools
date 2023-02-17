#include "data_robustanalyzer.hh"
#include <iostream>
#include <sstream>
#include <numeric>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

const float h0 = 16.; // 16.; //constant parameter for the rho-dep cut, and the H<.2E+h0 cut
const float h1 = 0.2; //factor on E for the rho-dependent cut
const float Crho = 0.07; //factor on rho for the rho-dependent cut
//set issig to false to skip genMatching (only bkg will have the data).
const bool issig = true; // false; //true; // false;
int ngm = 0;
int npassed = 0;
int noldpassed = 0;
//true for ElectronGun, false for real simulations like DY
bool isGun = true;
// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool isDoubleElectron){

  isDiEl = isDoubleElectron;
  
  TFile *inpfile = TFile::Open(filename,"READ");
  cout<<"Initializing for file: "<<filename<<endl;

  tree = new TTreeReader("mmtree/tree",inpfile);
  //bsx = new TTreeReaderArray<float>((*tree), "beamspot_x");
  //bsy = new TTreeReaderArray<float>((*tree), "beamspot_y");
  //bsz = new TTreeReaderArray<float>((*tree), "beamspot_z");
  n_ele = new TTreeReaderValue<UInt_t>((*tree), "n_ele");
  n_mu = new TTreeReaderValue<UInt_t>((*tree), "n_mu");
  n_rho = new TTreeReaderValue<UInt_t>((*tree), "n_rhoval");
  ele_pt = new TTreeReaderValue<vector<float>>((*tree), "Electron_pt");
  ele_eta = new TTreeReaderValue<vector<float>>((*tree), "Electron_eta");
  ele_phi = new TTreeReaderValue<vector<float>>((*tree), "Electron_phi");
//  ele_m = new TTreeReaderValue<vector<float>>((*tree), "Electron_m");
  ele_d0 = new TTreeReaderValue<vector<float>>((*tree), "Electron_d0");
  ele_dz = new TTreeReaderValue<vector<float>>((*tree), "Electron_dz");
  ele_detain = new TTreeReaderValue<vector<float>>((*tree), "Electron_detain");
  ele_dphiin = new TTreeReaderValue<vector<float>>((*tree), "Electron_dphiin");
  ele_sigmaietaieta = new TTreeReaderValue<vector<float>>((*tree), "Electron_sigmaietaieta");
  ele_hoe = new TTreeReaderValue<vector<float>>((*tree), "Electron_hoe");
  ele_ooemoop = new TTreeReaderValue<vector<float>>((*tree), "Electron_ooemoop");
  ele_mhits = new TTreeReaderValue<vector<int>>((*tree), "Electron_missinghits");
  ele_charge = new TTreeReaderValue<vector<int>>((*tree), "Electron_charge");
  ele_ecaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_ecaliso");
  ele_hcaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_hcaliso");
  ele_tkiso = new TTreeReaderValue<vector<float>>((*tree), "Electron_tkiso");
  ele_r9 = new TTreeReaderValue<vector<float>>((*tree), "Electron_r9");
  ele_smin = new TTreeReaderValue<vector<float>>((*tree), "Electron_smin");
  ele_smaj = new TTreeReaderValue<vector<float>>((*tree), "Electron_smaj");
  rho      = new TTreeReaderValue<vector<float>>((*tree), "rho");
  ele_seedid = new TTreeReaderValue<vector<unsigned int>>((*tree), "Electron_seedid");
  ele_enemat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_energymatrix");
  ele_timmat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_timingmatrix");

  //Muon info is needed for computing a fairer efficiency.
  mu_pt = new TTreeReaderValue<vector<float>>((*tree), "Muon_pt");  

  //now for the gen branches
  genpart_pt = new TTreeReaderValue<vector<float>>((*tree), "genpart_pt");
  genpart_eta = new TTreeReaderValue<vector<float>>((*tree), "genpart_eta");
  genpart_phi = new TTreeReaderValue<vector<float>>((*tree), "genpart_phi");
  genpart_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_pdg");

  //need this for real (DY) samples, but NOT for gun samples
  if(!isGun) {
      genpart_isFinalState = new TTreeReaderValue<vector<bool>>((*tree), "genpart_fromHardProcessFS");
  }
  else {
      genpart_isFinalState = new TTreeReaderValue<vector<bool>>((*tree), "genpart_isLastCopy"); 
  }

  n_genpart = new TTreeReaderValue<UInt_t>((*tree), "n_genpart");
  
  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
data_robustanalyzer::~data_robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void data_robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = tree->GetEntries();
  //cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = 0; //splitCnt*(totEntries/nCores);
  int endevent = totEntries; //(splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  //endevent = endevent<totEntries?endevent:totEntries;
  //tree->SetEntriesRange(beginevent, endevent);
  int event = beginevent-1;

  // Count events passing certain selections
  int nosel=0, noselZwind=0;

  // Define the histograms
  addhist("nosel");
  addhist("noselZwind");

  // vector of electron indices
  vector<int> noselelidx;
  vector<int> noselZwindelidx;
  
  //make 2D plot to compare H/E agains pt of electrons in the Z window.
  TH2F* hoevpt = new TH2F("hoevspt", "H/E vs. p_T for genmatched e^-s", 1000, 0, 10, 4000, 0, 4.0);
  TH2F* hoevptBkg = new TH2F("hoevsptBkg", "H/E vs. p_T for e^-s NOT genmatched", 1000, 0, 10, 4000, 0, 4.0);
  TH2F* EvsE = new TH2F("EvsE", "Sum of ECAL energy matrix vs. Reconstructed E", 15000, 0, 1500, 15000, 0, 1500);
    //make this to find reasonable starting cut for rho ?
  TH2F* RhovsH = new TH2F("RhovsH", "Hcal energy vs. Pileup #rho", 600, 0., 60., 10000, 0., 100.); 
    //make this to find reasonable values for the other parameters in the equation
  TH2F* HmRhoVsE = new TH2F("HmRhoVsE", "H/E - #rho/E Vs. E", 500, 0., 500., 300, -1.5, 1.5); 
  //hists for efficiency as function of pT
    //all gen electrons
  TH1F* all_gen = new TH1F("allgen", "", 500, 0., 50.);
    //electrons in events with 2 muons of pT>1.9 GeV (hence should pass L1 trigger)
  TH1F* all_l1 = new TH1F("alll1", "", 500, 0., 50.);
    //electrons in events with 2 muons of pT>1.9 GeV (hence should pass L1 trigger), vs. Rho instead of pT
  TH1F* all_l1_vRho = new TH1F("alll1vRho", "", 100, 0., 100.);
    //gen electrons in the middle barrel region
  TH1F* mid_gen = new TH1F("midgen", "", 500, 0., 50.);
    //gen electrons in the outer barrel region
  TH1F* out_gen = new TH1F("outgen", "", 500, 0., 50.);
    //gen electrons in the endcaps
  TH1F* ecs_gen = new TH1F("ecsgen", "", 500, 0., 50.);
    //all scouting electrons
  TH1F* all_sct = new TH1F("allsct", "", 500, 0., 50.);
    //all scouting electrons that pass the new cut
  TH1F* pass_sct = new TH1F("passsct", "", 500, 0., 50.);
    //all scouting electrons that pass the old cut
  TH1F* oldpass_sct = new TH1F("oldpasssct", "", 500, 0., 50.);
    //lead scouting electron only (all)
  TH1F* alllead_sct = new TH1F("allleadsct", "", 500, 0., 50.);
    //lead scouting electron only (passing new cut)
  TH1F* passlead_sct = new TH1F("passleadsct", "", 500, 0., 50.);
    //lead scouting electron only (passing old cut)
  TH1F* oldpasslead_sct = new TH1F("oldpassleadsct", "", 500, 0., 50.);
    //resolution=gen pt - scouting pt
  TH1F* pt_res_sct = new TH1F("ptressct", "", 1000, -50., 50.);
    //ptratio = (scouting e pt) / (gen e pt)
  //TH1F* pt_ratio_sctgen = new TH1F("ptratiosctgen", "", 200, 0., 2.);
  TH1F* e_ratio_sctgen = new TH1F("eratiosctgen", "", 200, 0., 2.);
    //energy ratio, except only for electrons with pt < 20
  TH1F* e_ratio_sctgenl20 = new TH1F("eratiosctgenl20", "", 200, 0., 2.);
    //scouting electrons in the middle region (|eta|<1)
  TH1F* mid_sct = new TH1F("midsct", "", 500, 0., 50.);
    //scouting electrons in the outer barrel (1<|eta|<1.479)
  TH1F* out_sct = new TH1F("outsct", "", 500, 0., 50.);
    //scouting electrons in the endcaps (|eta|>1.479)
  TH1F* ecs_sct = new TH1F("ecssct", "", 500, 0., 50.);
    //all gen-matched electrons
  TH1F* all_sig = new TH1F("allsig", "", 500, 0., 50.);
    //all gen-matched electrons, vs. Rho instead of pT
  TH1F* all_sig_vRho = new TH1F("allsigvRho", "", 100, 0., 100.);
    //gen-matched electrons passing the new H/E cut
  TH1F* passed_sig = new TH1F("passedsig", "", 500, 0., 50.);
    //gen-matched electrons passing the rho-dependent H/E cut
  TH1F* rho_passed_sig = new TH1F("rhopassedsig", "", 500, 0., 50.);
    //gen-matched electrons passing the old H/E cut H/E<.2
  TH1F* old_passed_sig = new TH1F("oldpassedsig", "", 500, 0., 50.);
    //gen-matched electrons passing the old H/E cut H/E<.2, vs. Rho instead of pT
  TH1F* old_passed_sig_vRho = new TH1F("oldpassedsigvRho", "", 100, 0., 100.);
    //all non-gen-matched electrons
  TH1F* all_bkg = new TH1F("allbkg", "", 500, 0., 50.);
    //non-gen-matched electrons passing the new H/E cut
  TH1F* passed_bkg = new TH1F("passedbkg", "", 500, 0., 50.);
    //non-gen-matched electrons passing the old H/E cut H/E<.2
  TH1F* old_passed_bkg = new TH1F("oldpassedbkg", "", 500, 0., 50.);
    //number of gen electrons in the event
  TH1D* n_gen_ele = new TH1D("ngenele", "", 15, 0, 15);
  TH1D* n_reco_ele = new TH1D("nrecoele", "", 15, 0, 15);

  //now hists as function of eta instead of pt
    //gen electrons
  TH1F* eta_gen = new TH1F("etagen", "", 2000, -10.0, 10.);
    //all scouting electrons
  TH1F* eta_sct = new TH1F("etasct", "", 2000, -10.0, 10.);
    //all genmatched scouting electrons
  TH1F* eta_sig = new TH1F("etasig", "", 2000, -10.0, 10.);
    //all non-genmatched scouting electrons
  TH1F* eta_bkg = new TH1F("etabkg", "", 2000, -10.0, 10.);
    //gen-matched scouting electrons passing old H/E<0.2 cut
  TH1F* eta_old_passed_sig = new TH1F("etaoldpassedsig", "", 2000, -10.0, 10.);
   //gen-matched scouting electrons passing new H<0.2E+p_1 cut
  TH1F* eta_new_passed_sig = new TH1F("etanewpassedsig", "", 2000, -10.0, 10.);
  //vector of TH2 plots comparing H and E
  // plot 0 is all pt, plot 1 is 2<pt<5, plot2 is 5<pt<8, plot3 is 8<pt<10 GeV.
  //hVe0: all
  //hVe1: pt<5; hVe2: 5<pt<8; hVe3: 8<pt<10; hVe4: 10<pt<15; hVe5:15<pt<20; hVe6:pt>20 GeV
   //GENMATCHED ONLY!
  vector<TH2F*> hVe(8);
  vector<TH1D*> deta_all(8);
  vector<TH1D*> deta_gmd(8);
  vector<TH1D*> qdphi_all(8);
  vector<TH1D*> qdphi_gmd(8);
  for(uint32_t i = 0; i < hVe.size(); i++) {
      stringstream hname;
      hname << "hVe" << i;
      //hVe[i] = new TH2F(hname.str().c_str(), "H vs. E", 2000, 0, 200, 15000, 0, 1500.0);
      hVe[i] = new TH2F(hname.str().c_str(), "H vs. E", 15000, 0, 1500, 2000, 0, 200);
      hVe[i]->GetYaxis()->SetTitle("HCal energy deposits (GeV)");
      hVe[i]->GetXaxis()->SetTitle("ECal energy deposits (GeV)");
      hVe[i]->SetMarkerColor(2);
        
      //plot of deltaEta for all possible scoutingElectron-genElectron pairs
      stringstream detaname; detaname << "detaall" << i;
      deta_all[i] = new TH1D(detaname.str().c_str(), "", 300, 0, 3.0 );
      //plot of deltaEta for only genMatched scoutingElectrons-genElectron pairs
      stringstream detagname; detagname << "detagmd" << i;
      deta_gmd[i] = new TH1D(detagname.str().c_str(), "", 300, 0, 3.0);
      //plot of deltaPhi for all possible genMatched scoutingElectrons-genElectron pairs
      stringstream qdphiname; qdphiname << "qdphiall" << i;
      qdphi_all[i] = new TH1D(qdphiname.str().c_str(), "", 630, -3.15, 3.15);
      //plot of deltaPhi for only genMatched scoutingElectrons-genElectron pairs
      stringstream qdphigname; qdphigname << "qdphigmd" << i;
      qdphi_gmd[i] = new TH1D(qdphigname.str().c_str(), "", 630, -3.15, 3.15);
  } //end i loop instantiating hVe.

    int nLastCut = 0; //keep track of how many gen electrons are cut by the 'isLastCopy' requirement
  hoevpt->GetXaxis()->SetTitle("#p_T (GeV)");
  hoevpt->GetYaxis()->SetTitle("H/E");
  hoevpt->SetMarkerColor(2);
  hoevptBkg->SetMarkerColor(3);
  RhovsH->GetYaxis()->SetTitle("Scouting electron H (GeV)");
  RhovsH->GetXaxis()->SetTitle("Pileup #rho");
  HmRhoVsE->GetXaxis()->SetTitle("Scouting electron E (GeV)");
  HmRhoVsE->GetYaxis()->SetTitle("Scouting electron (H-C_{#rho}*#rho)/E");
  int nZeroe = 0;
  int nNonzeroe = 0;
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    //if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if((*(*n_ele))<0) cout<<"Error!!! Wrong technical event processing. Negative number of electrons in event."<<endl;;
    if((*(*n_ele))==0) { 
        nZeroe++;
        //continue;
    }
    nNonzeroe++;
    
    //count how many muons there are in this event with pT>1.9 GeV (event should pass L1 trig if goodmucount>=2)
    int goodmucount = 0;
    for(unsigned int mucount=0; mucount < *(*n_mu); mucount++) {
        if((*mu_pt)->at(mucount) > 1.9){
            goodmucount++;
        }
    }

    int ngenele = 0;
    int ngenNotele = 0;
    //fill in the gen electrons histogram
    //for(unsigned int genctr=0; genctr < (*genpart_pt)->size(); genctr++) {
    for(unsigned int genctr=0; genctr < *(*n_genpart); genctr++) {
        //fill only electrons! AND only if is in final state! AND within the ECAL!
        float gen_eta = (*genpart_eta)->at(genctr);
        float gen_pt = (*genpart_pt)->at(genctr);
        //ALSO going to exclude the EB/EE gap region, where reconstruction efficiency is low.
        if(abs((*genpart_pdg)->at(genctr)) == 11 && (*genpart_isFinalState)->at(genctr) && fabs(gen_eta) < 2.5 && !(fabs(gen_eta) > 1.44 && fabs(gen_eta)<1.56)) {
            all_gen->Fill(gen_pt);
            eta_gen->Fill(gen_eta);
            ngenele++;
        }
        if(abs((*genpart_pdg)->at(genctr)) == 11 && (*genpart_isFinalState)->at(genctr) && fabs(gen_eta) < 2.5 && !(fabs(gen_eta) > 1.44 && fabs(gen_eta)<1.56) && goodmucount > 1) {
            all_l1->Fill(gen_pt);          
            if(gen_pt > 5.0) {
                float my_rho = (*rho)->at(0);
                all_l1_vRho->Fill(my_rho);
            }
        }
        else {
            if(abs((*genpart_pdg)->at(genctr)) == 11) {
                nLastCut++;
            }
            ngenNotele++;
        }
    } //end filling in the gen electrons hist
    //cout << "ngenele: " << ngenele << ", ngenNotele: " << ngenNotele << endl;
    n_gen_ele->Fill(ngenele);    
    n_reco_ele->Fill((*(*n_ele)));

    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine

    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];

      float this_eta = (*ele_eta)->at(elidx);

      //barrel ONLY at this time.
      //if(abs(this_eta) > 1.479) continue;
      //endcap only at THIS time.
      //if(abs(this_eta) < 1.479) continue;

      noselelidx.push_back(elidx);

      
    }// End of loop on electrons in the event loop

    // Event level selection on the electrons
    //pair<int,int> noselZwindels = inZwindow(noselelidx);
    //cout << "starting genMatch." << endl;
    //pair<int,int> noselZwindels = genMatch(noselelidx, issig);
    vector<int> noselZwindels = genMatch(noselelidx, issig, deta_all, qdphi_all, deta_gmd, qdphi_gmd, pt_res_sct, e_ratio_sctgen, e_ratio_sctgenl20);
    for(int matched_el : noselZwindels) {
        noselZwindelidx.push_back(matched_el);
    }
    //if(noselZwindels.first!=-1) {
    //  noselZwindelidx.push_back(noselZwindels.first);
    //}
    //else {
    //  noselZwindelidx.push_back(-1);
    //}
    //if(noselZwindels.second != -1) {
    //  noselZwindelidx.push_back(noselZwindels.second);
    //}
    //else {
    //  noselZwindelidx.push_back(-1);
    //}
    //cout << "done with genMatch." << endl;
    //if(noselZwindels.first == -1 || noselZwindels.second == -1) {
    //    cout << "*****WARNING: No genMatch found for event with the following electrons:*****" << endl;
    //    if(noselZwindels.first == -1) cout << "(first not found.)" << endl;
    //    if(noselZwindels.second == -1) cout << "(second not found.)" << endl;
    //    cout << "pt, eta, phi:" << endl;
    //    for( auto el: noselelidx) {
    //        cout << (*ele_pt)->at(el) << ", " << (*ele_eta)->at(el) << ", " << (*ele_phi)->at(el) << endl;
    //    }
    //    cout << "And gen particles: " << endl;
    //    for(int i=0; i < 2; i++) {
    //        cout << "pt, eta, phi, pdg:" << endl;
    //        cout << i << ".) " << (*genpart_pt)->at(i) << ", " << (*genpart_eta)->at(i) << ", " << (*genpart_phi)->at(i) 
    //             << ", " << (*genpart_phi)->at(i) << ", " << (*genpart_pdg)->at(i) << endl;
    //    }
    //    cout << "******************" << endl;
    //}
    
    //now add all the other electrons (not in the Z window) to the notZwindow array.
    vector<int> notZwindelidx;
    for(auto elidx : noselelidx) {
        bool isMatched = false;
        //if( elidx != noselZwindels.first && elidx != noselZwindels.second ){
        for(auto matchedidx : noselZwindels) {
            if(elidx == matchedidx) {
                isMatched = true;
                break;
            }
        }
        if(!isMatched) {
            notZwindelidx.push_back(elidx);
        }
    } //end loop over all electrons

    if(noselelidx.size()>0) nosel++;
    //fillhistinevent("nosel", noselelidx, hVe);
    //cout << "starting fillhist." << endl;
    fillhistinevent("nosel", noselelidx, hVe, noselZwindelidx, notZwindelidx, all_sig, all_sig_vRho, all_sct, mid_sct, out_sct, ecs_sct, passed_sig, rho_passed_sig, all_bkg, passed_bkg, old_passed_sig, old_passed_sig_vRho, old_passed_bkg, eta_sig, eta_old_passed_sig, eta_new_passed_sig, eta_sct, eta_bkg, mid_gen, out_gen, ecs_gen, pass_sct, oldpass_sct, alllead_sct, passlead_sct, oldpasslead_sct, EvsE, RhovsH, HmRhoVsE, goodmucount);
    //cout << "done with fillhistinevent." << endl;
    //if(noselZwindelidx.size()>0) noselZwind++;
    if(noselZwindelidx.size()>1 && (noselZwindelidx[0] > -1 || noselZwindelidx[1] > -1) ) noselZwind++;
    //fillhistinevent("noselZwind", noselZwindelidx, vector<TH2F*>(0), vector<int>(0));

    //bpgadding
    fillHoEvsPt(hoevpt, noselZwindelidx);
    fillHoEvsPt(hoevptBkg, notZwindelidx);

    //now find the efficiencies, plot them!
    ////jk DON't (need to do that after adding all files together)
    //passed_sig->Divide(all_sig);
    //passed_bkg->Divide(all_bkg);

    // Clear all vector
    noselelidx.clear();
    noselZwindelidx.clear();
    notZwindelidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<endl;
  cout << "Events with zero electrons: " << nZeroe << "; events with nonzero electrons: " << nNonzeroe << endl;
  //cout << "n gen electrons cut by isLast requirement: " << nLastCut << endl;
  cout << "Number of gen matched electrons: " << ngm << "; Number of passing electrons: " << npassed << "; Noldpassing: " << noldpassed << endl;
}

// Function to fill a set of histograms for gen particles
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> elidx, vector<TH2F*> hVe, vector<int> elidx_gm, vector<int> elidx_notgm, TH1F* all_sig, TH1F* all_sig_vRho, TH1F* all_sct, TH1F* mid_sct, TH1F* out_sct, TH1F* ecs_sct, TH1F* passed_sig, TH1F* rho_passed_sig, TH1F* all_bkg, TH1F* passed_bkg, TH1F* old_passed_sig, TH1F* old_passed_sig_vRho, TH1F* old_passed_bkg, TH1F* eta_sig, TH1F* eta_old_passed_sig, TH1F* eta_new_passed_sig, TH1F* eta_sct, TH1F* eta_bkg, TH1F* mid_gen, TH1F* out_gen, TH1F* ecs_gen, TH1F* pass_sct, TH1F* oldpass_sct, TH1F* alllead_sct, TH1F* passlead_sct, TH1F* oldpasslead_sct, TH2F* EvsE, TH2F* RhovsH, TH2F* HmRhoVsE, int goodmucount) {

  if(elidx.size()==0) return;

  int npassedevt = 0;
  int noldpassedevt = 0;
  int ngmevt = 0;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  //new plots for comparing el multiplicities
  TH1F* elmult_gm = (TH1F*) outfile->Get(selection+"sct_elmult_gm");
  //for these include ALL scouting electrons, not just the ones that are genmatched
  TH1F* elmult_passed = (TH1F*) outfile->Get(selection+"sct_elmult_passed");
  TH1F* elmult_oldpassed = (TH1F*) outfile->Get(selection+"sct_elmult_oldpassed");
  //
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sctbar_elpt");
  //TH1F* barelm = (TH1F*) outfile->Get(selection+"sctbar_elm");
  TH1F* bareld0 = (TH1F*) outfile->Get(selection+"sctbar_eld0");
  TH1F* barellog10d0 = (TH1F*) outfile->Get(selection+"sctbar_ellog10d0");
  TH1F* bareldz = (TH1F*) outfile->Get(selection+"sctbar_eldz");
  TH1F* bareldetain = (TH1F*) outfile->Get(selection+"sctbar_eldetain");
  TH1F* bareldphiin = (TH1F*) outfile->Get(selection+"sctbar_eldphiin");
  TH1F* barelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctbar_elsigmaietaieta");
  TH1F* barelhoe = (TH1F*) outfile->Get(selection+"sctbar_elhoe");
  TH1F* barelooemoop = (TH1F*) outfile->Get(selection+"sctbar_elooemoop");
  TH1F* barelmhits = (TH1F*) outfile->Get(selection+"sctbar_elmhits");
  TH1F* barelcharge = (TH1F*) outfile->Get(selection+"sctbar_elcharge");
  TH1F* barelecaliso = (TH1F*) outfile->Get(selection+"sctbar_elecaliso");
  TH1F* barelhcaliso = (TH1F*) outfile->Get(selection+"sctbar_elhcaliso");
  TH1F* bareltkiso = (TH1F*) outfile->Get(selection+"sctbar_eltkiso");
  TH1F* barelr9 = (TH1F*) outfile->Get(selection+"sctbar_elr9");
  TH1F* barelsmin = (TH1F*) outfile->Get(selection+"sctbar_elsmin");
  TH1F* barelsmaj = (TH1F*) outfile->Get(selection+"sctbar_elsmaj");
  
  TH1F* ecelpt = (TH1F*) outfile->Get(selection+"sctec_elpt");
  TH1F* ecelm = (TH1F*) outfile->Get(selection+"sctec_elm");
  TH1F* eceld0 = (TH1F*) outfile->Get(selection+"sctec_eld0");
  TH1F* ecellog10d0 = (TH1F*) outfile->Get(selection+"sctec_ellog10d0");
  TH1F* eceldz = (TH1F*) outfile->Get(selection+"sctec_eldz");
  TH1F* eceldetain = (TH1F*) outfile->Get(selection+"sctec_eldetain");
  TH1F* eceldphiin = (TH1F*) outfile->Get(selection+"sctec_eldphiin");
  TH1F* ecelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctec_elsigmaietaieta");
  TH1F* ecelhoe = (TH1F*) outfile->Get(selection+"sctec_elhoe");
  TH1F* ecelooemoop = (TH1F*) outfile->Get(selection+"sctec_elooemoop");
  TH1F* ecelmhits = (TH1F*) outfile->Get(selection+"sctec_elmhits");
  TH1F* ecelcharge = (TH1F*) outfile->Get(selection+"sctec_elcharge");
  TH1F* ecelecaliso = (TH1F*) outfile->Get(selection+"sctec_elecaliso");
  TH1F* ecelhcaliso = (TH1F*) outfile->Get(selection+"sctec_elhcaliso");
  TH1F* eceltkiso = (TH1F*) outfile->Get(selection+"sctec_eltkiso");
  TH1F* ecelr9 = (TH1F*) outfile->Get(selection+"sctec_elr9");
  TH1F* ecelsmin = (TH1F*) outfile->Get(selection+"sctec_elsmin");
  TH1F* ecelsmaj = (TH1F*) outfile->Get(selection+"sctec_elsmaj");


  elmult->Fill(elidx.size());
    //highest pt amongst all scouting electrons in this evt
  float lead_pt = 0.0;
    //highest pt amongst'd all scouting electrons passing the new cut in this event
  float passlead_pt = 0.0;
    //highest pt amongst'd've all scouting electrons whomst'd've passed the old cut in this event
  float oldpasslead_pt = 0.0;
  float my_rho = (*rho)->at(0);
  for(unsigned int ctr=0; ctr<elidx.size(); ctr++) {
    float this_eta = (*ele_eta)->at(elidx[ctr]);
    float this_pt = (*ele_pt)->at(elidx[ctr]);
    elpt->Fill(this_pt);
    //update leading pt
    lead_pt = this_pt>lead_pt? this_pt : lead_pt;
    eleta->Fill(this_eta);
    elphi->Fill((*ele_phi)->at(elidx[ctr]));
    //fill the histogram for all scouting electrons (regardless if genmatched or nah) as a function of eta 
    all_sct->Fill(this_pt);  
    eta_sct->Fill(this_eta); 
    //get energy of the electron (deposited in the ECAL)
    vector<float> energyMatrix = (*ele_enemat)->at(elidx[ctr]);
    float energy_sum = 0.;
    for(float engy : energyMatrix) {
        energy_sum += engy;
    }
    //trying to use reconstructed particle E instead of ECAL deposits??
    TLorentzVector tlv;
    tlv.SetPtEtaPhiM( this_pt , this_eta, (*ele_phi)->at(elidx[ctr]) , .0005 );
    float energy_reco = tlv.E();
    EvsE->Fill(energy_sum, energy_reco);
    if (*(*n_rho) < 1 ) {
        cout << "Error: nrho 0 but electrons??" << endl;
        continue;
    }
    float energy_SC = energy_sum;
    //apparently we need to use the reco energy instead of the ECAL one :/
    energy_sum = energy_reco;
    float had_energy = energy_sum * (*ele_hoe)->at(elidx[ctr]);

    //fill RhovsH only in the barrel for now!
    //  also exclude any 0 values!!
    if(fabs(this_eta)<2.4 && had_energy > 0.05){
        RhovsH->Fill(my_rho, had_energy);
    }
    //fill in the HminusRho/E vs. E 2d histogram
    if(fabs(this_eta)<2.4 && had_energy > 0.05){
        float hmrho = (had_energy - Crho*my_rho)/energy_sum;
        //cout << "HmRho=" << hmrho << "; E= " << energy_sum << endl;
        //HmRhoVsE->Fill(energy_sum, hmrho);
        HmRhoVsE->Fill(energy_SC, hmrho);
    }
    //see if it passes the old/new cut
    if( had_energy / energy_sum < 0.2) {
        oldpass_sct->Fill(this_pt);
        oldpasslead_pt = this_pt>oldpasslead_pt? this_pt : oldpasslead_pt;
    }
    if( had_energy < 0.2*energy_sum + h0 ) {
        pass_sct->Fill(this_pt);
        passlead_pt = this_pt>passlead_pt? this_pt : oldpasslead_pt;
    }
    //fill these histos with energy (not pt)
    if( fabs(this_eta) < 1 ) {
        mid_sct->Fill(energy_sum);
    }
    else if( fabs(this_eta) < 1.479) {
        out_sct->Fill(energy_sum);
    }
    else {
        ecs_sct->Fill(energy_sum);
    }
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el, el2;
      el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M());
    }

    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      //cout << "ERROR: barrel electrons still here!!!!!" << endl;
      //exit(0);
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      //barelm->Fill((*ele_m)->at(elidx[ctr]));
      bareld0->Fill((*ele_d0)->at(elidx[ctr]));
      barellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      bareldz->Fill((*ele_dz)->at(elidx[ctr]));
      bareldetain->Fill((*ele_detain)->at(elidx[ctr]));
      bareldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      barelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      barelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      barelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      barelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      barelcharge->Fill((*ele_charge)->at(elidx[ctr]));
      barelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      barelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      bareltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      barelr9->Fill((*ele_r9)->at(elidx[ctr]));
      barelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      barelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }

    else {
      //cout << "ERROR: endcap electrons still here!!!!!" << endl;
      //exit(0);
      ecelpt->Fill((*ele_pt)->at(elidx[ctr]));
      //ecelm->Fill((*ele_m)->at(elidx[ctr]));
      eceld0->Fill((*ele_d0)->at(elidx[ctr]));
      ecellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      eceldz->Fill((*ele_dz)->at(elidx[ctr]));
      eceldetain->Fill((*ele_detain)->at(elidx[ctr]));
      eceldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      ecelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      ecelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      ecelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      ecelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      ecelcharge->Fill((*ele_charge)->at(elidx[ctr]));
      ecelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      ecelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      eceltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      ecelr9->Fill((*ele_r9)->at(elidx[ctr]));
      ecelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      ecelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }
  } // End of main electron for loop
  alllead_sct->Fill(lead_pt);
  oldpasslead_sct->Fill(oldpasslead_pt);
  passlead_sct->Fill(passlead_pt);
  
//eta_sig, eta_old_passed_sig, eta_new_passed_sig, eta_sct, eta_bkg
  for(unsigned int ctr=0; ctr<elidx_gm.size(); ctr++) {
    //do I need this check??
    if(goodmucount < 2) break;
    if(elidx_gm[ctr] == -1) continue;
    float gen_eta = (*ele_eta)->at(elidx_gm[ctr]);
    //exclude electrons in the low-efficiency region b/t EB and EE.
    if(fabs(gen_eta) > 1.44 && fabs(gen_eta) < 1.56) {
        continue;
    }
    ++ngm;
    ++ngmevt;
    vector<float> energyMatrix = (*ele_enemat)->at(elidx_gm[ctr]);
    float Electron_pt = (*ele_pt)->at(elidx_gm[ctr]);
    //float energy_sum = 0.;
    //for(float engy : energyMatrix) {
    //    energy_sum += engy;
    //}
    //trying reco e E instead of ECAL deposits to see if this works...
    TLorentzVector tlv;
    tlv.SetPtEtaPhiM( Electron_pt , (*ele_eta)->at(elidx_gm[ctr]) , (*ele_phi)->at(elidx_gm[ctr]) , .0005 );
    float energy_sum = tlv.E();
    //H is H/E * E
    float had_energy = energy_sum * (*ele_hoe)->at(elidx_gm[ctr]);
    if(hVe.size() > 0 && hVe[0] != 0) {
        hVe[0]->Fill(energy_sum, had_energy);
    }
    float genElectron_pt = (*genpart_pt)->at(ctr);
    if(hVe.size() > 1 && Electron_pt < 5) {
        hVe[1]->Fill(energy_sum, had_energy);
    }
    else if(hVe.size() > 2 && Electron_pt < 8) {
        hVe[2]->Fill(energy_sum, had_energy);
    }
    else if(hVe.size() > 3 && Electron_pt < 10) {
        hVe[3]->Fill(energy_sum, had_energy);
    }
    else if(hVe.size() > 4 && Electron_pt < 15) {
        hVe[4]->Fill(energy_sum, had_energy);
    }
    else if(hVe.size() > 5 && Electron_pt < 20) {
        hVe[5]->Fill(energy_sum, had_energy);
    }
    else if(hVe.size() > 6) {
        hVe[6]->Fill(energy_sum, had_energy);
    }
    if(hVe.size() > 7 && Electron_pt < 10) {
        hVe[7]->Fill(energy_sum, had_energy);
    }

    //float newcut = 23 + 0.2*energy_sum;
    float newcut = h0 + 0.2*energy_sum;
    if(had_energy < newcut) {
        //passed the cut
        ++npassed;
        ++npassedevt;
        passed_sig->Fill(genElectron_pt);
        eta_new_passed_sig->Fill(gen_eta);
    }
    //now that this cut is implemented in the menu, the else *should* never be executed.
    else {
        //cout << endl << endl << endl;
        //cout << "***FAILED H<.2E+16 cut*** pt: " << Electron_pt << "; eta: " << gen_eta << "; phi: " << (*ele_phi)->at(elidx_gm[ctr]) << "; H: " << had_energy << "; E: " << energy_sum << endl;
        //cout << endl << endl << endl;
    }
    float rhocut = h0 + h1*energy_sum + my_rho*Crho;
    if(had_energy < rhocut){
        rho_passed_sig->Fill(genElectron_pt);
    }

    if(had_energy / energy_sum < 0.2) {
        old_passed_sig->Fill(genElectron_pt);
        if(genElectron_pt > 5.0) {
            old_passed_sig_vRho->Fill(my_rho);
        }
        eta_old_passed_sig->Fill(gen_eta);
        ++noldpassed;
        ++noldpassedevt;
    }
    all_sig->Fill(genElectron_pt);
    if(genElectron_pt > 5.0) {
        all_sig_vRho->Fill(my_rho);
    }
    eta_sig->Fill(gen_eta);
    
    float sct_eta = (*ele_eta)->at(elidx_gm[ctr]);
    if( fabs(sct_eta) < 1 ) {
        mid_gen->Fill(genElectron_pt);
    }
    else if( fabs(sct_eta) < 1.479) {
        out_gen->Fill(genElectron_pt);
    }
    else {
        ecs_gen->Fill(genElectron_pt);
    }
  } //end genmatched electron loop
  //cout << "elidx_notgm size: " << elidx_notgm.size() << endl;
  for(unsigned int ctr=0; ctr<elidx_notgm.size(); ctr++) {
    float gen_eta = (*ele_eta)->at(elidx_notgm[ctr]);
    //exclude electrons in the low-efficiency region b/t EB,EE
    if(fabs(gen_eta) > 1.44 && fabs(gen_eta) < 1.56) {
        continue;
    }
    vector<float> energyMatrix = (*ele_enemat)->at(elidx_notgm[ctr]);
    float Electron_pt = (*ele_pt)->at(elidx_notgm[ctr]);
    //float energy_sum = 0.;
    //for(float engy : energyMatrix) {
    //    energy_sum += engy;
    //}
    //H is H/E * E
    TLorentzVector tlv;
    tlv.SetPtEtaPhiM( Electron_pt , (*ele_eta)->at(elidx_notgm[ctr]) , (*ele_phi)->at(elidx_notgm[ctr]) , .0005 );
    float energy_sum = tlv.E();
    //
    float had_energy = energy_sum * (*ele_hoe)->at(elidx_notgm[ctr]);
    float newcut = h0 + 0.2*energy_sum;
    //cout << "notgenmatched " << ctr << ": E, H, pt: " << energy_sum << ", " << had_energy << ", " << Electron_pt << endl;
    if(had_energy < newcut) {
        //passed the cut
        passed_bkg->Fill(Electron_pt);
        npassedevt++;
    }
    if(had_energy / energy_sum < 0.2) {
        old_passed_bkg->Fill(Electron_pt);
        noldpassedevt++;
    }
    all_bkg->Fill(Electron_pt);
    eta_bkg->Fill(gen_eta);
  } //end notgenmatched electron loop 

  elmult_gm->Fill(ngmevt);
  elmult_passed->Fill(npassedevt);
  elmult_oldpassed->Fill(noldpassedevt);
  
}  

// Function to add a set of histograms for gen particles
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  //new hists
  all1dhists.push_back(new TH1F(selection+"sct_elmult_gm","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elmult_passed","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elmult_oldpassed","N e",50,-5,45));
  //
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",600,-3,3));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",10000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  //all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",100,0,1.0));
  all1dhists.push_back(new TH1F(selection+"sctbar_elooemoop","1/E-1/p",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctec_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctec_eldetain","#Delta#eta_{in}",10000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_eldphiin","#Delta#phi_{in}",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  //all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",100,0,1.0));
  all1dhists.push_back(new TH1F(selection+"sctec_elooemoop","1/E-1/p",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctec_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctec_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmaj","smaj",1500,-0.1,1.4));
}

// Function to sort the indices based on a factor (Usually pT)
void data_robustanalyzer::sort(int* idx, TTreeReaderValue<std::vector<float>> *factor, int n) {
  for(unsigned int i=0; i<n; i++) {
    for(unsigned int j=i+1; j<n; j++) {
      if((*factor)->at((*(idx+j)))>(*factor)->at((*(idx+i)))) { // Sort in decreasing value of factor
	double temp = *(idx+i);
	*(idx+i) = *(idx+j);
	*(idx+j) = temp;
      }
    }
  }
}

// Function to find an electron pair in the Z mass window 80<M<100
pair<int,int> data_robustanalyzer::inZwindow(vector<int> elidx) {
  
  if(elidx.size()<2) { // Require atleast two electrons
    return make_pair(-1,-1);
  }

  pair<int,int> foundZels = make_pair(-1,-1);
  TLorentzVector el1, el2;
  for(unsigned int ctr1=0; ctr1<elidx.size(); ctr1++) {
    for(unsigned int ctr2=ctr1+1; ctr2<elidx.size(); ctr2++) {
      el1.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr1]),(*ele_eta)->at(elidx[ctr1]),(*ele_phi)->at(elidx[ctr1]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      if((el1+el2).M()>75 && (el1+el2).M()<100) {
	foundZels = make_pair(elidx[ctr1],elidx[ctr2]);
	break;
      }
    }
  }

  return foundZels;
}

//bpg trying to add this plot.
//void plotHoEvsPt(TString selection, vector<int> signalElectrons) {
void data_robustanalyzer::fillHoEvsPt(TH2F* hoevpt, vector<int> signalElectrons) {
   // TH2F* hoevpt = new TH2F("hoevspt", "H/E vs. #p_T", 1000, 0, 1000, 100, 0, 1.0);
    for(int j: signalElectrons) {
        if(j == -1) continue;
        //cout << "Filling in electron pt and H/E: " << (*ele_pt)->at(j) << ", " << (*ele_hoe)->at(j) << endl;
        hoevpt->Fill((*ele_pt)->at(j), (*ele_hoe)->at(j), 1.0);
    }
} //end fillHoEvsPt

//find the best match for the gen particles
//pair<int,int> data_robustanalyzer::genMatch(vector<int> electronIdxs, bool domatch) {
vector<int> data_robustanalyzer::genMatch(vector<int> electronIdxs, bool domatch, vector<TH1D*> deta_all, vector<TH1D*> qdphi_all, vector<TH1D*> deta_gmd, vector<TH1D*> qdphi_gmd, TH1F* pt_res_sct, TH1F* e_ratio_sctgen, TH1F* e_ratio_sctgenl20) {
    //pair for the electrons matched to the gen electrons
    //pair<int,int> matchedPair;
    int nparts = (int)((*genpart_pt)->size());
    int nele = (int)((*ele_pt)->size());
    //cout << "ngenparts: " << nparts << ", nele: " << nele << endl;
    vector<int> matchedPair(nparts, -1);
    //matchedPair.first = -1;
    //matchedPair.second = -1;
    if(!domatch) return matchedPair;
    //first make sure we have the needed gen information
    //if((*genpart_pt)->size() != 2) {
    //    cout << "Error: need exactly 2 gen particles but we have " << (*genpart_pt)->size() << endl;
    //    //throw "NGenParticleError";
    //    return matchedPair;
    //}
    //vector of (2) gen particles--
    //set each to -1 once it's matched to make sure we can't match to the same one twice!
    vector<int> genParts(nparts, 0);
    //genParts[0] = 0;
    //genParts[1] = 1;
    for(int j=0; j < nparts; j++) {
        genParts[j] = j;
    }
    for(int el : electronIdxs) {
        double my_eta = (*ele_eta)->at(el);
        //for(int gp = 0; gp < 2; gp++) {
        for(int gp = 0; gp < nparts; gp++) {
            //electrons only!!
            if(abs((*genpart_pdg)->at(gp)) != 11) continue;
            if(!(*genpart_isFinalState)->at(gp)) continue;
            if(genParts[gp] == -1) continue;
            //cout << "el: " << el << "; gp: " << gp << endl;
            double diffeta = abs( my_eta - (*genpart_eta)->at(gp) );
            double diffphi = abs( (*ele_phi)->at(el) - (*genpart_phi)->at(gp));
            //fill the deta before hist
            deta_all[0]->Fill(diffeta);
            
            TLorentzVector vec_el, vec_gen;
            vec_gen.SetPtEtaPhiM((*genpart_pt)->at(gp),(*genpart_eta)->at(gp),(*genpart_phi)->at(gp),0.000511);
            vec_el.SetPtEtaPhiM((*ele_pt)->at(el),(*ele_eta)->at(el),(*ele_phi)->at(el),0.000511);
            double qdiffphi = ((*genpart_pdg)->at(gp)/fabs((*genpart_pdg)->at(gp)))*(vec_gen.DeltaPhi(vec_el));
            float pt = (*genpart_pt)->at(gp);
            //fill the dphi before hist
            qdphi_all[0]->Fill(qdiffphi); 
            if(pt>=2 && pt <= 5){
                deta_all[1]->Fill(diffeta);
                qdphi_all[1]->Fill(qdiffphi);
            }
            else if(pt > 5 && pt <= 8) { 
                deta_all[2]->Fill(diffeta);
                qdphi_all[2]->Fill(qdiffphi);
            }
            else if(pt > 8 && pt <= 10) {
                deta_all[3]->Fill(diffeta);
                qdphi_all[3]->Fill(qdiffphi);
            }
            else if(pt > 10 && pt <= 15) { 
                deta_all[4]->Fill(diffeta);
                qdphi_all[4]->Fill(qdiffphi);
            }
            else if(pt > 15 && pt <= 20){
                deta_all[5]->Fill(diffeta);
                qdphi_all[5]->Fill(qdiffphi);
            }
            else {
                //cout << "super big pt: " << pt << endl;
                deta_all[6]->Fill(diffeta);
                qdphi_all[6]->Fill(qdiffphi);
            }
            //cout << "pt: " << pt << endl;
            bool foundMatch = false;
            if( abs(my_eta) < 1.479 ) {
                //if(diffeta<0.1 && qdiffphi<0.15 && qdiffphi>-0.25) {
                //new cutoffs
                if(diffeta<0.2 && qdiffphi<0.05 && qdiffphi>-0.3) {
                    foundMatch = true;
                } //end found match 
            } //end central eta region
            else {
                //if(diffeta<0.05 && qdiffphi<0.1 && qdiffphi>-0.15) {
                //new cutoffs
                if(diffeta<0.1 && qdiffphi<0.05 && qdiffphi>-0.2) {
                    foundMatch = true;
                } //end found match
            } //end not central eta region
            if(foundMatch) {
                //cout << "found match!" << endl;
                //if(gp == 0) matchedPair.first = el;
                //else matchedPair.second = el;
                matchedPair[gp] = el;
                genParts[gp] = -1; 
                //fill the genmatched deta, qdphi hists
                deta_gmd[0]->Fill(diffeta);
                qdphi_gmd[0]->Fill(qdiffphi);
                if(pt>=2 && pt <= 5){
                    deta_gmd[1]->Fill(diffeta);
                    qdphi_gmd[1]->Fill(qdiffphi);
                }
                else if(pt > 5 && pt <= 8) { 
                    deta_gmd[2]->Fill(diffeta);
                    qdphi_gmd[2]->Fill(qdiffphi);
                }
                else if(pt > 8 && pt <= 10) {
                    deta_gmd[3]->Fill(diffeta);
                    qdphi_gmd[3]->Fill(qdiffphi);
                }
                else if(pt > 10 && pt <= 15) { 
                    deta_gmd[4]->Fill(diffeta);
                    qdphi_gmd[4]->Fill(qdiffphi);
                }
                else if(pt > 15 && pt <= 20){
                    deta_gmd[5]->Fill(diffeta);
                    qdphi_gmd[5]->Fill(qdiffphi);
                }
                else {
                    deta_gmd[6]->Fill(diffeta);
                    qdphi_gmd[6]->Fill(qdiffphi);
                }
                //fill the pt resolution histogram
                //cout << "filling the ptres hist." << endl;
                float pt_res = (*ele_pt)->at(el) - pt;
                pt_res_sct->Fill(pt_res);
                //now find the energy ratio
                vector<float> energyMatrix = (*ele_enemat)->at(el);
                float energy_sum = 0.;
                for(float engy : energyMatrix) {
                    energy_sum += engy;
                }
                TLorentzVector gen4vec;
                gen4vec.SetPtEtaPhiM( (*genpart_pt)->at(gp), (*genpart_eta)->at(gp), (*genpart_phi)->at(gp), 0.000511 );
                TLorentzVector reco4vec;
                reco4vec.SetPtEtaPhiM( (*ele_pt)->at(el), (*ele_eta)->at(el), (*ele_phi)->at(el), 0.000511 );
                float reco_energy = reco4vec.Energy();
                float gen_energy = gen4vec.Energy();
                //float pt_ratio = (*ele_pt)->at(el) / pt;
                //float e_ratio = energy_sum / gen_energy;
                //cout << "energy matrix sum: " << energy_sum << "; reco energy: " << reco_energy << endl;
                float e_ratio = reco_energy / gen_energy;
                //pt_ratio_sctgen->Fill(pt_ratio);
                e_ratio_sctgen->Fill(e_ratio);
                //ONLY fill the l20 if the pt of the (gen) electron is < 20 GeV!! (DoubleElectronGun pT goes too high).
                //if(pt < 20) {
                if(pt < 10) {
                    e_ratio_sctgenl20->Fill(e_ratio);
                }
            } //end foundMatch
        } //end loop over (2) gen particles
    } //end electrons loop
    //cout << "ngen: " << nparts << ", matched: ";
    //for(int j : matchedPair) cout << j << ", ";
    //cout << endl;
    return matchedPair;
} //end genMatch function

