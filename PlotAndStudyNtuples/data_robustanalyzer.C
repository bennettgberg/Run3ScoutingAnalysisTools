#include "data_robustanalyzer.hh"
#include <iostream>
#include <sstream>
#include <numeric>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool isDoubleElectron){

  isDiEl = isDoubleElectron;
  
  TFile *inpfile = TFile::Open(filename,"READ");
  cout<<"Initializing for file: "<<filename<<endl;

  tree = new TTreeReader("mmtree/tree",inpfile);
  bsx = new TTreeReaderArray<float>((*tree), "beamspot_x");
  bsy = new TTreeReaderArray<float>((*tree), "beamspot_y");
  bsz = new TTreeReaderArray<float>((*tree), "beamspot_z");
  n_ele = new TTreeReaderValue<UInt_t>((*tree), "n_ele");
  ele_pt = new TTreeReaderValue<vector<float>>((*tree), "Electron_pt");
  ele_eta = new TTreeReaderValue<vector<float>>((*tree), "Electron_eta");
  ele_phi = new TTreeReaderValue<vector<float>>((*tree), "Electron_phi");
  ele_m = new TTreeReaderValue<vector<float>>((*tree), "Electron_m");
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
  ele_seedid = new TTreeReaderValue<vector<unsigned int>>((*tree), "Electron_seedid");
  ele_enemat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_energymatrix");
  ele_timmat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_timingmatrix");

  //now for the gen branches
  genpart_pt = new TTreeReaderValue<vector<float>>((*tree), "genpart_pt");
  genpart_eta = new TTreeReaderValue<vector<float>>((*tree), "genpart_eta");
  genpart_phi = new TTreeReaderValue<vector<float>>((*tree), "genpart_phi");
  genpart_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_pdg");
  
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
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  endevent = endevent<totEntries?endevent:totEntries;
  tree->SetEntriesRange(beginevent, endevent);
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
  TH2F* hoevpt = new TH2F("hoevspt", "H/E vs. #p_T for genmatched e-s", 1000, 0, 10, 4000, 0, 4.0);
  TH2F* hoevptBkg = new TH2F("hoevsptBkg", "H/E vs. #p_T for e-s NOT genmatched", 1000, 0, 10, 4000, 0, 4.0);
  //hists for efficiency as function of pT
    //all gen electrons
  TH1F* all_gen = new TH1F("allgen", "", 500, 0., 50.);
    //all gen-matched electrons
  TH1F* all_sig = new TH1F("allsig", "", 500, 0., 50.);
    //gen-matched electrons passing the new H/E cut
  TH1F* passed_sig = new TH1F("passedsig", "", 500, 0., 50.);
    //gen-matched electrons passing the old H/E cut H/E<.2
  TH1F* old_passed_sig = new TH1F("oldpassedsig", "", 500, 0., 50.);
    //all non-gen-matched electrons
  TH1F* all_bkg = new TH1F("allbkg", "", 500, 0., 50.);
    //non-gen-matched electrons passing the new H/E cut
  TH1F* passed_bkg = new TH1F("passedbkg", "", 500, 0., 50.);
    //non-gen-matched electrons passing the old H/E cut H/E<.2
  TH1F* old_passed_bkg = new TH1F("oldpassedbkg", "", 500, 0., 50.);
  //vector of TH2 plots comparing H and E
  // plot 0 is all pt, plot 1 is 2<pt<5, plot2 is 5<pt<8, plot3 is 8<pt<10 GeV.
  //hVe0: all
  //hVe1: pt<5; hVe2: 5<pt<8; hVe3: 8<pt<10; hVe4: 10<pt<15; hVe5:15<pt<20; hVe6:pt>20 GeV
   //GENMATCHED ONLY!
  vector<TH2F*> hVe(7);
  for(uint32_t i = 0; i < hVe.size(); i++) {
      stringstream hname;
      hname << "hVe" << i;
      hVe[i] = new TH2F(hname.str().c_str(), "H vs. E", 2000, 0, 200, 15000, 0, 1500.0);
      hVe[i]->GetXaxis()->SetTitle("HCal energy deposits (GeV)");
      hVe[i]->GetYaxis()->SetTitle("ECal energy deposits (GeV)");
      hVe[i]->SetMarkerColor(2);
  } //end i loop instantiating hVe.

  hoevpt->GetXaxis()->SetTitle("#p_T (GeV)");
  hoevpt->GetYaxis()->SetTitle("H/E");
  hoevpt->SetMarkerColor(2);
  hoevptBkg->SetMarkerColor(3);
  int nZeroe = 0;
  int nNonzeroe = 0;
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if((*(*n_ele))<0) cout<<"Error!!! Wrong technical event processing. Negative number of electrons in event."<<endl;;
    if((*(*n_ele))==0) { 
        nZeroe++;
        continue;
    }
    nNonzeroe++;
    
    //fill in the gen electrons histogram
    for(unsigned int genctr=0; genctr < (*genpart_pt)->size(); genctr++) {
        all_gen->Fill( (*genpart_pt)->at(genctr) );
    }
    
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
      if(abs(this_eta) < 1.479) continue;

      noselelidx.push_back(elidx);

      
    }// End of loop on electrons in the event loop

    // Event level selection on the electrons
    //pair<int,int> noselZwindels = inZwindow(noselelidx);
    pair<int,int> noselZwindels = genMatch(noselelidx);
    if(noselZwindels.first!=-1) {
      noselZwindelidx.push_back(noselZwindels.first);
    }
    if(noselZwindels.second != -1) {
      noselZwindelidx.push_back(noselZwindels.second);
    }
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
        if( elidx != noselZwindels.first && elidx != noselZwindels.second ){
            notZwindelidx.push_back(elidx);
        }
    } //end loop over all electrons

    if(noselelidx.size()>0) nosel++;
    //fillhistinevent("nosel", noselelidx, hVe);
    fillhistinevent("nosel", noselelidx, hVe, noselZwindelidx, notZwindelidx, all_sig, passed_sig, all_bkg, passed_bkg, old_passed_sig, old_passed_bkg);
    if(noselZwindelidx.size()>0) noselZwind++;
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
}

// Function to fill a set of histograms for gen particles
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> elidx, vector<TH2F*> hVe, vector<int> elidx_gm, vector<int> elidx_notgm, TH1F* all_sig, TH1F* passed_sig, TH1F* all_bkg, TH1F* passed_bkg, TH1F* old_passed_sig, TH1F* old_passed_bkg) {

  if(elidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sctbar_elpt");
  TH1F* barelm = (TH1F*) outfile->Get(selection+"sctbar_elm");
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
  for(unsigned int ctr=0; ctr<elidx.size(); ctr++) {
    elpt->Fill((*ele_pt)->at(elidx[ctr]));
    eleta->Fill((*ele_eta)->at(elidx[ctr]));
    elphi->Fill((*ele_phi)->at(elidx[ctr]));
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el, el2;
      el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M());
    }

    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      cout << "ERROR: barrel electrons still here!!!!!" << endl;
      exit(0);
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      barelm->Fill((*ele_m)->at(elidx[ctr]));
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
      ecelm->Fill((*ele_m)->at(elidx[ctr]));
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
  
  for(unsigned int ctr=0; ctr<elidx_gm.size(); ctr++) {
    vector<float> energyMatrix = (*ele_enemat)->at(elidx_gm[ctr]);
    float energy_sum = 0.;
    for(float engy : energyMatrix) {
        energy_sum += engy;
    }
    //H is H/E * E
    float had_energy = energy_sum * (*ele_hoe)->at(elidx_gm[ctr]);
    if(hVe.size() > 0 && hVe[0] != 0) {
        hVe[0]->Fill(had_energy, energy_sum);
    }
    float Electron_pt = (*ele_pt)->at(elidx_gm[ctr]);
    if(hVe.size() > 1 && Electron_pt < 5) {
        hVe[1]->Fill(had_energy, energy_sum);
    }
    else if(hVe.size() > 2 && Electron_pt < 8) {
        hVe[2]->Fill(had_energy, energy_sum);
    }
    else if(hVe.size() > 3 && Electron_pt < 10) {
        hVe[3]->Fill(had_energy, energy_sum);
    }
    else if(hVe.size() > 4 && Electron_pt < 15) {
        hVe[4]->Fill(had_energy, energy_sum);
    }
    else if(hVe.size() > 5 && Electron_pt < 20) {
        hVe[5]->Fill(had_energy, energy_sum);
    }
    else if(hVe.size() > 6) {
        hVe[6]->Fill(had_energy, energy_sum);
    }

    float newcut = 23 + 0.2*energy_sum;
    if(had_energy < newcut) {
        //passed the cut
        passed_sig->Fill(Electron_pt);
    }
    if(had_energy / energy_sum < 0.2) {
        old_passed_sig->Fill(Electron_pt);
    }
    all_sig->Fill(Electron_pt);
    
  } //end genmatched electron loop
  //cout << "elidx_notgm size: " << elidx_notgm.size() << endl;
  for(unsigned int ctr=0; ctr<elidx_notgm.size(); ctr++) {
    vector<float> energyMatrix = (*ele_enemat)->at(elidx_notgm[ctr]);
    float energy_sum = 0.;
    for(float engy : energyMatrix) {
        energy_sum += engy;
    }
    //H is H/E * E
    float had_energy = energy_sum * (*ele_hoe)->at(elidx_notgm[ctr]);
    float newcut = 23 + 0.2*energy_sum;
    float Electron_pt = (*ele_pt)->at(elidx_notgm[ctr]);
    //cout << "notgenmatched " << ctr << ": E, H, pt: " << energy_sum << ", " << had_energy << ", " << Electron_pt << endl;
    if(had_energy < newcut) {
        //passed the cut
        passed_bkg->Fill(Electron_pt);
    }
    if(had_energy / energy_sum < 0.2) {
        old_passed_bkg->Fill(Electron_pt);
    }
    all_bkg->Fill(Electron_pt);
  } //end notgenmatched electron loop
  
}  

// Function to add a set of histograms for gen particles
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
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
        //cout << "Filling in electron pt and H/E: " << (*ele_pt)->at(j) << ", " << (*ele_hoe)->at(j) << endl;
        hoevpt->Fill((*ele_pt)->at(j), (*ele_hoe)->at(j), 1.0);
    }
} //end fillHoEvsPt

//find the best match for the gen particles
pair<int,int> data_robustanalyzer::genMatch(vector<int> electronIdxs) {
    //pair for the electrons matched to the gen electrons
    pair<int,int> matchedPair;
    matchedPair.first = -1;
    matchedPair.second = -1;
    //first make sure we have the needed gen information
    if((*genpart_pt)->size() != 2) {
        cout << "Error: need exactly 2 gen particles but we have " << (*genpart_pt)->size() << endl;
        //throw "NGenParticleError";
        return matchedPair;
    }
    //vector of (2) gen particles--
    //set each to -1 once it's matched to make sure we can't match to the same one twice!
    vector<int> genParts(2);
    genParts[0] = 0;
    genParts[1] = 1;
    for(int el : electronIdxs) {
        double my_eta = (*ele_eta)->at(el);
        for(int gp = 0; gp < 2; gp++) {
            if(genParts[gp] == -1) continue;
            double diffeta = abs( my_eta - (*genpart_eta)->at(gp) );
            double diffphi = abs( (*ele_phi)->at(el) - (*genpart_phi)->at(gp));
            TLorentzVector vec_el, vec_gen;
            vec_gen.SetPtEtaPhiM((*genpart_pt)->at(gp),(*genpart_eta)->at(gp),(*genpart_phi)->at(gp),0.000511);
            vec_el.SetPtEtaPhiM((*ele_pt)->at(el),(*ele_eta)->at(el),(*ele_phi)->at(el),0.000511);
            double qdiffphi = ((*genpart_pdg)->at(gp)/(*genpart_pdg)->at(gp))*(vec_gen.DeltaPhi(vec_el));
        
            if( abs(my_eta) < 1.479 ) {
                if(diffeta<0.1 && qdiffphi<0.15 && qdiffphi>-0.25) {
                    if(gp == 0) matchedPair.first = el;
                    else matchedPair.second = el;
                    genParts[gp] = -1; 
                } //end found match 
            } //end central eta region
            else {
                if(diffeta<0.05 && qdiffphi<0.1 && qdiffphi>-0.15) {
                    if(gp == 0) matchedPair.first = el;
                    else matchedPair.second = el;
                    genParts[gp] = -1;  
                } //end found match
            } //end not central eta region
        } //end loop over (2) gen particles
    } //end electrons loop
    return matchedPair;
} //end genMatch function

