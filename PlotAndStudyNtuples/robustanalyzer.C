#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores, bool isDrellYan, bool isMonteCarlo){

  nC = numCores;
  isDY = isDrellYan;
  isMC = isMonteCarlo;

  cout<<"Initializing for file: "<<filename<<endl;
  TChain* chain = new TChain("mmtree/tree");
  chain->Add(filename);

  tree = new TTreeReader(chain);
  if(isMC) {
    n_gen = new TTreeReaderValue<unsigned int>((*tree), "n_genpart");
    gen_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_pdg");
    gen_pt = new TTreeReaderValue<vector<float>>((*tree), "genpart_pt");
    gen_eta = new TTreeReaderValue<vector<float>>((*tree), "genpart_eta");
    gen_phi = new TTreeReaderValue<vector<float>>((*tree), "genpart_phi");
    gen_vx = new TTreeReaderValue<vector<float>>((*tree), "genpart_vx");
    gen_vy = new TTreeReaderValue<vector<float>>((*tree), "genpart_vy");
    gen_vz = new TTreeReaderValue<vector<float>>((*tree), "genpart_vz");
    gen_nmoms = new TTreeReaderValue<vector<int>>((*tree), "genpart_nmoms");
    gen_mompdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_mompdg");
    gen_islastcopy = new TTreeReaderValue<vector<bool>>((*tree), "genpart_isLastCopy");
  }
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
  n_rho = new TTreeReaderValue<UInt_t>((*tree), "n_rhoval");
  rho = new TTreeReaderValue<vector<float>>((*tree), "rho");

  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
robustanalyzer::~robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = tree->GetEntries(true);
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = nC;
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  endevent = endevent<totEntries?endevent:totEntries;
  tree->SetEntriesRange(beginevent, endevent);
  cout<<"Processing events in range: [ "<<beginevent<<" , "<<endevent<<" )"<<endl;
  int event = beginevent-1;

  // Count events passing certain selections
  int nosel=0, vetoidsel=0, looseidsel=0, mediumidsel=0, tightidsel=0, tightnoneidsel=0, nonetightidsel=0, tightlooseidsel=0, loosetightidsel=0, tightmediumidesl=0, mediumtightidsel=0;

  // Define the histograms
  if(isMC) {
    addgenhist("noselgen");
    //addgenmchhist("noselgenAnosel");
  }
  addhist("nosel");
  addhist("nosel_Zwind_");
  addhist("leadepemsel");
  addhist("leadepem_only_sel");
  addhist("leadepem_ptgt20_sel");
  addhist("leadepem_Zwind_sel");
  addhist("leadepem_Zwindptgt20_sel");
  addhist("vetosel");
  addhist("loosesel");
  addhist("loosesel_Zwind_");
  addhist("mediumsel");
  addhist("mediumsel_Zwind_");
  addhist("tightsel");
  addhist("tightsel_Zwind_");
  addhist("tightsel_SideBand1_");
  addhist("tightsel_SideBand2_");
  addhist("tightsel_SideBand3_");
  addhist("tightnonesel");
  addhist("tightnonesel_Zwind_");
  addhist("tightnonesel_SideBand1_");
  addhist("tightnonesel_SideBand2_");
  addhist("tightnonesel_SideBand3_");
  addhist("nonetightsel");
  addhist("nonetightsel_Zwind_");
  addhist("nonetightsel_SideBand1_");
  addhist("nonetightsel_SideBand2_");
  addhist("nonetightsel_SideBand3_");
  addhist("tightloosesel");
  addhist("tightloosesel_Zwind_");
  addhist("tightloosesel_SideBand1_");
  addhist("tightloosesel_SideBand2_");
  addhist("tightloosesel_SideBand3_");
  addhist("loosetightsel");
  addhist("loosetightsel_Zwind_");
  addhist("loosetightsel_SideBand1_");
  addhist("loosetightsel_SideBand2_");
  addhist("loosetightsel_SideBand3_");
  addhist("tightmediumsel");
  addhist("tightmediumsel_Zwind_");
  addhist("tightmediumsel_SideBand1_");
  addhist("tightmediumsel_SideBand2_");
  addhist("tightmediumsel_SideBand3_");
  addhist("mediumtightsel");
  addhist("mediumtightsel_Zwind_");
  addhist("mediumtightsel_SideBand1_");
  addhist("mediumtightsel_SideBand2_");
  addhist("mediumtightsel_SideBand3_");
  
  // vector of electron indices
  vector<int> noselgenidx;
  vector<int> noselelidx;
  vector<int> noselZwindelidx;
  vector<int> leadepemselelidx;
  vector<int> leadepemonlyselelidx;
  vector<int> leadepemptgt20selelidx;
  vector<int> leadepemselZwindelidx;
  vector<int> leadepemZwindptgt20selelidx;
  vector<int> vetoselelidx;
  vector<int> looseselelidx;
  vector<int> looseselZwindelidx;
  vector<int> mediumselelidx;
  vector<int> mediumselZwindelidx;
  vector<int> tightselelidx;
  vector<int> tightselZwindelidx;
  vector<int> tightselSideBand1elidx;
  vector<int> tightselSideBand2elidx;
  vector<int> tightselSideBand3elidx;
  vector<int> tightnoneselelidx;
  vector<int> tightnoneselZwindelidx;
  vector<int> tightnoneselSideBand1elidx;
  vector<int> tightnoneselSideBand2elidx;
  vector<int> tightnoneselSideBand3elidx;
  vector<int> nonetightselelidx;
  vector<int> nonetightselZwindelidx;
  vector<int> nonetightselSideBand1elidx;
  vector<int> nonetightselSideBand2elidx;
  vector<int> nonetightselSideBand3elidx;
  vector<int> tightlooseselelidx;
  vector<int> tightlooseselZwindelidx;
  vector<int> tightlooseselSideBand1elidx;
  vector<int> tightlooseselSideBand2elidx;
  vector<int> tightlooseselSideBand3elidx;
  vector<int> loosetightselelidx;
  vector<int> loosetightselZwindelidx;
  vector<int> loosetightselSideBand1elidx;
  vector<int> loosetightselSideBand2elidx;
  vector<int> loosetightselSideBand3elidx;
  vector<int> tightmediumselelidx;
  vector<int> tightmediumselZwindelidx;
  vector<int> tightmediumselSideBand1elidx;
  vector<int> tightmediumselSideBand2elidx;
  vector<int> tightmediumselSideBand3elidx;
  vector<int> mediumtightselelidx;
  vector<int> mediumtightselZwindelidx;
  vector<int> mediumtightselSideBand1elidx;
  vector<int> mediumtightselSideBand2elidx;
  vector<int> mediumtightselSideBand3elidx;
    
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%100000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    bool noselcond = false;
    // Loop on Gen particles to select good gen electrons
    if(isMC) {
      
      for(unsigned int gen_ctr=0; gen_ctr<(*(*n_gen)); gen_ctr++) {

	noselcond = true;
	noselcond *= isDY?(abs((*gen_pdg)->at(gen_ctr))==11 && (*gen_mompdg)->at(gen_ctr)==23):((*gen_islastcopy)->at(gen_ctr));
	if(noselcond) noselgenidx.push_back(gen_ctr);
	
      } // End of loop on gen electrons
      
      if(noselgenidx.size()>=2) fillgenhistinevent("noselgen", noselgenidx);
    }
    
    if(isMC && isDY) {
      if(noselgenidx.size()<2) continue;
    }
    
    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine

    if((*(*n_ele))<0) throw "Error!!! Wrong technical event processing. Negative number of electrons in event.";

    bool vetoselcond = false;
    bool looseselcond = false;
    bool mediumselcond = false;
    bool tightselcond = false;
    
    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];

      noselelidx.push_back(elidx);

      // Get the energy of this electron
      double ele_energy=0;
      for(unsigned int ctr=0; ctr<((*ele_enemat)->at(elidx)).size(); ctr++) {
	//if(((*ele_enemat)->at(elidx))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	if(((*ele_enemat)->at(elidx))[ctr]>0) {
	  ele_energy += ((*ele_enemat)->at(elidx))[ctr];
	}
      }
      
      vetoselcond = true;
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_sigmaietaieta)->at(elidx)<0.0126):((*ele_sigmaietaieta)->at(elidx)<0.0457);
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_detain)->at(elidx)<0.00463):((*ele_detain)->at(elidx)<0.00814);
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_dphiin)->at(elidx)<0.148):((*ele_dphiin)->at(elidx)<0.19);
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_hoe)->at(elidx)<(0.05+(1.16/ele_energy)/*+(0.0324*(*(*rho))/ele_energy)*/))
	:((*ele_hoe)->at(elidx)<(0.05+(2.54/ele_energy)/*+(0.183*(*(*rho))/ele_energy)*/));
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_tkiso)->at(elidx)<(0.198+(0.506/(*ele_pt)->at(elidx))))
	:((*ele_tkiso)->at(elidx)<(0.203+(0.963/(*ele_pt)->at(elidx))));
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_ooemoop)->at(elidx)<0.209):((*ele_detain)->at(elidx)<0.132);
      vetoselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_mhits)->at(elidx)<=2):((*ele_mhits)->at(elidx)<=3);
      if(vetoselcond) vetoselelidx.push_back(elidx);
      
      looseselcond = true;
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_sigmaietaieta)->at(elidx)<0.0112):((*ele_sigmaietaieta)->at(elidx)<0.0425);
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_detain)->at(elidx)<0.00377):((*ele_detain)->at(elidx)<0.00674);
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_dphiin)->at(elidx)<0.0884):((*ele_dphiin)->at(elidx)<0.169);
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_hoe)->at(elidx)<(0.05+(1.16/ele_energy)/*+(0.0324*(*(*rho))/ele_energy)*/))
	:((*ele_hoe)->at(elidx)<(0.0441+(2.54/ele_energy)/*+(0.183*(*(*rho))/ele_energy)*/));
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_tkiso)->at(elidx)<(0.112+(0.506/(*ele_pt)->at(elidx))))
	:((*ele_tkiso)->at(elidx)<(0.108+(0.963/(*ele_pt)->at(elidx))));
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_ooemoop)->at(elidx)<0.193):((*ele_detain)->at(elidx)<0.111);
      looseselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_mhits)->at(elidx)<=1):((*ele_mhits)->at(elidx)<=1);
      if(looseselcond) looseselelidx.push_back(elidx);
      
      mediumselcond = true;
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_sigmaietaieta)->at(elidx)<0.0106):((*ele_sigmaietaieta)->at(elidx)<0.0387);
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_detain)->at(elidx)<0.0032):((*ele_detain)->at(elidx)<0.00632);
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_dphiin)->at(elidx)<0.0547):((*ele_dphiin)->at(elidx)<0.0394);
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_hoe)->at(elidx)<(0.046+(1.16/ele_energy)/*+(0.0324*(*(*rho))/ele_energy)*/))
	:((*ele_hoe)->at(elidx)<(0.0275+(2.52/ele_energy)/*+(0.183*(*(*rho))/ele_energy)*/));
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_tkiso)->at(elidx)<(0.0478+(0.506/(*ele_pt)->at(elidx))))
	:((*ele_tkiso)->at(elidx)<(0.0658+(0.963/(*ele_pt)->at(elidx))));
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_ooemoop)->at(elidx)<0.184):((*ele_detain)->at(elidx)<0.0721);
      mediumselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_mhits)->at(elidx)<=1):((*ele_mhits)->at(elidx)<=1);
      if(mediumselcond) mediumselelidx.push_back(elidx);
      
      tightselcond = true;
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_sigmaietaieta)->at(elidx)<0.0104):((*ele_sigmaietaieta)->at(elidx)<0.0353);
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_detain)->at(elidx)<0.00255):((*ele_detain)->at(elidx)<0.00501);
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_dphiin)->at(elidx)<0.022):((*ele_dphiin)->at(elidx)<0.0236);
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_hoe)->at(elidx)<(0.026+(1.15/ele_energy)/*+(0.0324*(*(*rho))/ele_energy)*/))
	:((*ele_hoe)->at(elidx)<(0.0188+(2.06/ele_energy)/*+(0.183*(*(*rho))/ele_energy)*/));
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_tkiso)->at(elidx)<(0.0287+(0.506/(*ele_pt)->at(elidx))))
	:((*ele_tkiso)->at(elidx)<(0.0455+(0.963/(*ele_pt)->at(elidx))));
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_ooemoop)->at(elidx)<0.159):((*ele_detain)->at(elidx)<0.0197);
      tightselcond *= (abs((*ele_eta)->at(elidx))<1.479)?((*ele_mhits)->at(elidx)<=1):((*ele_mhits)->at(elidx)<=1);
      if(tightselcond) tightselelidx.push_back(elidx);
      
    }// End of loop on electrons in the event loop

    // At least 2 opp. charged electrons in event
    if(noselelidx.size()>=2) {
      if((*ele_charge)->at(noselelidx[0])*(*ele_charge)->at(noselelidx[1])<0) {

        leadepemselelidx = noselelidx;
	leadepemonlyselelidx.push_back(leadepemselelidx[0]);
	leadepemonlyselelidx.push_back(leadepemselelidx[1]);
	if((*ele_pt)->at(leadepemselelidx[0])>20 && (*ele_pt)->at(leadepemselelidx[1])>20) {
	  leadepemptgt20selelidx.push_back(leadepemselelidx[0]);
	  leadepemptgt20selelidx.push_back(leadepemselelidx[1]);
	}
	
	bool leadloosecond=false, leadmediumcond=false, leadtightcond=false;
	bool subleadloosecond=false, subleadmediumcond=false, subleadtightcond=false;
	
	double lead_ele_energy=0;
	for(unsigned int ctr=0; ctr<((*ele_enemat)->at(noselelidx[0])).size(); ctr++) {
	  //if(((*ele_enemat)->at(noselelidx[0]))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	  lead_ele_energy += ((*ele_enemat)->at(noselelidx[0]))[ctr];
	}
	
	double sublead_ele_energy=0;
	for(unsigned int ctr=0; ctr<((*ele_enemat)->at(noselelidx[1])).size(); ctr++) {
	  //if(((*ele_enemat)->at(noselelidx[1]))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	  sublead_ele_energy += ((*ele_enemat)->at(noselelidx[1]))[ctr];
	}
	
	leadloosecond = true;
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[0])<0.0112):((*ele_sigmaietaieta)->at(noselelidx[0])<0.0425);
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_detain)->at(noselelidx[0])<0.00377):((*ele_detain)->at(noselelidx[0])<0.00674);
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_dphiin)->at(noselelidx[0])<0.0884):((*ele_dphiin)->at(noselelidx[0])<0.169);
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_hoe)->at(noselelidx[0])<(0.05+(1.16/lead_ele_energy)/*+(0.0324*(*(*rho))/lead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[0])<(0.0441+(2.54/lead_ele_energy)/*+(0.183*(*(*rho))/lead_ele_energy)*/));
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_tkiso)->at(noselelidx[0])<(0.112+(0.506/(*ele_pt)->at(noselelidx[0]))))
	  :((*ele_tkiso)->at(noselelidx[0])<(0.108+(0.963/(*ele_pt)->at(noselelidx[0]))));
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_ooemoop)->at(noselelidx[0])<0.193):((*ele_detain)->at(noselelidx[0])<0.111);
	leadloosecond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_mhits)->at(noselelidx[0])<=1):((*ele_mhits)->at(noselelidx[0])<=1);
      
	leadmediumcond = true;
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[0])<0.0106):((*ele_sigmaietaieta)->at(noselelidx[0])<0.0387);
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_detain)->at(noselelidx[0])<0.0032):((*ele_detain)->at(noselelidx[0])<0.00632);
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_dphiin)->at(noselelidx[0])<0.0547):((*ele_dphiin)->at(noselelidx[0])<0.0394);
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_hoe)->at(noselelidx[0])<(0.046+(1.16/lead_ele_energy)/*+(0.0324*(*(*rho))/lead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[0])<(0.0275+(2.52/lead_ele_energy)/*+(0.183*(*(*rho))/lead_ele_energy)*/));
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_tkiso)->at(noselelidx[0])<(0.0478+(0.506/(*ele_pt)->at(noselelidx[0]))))
	  :((*ele_tkiso)->at(noselelidx[0])<(0.0658+(0.963/(*ele_pt)->at(noselelidx[0]))));
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_ooemoop)->at(noselelidx[0])<0.184):((*ele_detain)->at(noselelidx[0])<0.0721);
	leadmediumcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_mhits)->at(noselelidx[0])<=1):((*ele_mhits)->at(noselelidx[0])<=1);
      
	leadtightcond = true;
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[0])<0.0104):((*ele_sigmaietaieta)->at(noselelidx[0])<0.0353);
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_detain)->at(noselelidx[0])<0.00255):((*ele_detain)->at(noselelidx[0])<0.00501);
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_dphiin)->at(noselelidx[0])<0.022):((*ele_dphiin)->at(noselelidx[0])<0.0236);
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_hoe)->at(noselelidx[0])<(0.026+(1.15/lead_ele_energy)/*+(0.0324*(*(*rho))/lead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[0])<(0.0188+(2.06/lead_ele_energy)/*+(0.183*(*(*rho))/lead_ele_energy)*/));
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_tkiso)->at(noselelidx[0])<(0.0287+(0.506/(*ele_pt)->at(noselelidx[0]))))
	  :((*ele_tkiso)->at(noselelidx[0])<(0.0455+(0.963/(*ele_pt)->at(noselelidx[0]))));
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_ooemoop)->at(noselelidx[0])<0.159):((*ele_detain)->at(noselelidx[0])<0.0197);
	leadtightcond *= (abs((*ele_eta)->at(noselelidx[0]))<1.479)?((*ele_mhits)->at(noselelidx[0])<=1):((*ele_mhits)->at(noselelidx[0])<=1);

	subleadloosecond = true;
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[1])<0.0112):((*ele_sigmaietaieta)->at(noselelidx[1])<0.0425);
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_detain)->at(noselelidx[1])<0.00377):((*ele_detain)->at(noselelidx[1])<0.00674);
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_dphiin)->at(noselelidx[1])<0.0884):((*ele_dphiin)->at(noselelidx[1])<0.169);
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_hoe)->at(noselelidx[1])<(0.05+(1.16/sublead_ele_energy)/*+(0.0324*(*(*rho))/sublead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[1])<(0.0441+(2.54/sublead_ele_energy)/*+(0.183*(*(*rho))/sublead_ele_energy)*/));
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_tkiso)->at(noselelidx[1])<(0.112+(0.506/(*ele_pt)->at(noselelidx[1]))))
	  :((*ele_tkiso)->at(noselelidx[1])<(0.108+(0.963/(*ele_pt)->at(noselelidx[1]))));
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_ooemoop)->at(noselelidx[1])<0.193):((*ele_detain)->at(noselelidx[1])<0.111);
	subleadloosecond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_mhits)->at(noselelidx[1])<=1):((*ele_mhits)->at(noselelidx[1])<=1);
      
	subleadmediumcond = true;
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[1])<0.0106):((*ele_sigmaietaieta)->at(noselelidx[1])<0.0387);
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_detain)->at(noselelidx[1])<0.0032):((*ele_detain)->at(noselelidx[1])<0.00632);
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_dphiin)->at(noselelidx[1])<0.0547):((*ele_dphiin)->at(noselelidx[1])<0.0394);
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_hoe)->at(noselelidx[1])<(0.046+(1.16/sublead_ele_energy)/*+(0.0324*(*(*rho))/sublead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[1])<(0.0275+(2.52/sublead_ele_energy)/*+(0.183*(*(*rho))/sublead_ele_energy)*/));
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_tkiso)->at(noselelidx[1])<(0.0478+(0.506/(*ele_pt)->at(noselelidx[1]))))
	  :((*ele_tkiso)->at(noselelidx[1])<(0.0658+(0.963/(*ele_pt)->at(noselelidx[1]))));
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_ooemoop)->at(noselelidx[1])<0.184):((*ele_detain)->at(noselelidx[1])<0.0721);
	subleadmediumcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_mhits)->at(noselelidx[1])<=1):((*ele_mhits)->at(noselelidx[1])<=1);
      
	subleadtightcond = true;
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_sigmaietaieta)->at(noselelidx[1])<0.0104):((*ele_sigmaietaieta)->at(noselelidx[1])<0.0353);
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_detain)->at(noselelidx[1])<0.00255):((*ele_detain)->at(noselelidx[1])<0.00501);
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_dphiin)->at(noselelidx[1])<0.022):((*ele_dphiin)->at(noselelidx[1])<0.0236);
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_hoe)->at(noselelidx[1])<(0.026+(1.15/sublead_ele_energy)/*+(0.0324*(*(*rho))/sublead_ele_energy)*/))
	  :((*ele_hoe)->at(noselelidx[1])<(0.0188+(2.06/sublead_ele_energy)/*+(0.183*(*(*rho))/sublead_ele_energy)*/));
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_tkiso)->at(noselelidx[1])<(0.0287+(0.506/(*ele_pt)->at(noselelidx[1]))))
	  :((*ele_tkiso)->at(noselelidx[1])<(0.0455+(0.963/(*ele_pt)->at(noselelidx[1]))));
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_ooemoop)->at(noselelidx[1])<0.159):((*ele_detain)->at(noselelidx[1])<0.0197);
	subleadtightcond *= (abs((*ele_eta)->at(noselelidx[1]))<1.479)?((*ele_mhits)->at(noselelidx[1])<=1):((*ele_mhits)->at(noselelidx[1])<=1);

	if(leadtightcond) {
	  tightnoneselelidx.push_back(noselelidx[0]);
	  tightnoneselelidx.push_back(noselelidx[1]);
	  if(subleadloosecond) {
	    tightlooseselelidx.push_back(noselelidx[0]);
	    tightlooseselelidx.push_back(noselelidx[1]);
	  }
	  if(subleadmediumcond) {
	    tightmediumselelidx.push_back(noselelidx[0]);
	    tightmediumselelidx.push_back(noselelidx[1]);
	  }
	}

	if(subleadtightcond) {
	  nonetightselelidx.push_back(noselelidx[0]);
	  nonetightselelidx.push_back(noselelidx[1]);
	  if(leadloosecond) {
	    loosetightselelidx.push_back(noselelidx[0]);
	    loosetightselelidx.push_back(noselelidx[1]);
	  }
	  if(leadmediumcond) {
	    mediumtightselelidx.push_back(noselelidx[0]);
	    mediumtightselelidx.push_back(noselelidx[1]);
	  }
	}
      }	
    } // Atleast 2 opp. charged e in event

    // At least 2 opp. charged electrons with no sel in event
    if(noselelidx.size()>=2) {
      if((*ele_charge)->at(noselelidx[0])*(*ele_charge)->at(noselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(noselelidx[0]),(*ele_eta)->at(noselelidx[0]),(*ele_phi)->at(noselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(noselelidx[1]),(*ele_eta)->at(noselelidx[1]),(*ele_phi)->at(noselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  noselZwindelidx.push_back(noselelidx[0]);
	  noselZwindelidx.push_back(noselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with nosel in event

    if(leadepemonlyselelidx.size()>=2) {
      TLorentzVector leadel, subleadel;
      leadel.SetPtEtaPhiM((*ele_pt)->at(leadepemonlyselelidx[0]),(*ele_eta)->at(leadepemonlyselelidx[0]),(*ele_phi)->at(leadepemonlyselelidx[0]),0.0005);
      subleadel.SetPtEtaPhiM((*ele_pt)->at(leadepemonlyselelidx[1]),(*ele_eta)->at(leadepemonlyselelidx[1]),(*ele_phi)->at(leadepemonlyselelidx[1]),0.0005);
      if(inZwind(leadel, subleadel)) {
	leadepemselZwindelidx.push_back(noselelidx[0]);
	leadepemselZwindelidx.push_back(noselelidx[1]);
	if((*ele_pt)->at(leadepemonlyselelidx[0])>20 && (*ele_pt)->at(leadepemonlyselelidx[1])>20) {
	  leadepemZwindptgt20selelidx.push_back(noselelidx[0]);
	  leadepemZwindptgt20selelidx.push_back(noselelidx[1]);
	}
      }
    }
    
    // At least 2 opp. charged electrons with loose ID in event
    if(looseselelidx.size()>=2) {
      if((*ele_charge)->at(looseselelidx[0])*(*ele_charge)->at(looseselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(looseselelidx[0]),(*ele_eta)->at(looseselelidx[0]),(*ele_phi)->at(looseselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(looseselelidx[1]),(*ele_eta)->at(looseselelidx[1]),(*ele_phi)->at(looseselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  looseselZwindelidx.push_back(looseselelidx[0]);
	  looseselZwindelidx.push_back(looseselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with loose ID in event
    
    // At least 2 opp. charged electrons with medium ID in event
    if(mediumselelidx.size()>=2) {
      if((*ele_charge)->at(mediumselelidx[0])*(*ele_charge)->at(mediumselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(mediumselelidx[0]),(*ele_eta)->at(mediumselelidx[0]),(*ele_phi)->at(mediumselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(mediumselelidx[1]),(*ele_eta)->at(mediumselelidx[1]),(*ele_phi)->at(mediumselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  mediumselZwindelidx.push_back(mediumselelidx[0]);
	  mediumselZwindelidx.push_back(mediumselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with medium ID in event
    
    // At least 2 opp. charged electrons with tight ID in event
    if(tightselelidx.size()>=2) {
      if((*ele_charge)->at(tightselelidx[0])*(*ele_charge)->at(tightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightselelidx[0]),(*ele_eta)->at(tightselelidx[0]),(*ele_phi)->at(tightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightselelidx[1]),(*ele_eta)->at(tightselelidx[1]),(*ele_phi)->at(tightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightselZwindelidx.push_back(tightselelidx[0]);
	  tightselZwindelidx.push_back(tightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightselSideBand1elidx.push_back(tightselelidx[0]);
	  tightselSideBand1elidx.push_back(tightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightselSideBand2elidx.push_back(tightselelidx[0]);
	  tightselSideBand2elidx.push_back(tightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightselSideBand3elidx.push_back(tightselelidx[0]);
	  tightselSideBand3elidx.push_back(tightselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with tight ID in event
    
    // Lead tight Sub-lead none in event
    if(tightnoneselelidx.size()>=2) {
      if((*ele_charge)->at(tightnoneselelidx[0])*(*ele_charge)->at(tightnoneselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightnoneselelidx[0]),(*ele_eta)->at(tightnoneselelidx[0]),(*ele_phi)->at(tightnoneselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightnoneselelidx[1]),(*ele_eta)->at(tightnoneselelidx[1]),(*ele_phi)->at(tightnoneselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightnoneselZwindelidx.push_back(tightnoneselelidx[0]);
	  tightnoneselZwindelidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightnoneselSideBand1elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand1elidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightnoneselSideBand2elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand2elidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightnoneselSideBand3elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand3elidx.push_back(tightnoneselelidx[1]);
	}
      }
    } // Lead tight Sub-lead none in event

    // Lead none Sub-lead tight in event
    if(nonetightselelidx.size()>=2) {
      if((*ele_charge)->at(nonetightselelidx[0])*(*ele_charge)->at(nonetightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(nonetightselelidx[0]),(*ele_eta)->at(nonetightselelidx[0]),(*ele_phi)->at(nonetightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(nonetightselelidx[1]),(*ele_eta)->at(nonetightselelidx[1]),(*ele_phi)->at(nonetightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  nonetightselZwindelidx.push_back(nonetightselelidx[0]);
	  nonetightselZwindelidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  nonetightselSideBand1elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand1elidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  nonetightselSideBand2elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand2elidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  nonetightselSideBand3elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand3elidx.push_back(nonetightselelidx[1]);
	}
      }
    } // Lead none Sub-lead tight in event

    // Lead tight Sub-lead loose in event
    if(tightlooseselelidx.size()>=2) {
      if((*ele_charge)->at(tightlooseselelidx[0])*(*ele_charge)->at(tightlooseselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightlooseselelidx[0]),(*ele_eta)->at(tightlooseselelidx[0]),(*ele_phi)->at(tightlooseselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightlooseselelidx[1]),(*ele_eta)->at(tightlooseselelidx[1]),(*ele_phi)->at(tightlooseselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightlooseselZwindelidx.push_back(tightlooseselelidx[0]);
	  tightlooseselZwindelidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightlooseselSideBand1elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand1elidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightlooseselSideBand2elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand2elidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightlooseselSideBand3elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand3elidx.push_back(tightlooseselelidx[1]);
	}
      }
    } // Lead tight Sub-lead loose in event

    // Lead loose Sub-lead tight in event
    if(loosetightselelidx.size()>=2) {
      if((*ele_charge)->at(loosetightselelidx[0])*(*ele_charge)->at(loosetightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(loosetightselelidx[0]),(*ele_eta)->at(loosetightselelidx[0]),(*ele_phi)->at(loosetightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(loosetightselelidx[1]),(*ele_eta)->at(loosetightselelidx[1]),(*ele_phi)->at(loosetightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  loosetightselZwindelidx.push_back(loosetightselelidx[0]);
	  loosetightselZwindelidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  loosetightselSideBand1elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand1elidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  loosetightselSideBand2elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand2elidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  loosetightselSideBand3elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand3elidx.push_back(loosetightselelidx[1]);
	}
      }
    } // Lead loose Sub-lead tight in event

    // Lead tight Sub-lead medium in event
    if(tightmediumselelidx.size()>=2) {
      if((*ele_charge)->at(tightmediumselelidx[0])*(*ele_charge)->at(tightmediumselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightmediumselelidx[0]),(*ele_eta)->at(tightmediumselelidx[0]),(*ele_phi)->at(tightmediumselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightmediumselelidx[1]),(*ele_eta)->at(tightmediumselelidx[1]),(*ele_phi)->at(tightmediumselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightmediumselZwindelidx.push_back(tightmediumselelidx[0]);
	  tightmediumselZwindelidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightmediumselSideBand1elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand1elidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightmediumselSideBand2elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand2elidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightmediumselSideBand3elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand3elidx.push_back(tightmediumselelidx[1]);
	}
      }
    } // Lead tight Sub-lead medium in event

    // Lead medium Sub-lead tight in event
    if(mediumtightselelidx.size()>=2) {
      if((*ele_charge)->at(mediumtightselelidx[0])*(*ele_charge)->at(mediumtightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(mediumtightselelidx[0]),(*ele_eta)->at(mediumtightselelidx[0]),(*ele_phi)->at(mediumtightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(mediumtightselelidx[1]),(*ele_eta)->at(mediumtightselelidx[1]),(*ele_phi)->at(mediumtightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  mediumtightselZwindelidx.push_back(mediumtightselelidx[0]);
	  mediumtightselZwindelidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  mediumtightselSideBand1elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand1elidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  mediumtightselSideBand2elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand2elidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  mediumtightselSideBand3elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand3elidx.push_back(mediumtightselelidx[1]);
	}
      }
    } // Lead medium Sub-lead tight in event

    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel", noselelidx);
    fillhistinevent("nosel_Zwind_", noselZwindelidx);
    fillhistinevent("leadepemsel", leadepemselelidx);
    fillhistinevent("leadepem_only_sel", leadepemonlyselelidx);
    fillhistinevent("leadepem_ptgt20_sel", leadepemptgt20selelidx);
    fillhistinevent("leadepem_Zwind_sel", leadepemselZwindelidx);
    fillhistinevent("leadepem_Zwindptgt20_sel", leadepemZwindptgt20selelidx);
    fillhistinevent("vetosel", vetoselelidx);
    fillhistinevent("loosesel", looseselelidx);
    fillhistinevent("loosesel_Zwind_", looseselZwindelidx);
    fillhistinevent("mediumsel", mediumselelidx);
    fillhistinevent("mediumsel_Zwind_", mediumselZwindelidx);
    fillhistinevent("tightsel", tightselelidx);
    fillhistinevent("tightsel_Zwind_", tightselZwindelidx);
    fillhistinevent("tightsel_SideBand1_", tightselSideBand1elidx);
    fillhistinevent("tightsel_SideBand2_", tightselSideBand2elidx);
    fillhistinevent("tightsel_SideBand3_", tightselSideBand3elidx);
    fillhistinevent("tightnonesel", tightnoneselelidx);
    fillhistinevent("tightnonesel_Zwind_", tightnoneselZwindelidx);
    fillhistinevent("tightnonesel_SideBand1_", tightnoneselSideBand1elidx);
    fillhistinevent("tightnonesel_SideBand2_", tightnoneselSideBand2elidx);
    fillhistinevent("tightnonesel_SideBand3_", tightnoneselSideBand3elidx);
    fillhistinevent("nonetightsel", nonetightselelidx);
    fillhistinevent("nonetightsel_Zwind_", nonetightselZwindelidx);
    fillhistinevent("nonetightsel_SideBand1_", nonetightselSideBand1elidx);
    fillhistinevent("nonetightsel_SideBand2_", nonetightselSideBand2elidx);
    fillhistinevent("nonetightsel_SideBand3_", nonetightselSideBand3elidx);
    fillhistinevent("tightloosesel", tightlooseselelidx);
    fillhistinevent("tightloosesel_Zwind_", tightlooseselZwindelidx);
    fillhistinevent("tightloosesel_SideBand1_", tightlooseselSideBand1elidx);
    fillhistinevent("tightloosesel_SideBand2_", tightlooseselSideBand2elidx);
    fillhistinevent("tightloosesel_SideBand3_", tightlooseselSideBand3elidx);
    fillhistinevent("loosetightsel", loosetightselelidx);
    fillhistinevent("loosetightsel_Zwind_", loosetightselZwindelidx);
    fillhistinevent("loosetightsel_SideBand1_", loosetightselSideBand1elidx);
    fillhistinevent("loosetightsel_SideBand2_", loosetightselSideBand2elidx);
    fillhistinevent("loosetightsel_SideBand3_", loosetightselSideBand3elidx);
    fillhistinevent("tightmediumsel", tightmediumselelidx);
    fillhistinevent("tightmediumsel_Zwind_", tightmediumselZwindelidx);
    fillhistinevent("tightmediumsel_SideBand1_", tightmediumselSideBand1elidx);
    fillhistinevent("tightmediumsel_SideBand2_", tightmediumselSideBand2elidx);
    fillhistinevent("tightmediumsel_SideBand3_", tightmediumselSideBand3elidx);
    fillhistinevent("mediumtightsel", mediumtightselelidx);
    fillhistinevent("mediumtightsel_Zwind_", mediumtightselZwindelidx);
    fillhistinevent("mediumtightsel_SideBand1_", mediumtightselSideBand1elidx);
    fillhistinevent("mediumtightsel_SideBand2_", mediumtightselSideBand2elidx);
    fillhistinevent("mediumtightsel_SideBand3_", mediumtightselSideBand3elidx);
    //if(isMC && noselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAnosel", noselgenidx, noselelidx);

    // Clear all vector
    noselgenidx.clear();
    noselelidx.clear();
    noselZwindelidx.clear();
    leadepemselelidx.clear();
    leadepemonlyselelidx.clear();
    leadepemptgt20selelidx.clear();
    leadepemselZwindelidx.clear();
    leadepemZwindptgt20selelidx.clear();
    vetoselelidx.clear();
    looseselelidx.clear();
    looseselZwindelidx.clear();
    mediumselelidx.clear();
    mediumselZwindelidx.clear();
    tightselelidx.clear();
    tightselZwindelidx.clear();
    tightselSideBand1elidx.clear();
    tightselSideBand2elidx.clear();
    tightselSideBand3elidx.clear();
    tightnoneselelidx.clear();
    tightnoneselZwindelidx.clear();
    tightnoneselSideBand1elidx.clear();
    tightnoneselSideBand2elidx.clear();
    tightnoneselSideBand3elidx.clear();
    nonetightselelidx.clear();
    nonetightselZwindelidx.clear();
    nonetightselSideBand1elidx.clear();
    nonetightselSideBand2elidx.clear();
    nonetightselSideBand3elidx.clear();
    tightlooseselelidx.clear();
    tightlooseselZwindelidx.clear();
    tightlooseselSideBand1elidx.clear();
    tightlooseselSideBand2elidx.clear();
    tightlooseselSideBand3elidx.clear();
    loosetightselelidx.clear();
    loosetightselZwindelidx.clear();
    loosetightselSideBand1elidx.clear();
    loosetightselSideBand2elidx.clear();
    loosetightselSideBand3elidx.clear();
    tightmediumselelidx.clear();
    tightmediumselZwindelidx.clear();
    tightmediumselSideBand1elidx.clear();
    tightmediumselSideBand2elidx.clear();
    tightmediumselSideBand3elidx.clear();
    mediumtightselelidx.clear();
    mediumtightselZwindelidx.clear();
    mediumtightselSideBand1elidx.clear();
    mediumtightselSideBand2elidx.clear();
    mediumtightselSideBand3elidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<endl;
}

// Fill histograms for gen matching
void robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> elidx) {

  if(!isMC) throw "Error!!! Trying to fill gen info. for a non MC file.";
    
  // Variables before gen match
  TH1F* bardEta = (TH1F*) outfile->Get(selection+"genelsctbar_dEta");
  TH1F* barqdPhi = (TH1F*) outfile->Get(selection+"genelsctbar_qdPhi");
  TH1F* eedEta = (TH1F*) outfile->Get(selection+"genelsctee_dEta");
  TH1F* eeqdPhi = (TH1F*) outfile->Get(selection+"genelsctee_qdPhi");
  
  // Variables after gen match
  TH1F* mchbardEta = (TH1F*) outfile->Get(selection+"genelsctmchbar_dEta");
  TH1F* mchbarqdPhi = (TH1F*) outfile->Get(selection+"genelsctmchbar_qdPhi");
  TH1F* mcheedEta = (TH1F*) outfile->Get(selection+"genelsctmchee_dEta");
  TH1F* mcheeqdPhi = (TH1F*) outfile->Get(selection+"genelsctmchee_qdPhi");

  TH1F* genelmult = (TH1F*) outfile->Get(selection+"sctmchgenel_elmult");
  TH1F* genelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_elpt");
  TH1F* geneleta = (TH1F*) outfile->Get(selection+"sctmchgenel_eleta");
  TH1F* genelphi = (TH1F*) outfile->Get(selection+"sctmchgenel_elphi");
  TH1F* gendielM = (TH1F*) outfile->Get(selection+"sctmchgenel_dielM");

  TH1F* genmchsctelmult = (TH1F*) outfile->Get(selection+"genmchsct_elmult");
  TH1F* genmchsctelpt = (TH1F*) outfile->Get(selection+"genmchsct_elpt");
  TH1F* genmchscteleta = (TH1F*) outfile->Get(selection+"genmchsct_eleta");
  TH1F* genmchsctelphi = (TH1F*) outfile->Get(selection+"genmchsct_elphi");
  TH1F* genmchsctdielM = (TH1F*) outfile->Get(selection+"genmchsct_dielM");

  // Fill variables before gen match
  for(unsigned int sct=0; sct<elidx.size(); sct++) {
    for(unsigned int gen=0; gen<genidx.size(); gen++) {
      double charge = ((*gen_pdg)->at(genidx[gen]))/TMath::Abs((*gen_pdg)->at(genidx[gen]));
      if(TMath::Abs((*ele_eta)->at(elidx[sct]))<1.479) {
	bardEta->Fill((*ele_eta)->at(elidx[sct])-(*gen_eta)->at(genidx[gen]));
	barqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*gen_phi)->at(genidx[gen])));
      }
      else {
	eedEta->Fill((*ele_eta)->at(elidx[sct])-(*gen_eta)->at(genidx[gen]));
	eeqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*gen_phi)->at(genidx[gen])));
      }
    }
  }
  
  // Fill variables after gen match
  vector< pair<int,int> > sctelgenmch = diElecGenMatching(genidx, elidx);

  int gen1 = sctelgenmch[0].first;
  int el1 = sctelgenmch[0].second;
  int gen2 = sctelgenmch[1].first;
  int el2 = sctelgenmch[1].second;
  
  int genmchfound = 0;
  if(el1!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen1));
    geneleta->Fill((*gen_eta)->at(gen1));
    genelphi->Fill((*gen_phi)->at(gen1));
  }
  if(el2!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen2));
    geneleta->Fill((*gen_eta)->at(gen2));
    genelphi->Fill((*gen_phi)->at(gen2));
  }
  genelmult->Fill(genmchfound);
  if(genmchfound==2) {
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*gen_pt)->at(gen1),(*gen_eta)->at(gen1),(*gen_phi)->at(gen1),0.0005);
    elec2.SetPtEtaPhiM((*gen_pt)->at(gen2),(*gen_eta)->at(gen2),(*gen_phi)->at(gen2),0.0005);
    gendielM->Fill((elec1+elec2).M());
  }
  
  if(el1!=-1) {
    double charge = ((*gen_pdg)->at(gen1))/TMath::Abs((*gen_pdg)->at(gen1));
    if(TMath::Abs((*ele_eta)->at(el1))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el1)-(*gen_eta)->at(gen1));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*gen_phi)->at(gen1)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el1)-(*gen_eta)->at(gen1));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*gen_phi)->at(gen1)));
    }
    genmchsctelpt->Fill((*ele_pt)->at(el1));
    genmchscteleta->Fill((*ele_eta)->at(el1));
    genmchsctelphi->Fill((*ele_phi)->at(el1));
  }
  
  if(el2!=-1) {
    double charge = ((*gen_pdg)->at(gen2))/TMath::Abs((*gen_pdg)->at(gen2));
    if(TMath::Abs((*ele_eta)->at(el2))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el2)-(*gen_eta)->at(gen2));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*gen_phi)->at(gen2)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el2)-(*gen_eta)->at(gen2));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*gen_phi)->at(gen2)));
    }
    genmchsctelpt->Fill((*ele_pt)->at(el2));
    genmchscteleta->Fill((*ele_eta)->at(el2));
    genmchsctelphi->Fill((*ele_phi)->at(el2));
  }
  genmchsctelmult->Fill(genmchfound);
  if(genmchfound==2){
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*ele_pt)->at(el1),(*ele_eta)->at(el1),(*ele_phi)->at(el1),0.0005);
    elec2.SetPtEtaPhiM((*ele_pt)->at(el2),(*ele_eta)->at(el2),(*ele_phi)->at(el2),0.0005);
    genmchsctdielM->Fill((elec1+elec2).M());
  }
  
}

// Function to fill a set of histograms for gen particles
void robustanalyzer::fillgenhistinevent(TString selection, vector<int> genidx) {

  if(!isMC) throw "Error!!! Trying to fill gen info. for a non MC file.";
    
  if(genidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"gen_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"gen_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"gen_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"gen_elphi");

  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"gen_lead_elpt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"gen_lead_eleta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"gen_lead_elphi");

  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"gen_sublead_elpt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"gen_sublead_eleta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"gen_sublead_elphi");

  TH1F* dielM = (TH1F*) outfile->Get(selection+"gen_diel_M");
  TH1F* dieldeta = (TH1F*) outfile->Get(selection+"gen_diel_deta");
  TH1F* dieldphi = (TH1F*) outfile->Get(selection+"gen_diel_dphi");
  TH1F* dieldR = (TH1F*) outfile->Get(selection+"gen_diel_dR");

  unsigned int leadptpos=-1, subleadptpos=-1;
  
  elmult->Fill(genidx.size());
  for(unsigned int ctr=0; ctr<genidx.size(); ctr++) {
    elpt->Fill((*gen_pt)->at(genidx[ctr]));
    eleta->Fill((*gen_eta)->at(genidx[ctr]));
    elphi->Fill((*gen_phi)->at(genidx[ctr]));
    
    if(leadptpos==-1 || ( (*gen_pt)->at(genidx[ctr])>(*gen_pt)->at(leadptpos) )) {
      subleadptpos = leadptpos;
      leadptpos = genidx[ctr];
    }
    if(leadptpos!=genidx[ctr] &&
       (subleadptpos==-1 || ( (*gen_pt)->at(genidx[ctr])>(*gen_pt)->at(subleadptpos) ) )
       ) {
      subleadptpos = genidx[ctr];
    }
  }
  if(leadptpos!=-1) {
    leadelpt->Fill((*gen_pt)->at(leadptpos));
    leadeleta->Fill((*gen_eta)->at(leadptpos));
    leadelphi->Fill((*gen_phi)->at(leadptpos));
  }
  if(subleadptpos!=-1) {
    subleadelpt->Fill((*gen_pt)->at(subleadptpos));
    subleadeleta->Fill((*gen_eta)->at(subleadptpos));
    subleadelphi->Fill((*gen_phi)->at(subleadptpos));
  }
  if(leadptpos!=-1 && subleadptpos!=-1) {
    TLorentzVector el1, el2;
    el1.SetPtEtaPhiM((*gen_pt)->at(leadptpos),(*gen_eta)->at(leadptpos),(*gen_phi)->at(leadptpos),0.0005);
    el2.SetPtEtaPhiM((*gen_pt)->at(subleadptpos),(*gen_eta)->at(subleadptpos),(*gen_phi)->at(subleadptpos),0.0005);
    dielM->Fill((el1+el2).M());
    dieldeta->Fill((*gen_eta)->at(leadptpos)-(*gen_eta)->at(subleadptpos));
    dieldphi->Fill(el1.DeltaPhi(el2));
    dieldR->Fill(el1.DeltaR(el2));
  }
  else {
    cout<<"Warning!!! Unequal to two electron found in event."<<endl;
  }
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  if(elidx.size()==0) return;

  TH1F* rhohist = (TH1F*) outfile->Get(selection+"sct_rho");
  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  TH1F* leadsubleaddielM = (TH1F*) outfile->Get(selection+"sct_leadsublead_dielM");
  TH1F* leadsubleadqprod = (TH1F*) outfile->Get(selection+"sct_leadsublead_chargeprod");
  TH1F* leadbarsubleadbardielM = (TH1F*) outfile->Get(selection+"sct_leadbarsubleadbar_dielM");
  TH1F* leadecsubleadecdielM = (TH1F*) outfile->Get(selection+"sct_leadecsubleadec_dielM");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sctbar_elpt");
  TH1F* bareleta = (TH1F*) outfile->Get(selection+"sctbar_eleta");
  TH1F* barelphi = (TH1F*) outfile->Get(selection+"sctbar_elphi");
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
  TH1F* eceleta = (TH1F*) outfile->Get(selection+"sctec_eleta");
  TH1F* ecelphi = (TH1F*) outfile->Get(selection+"sctec_elphi");
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

  rhohist->Fill((*rho)->at(0));
  elmult->Fill(elidx.size());

  if(elidx.size()>=2) {
    TLorentzVector leadel, subleadel;
    leadel.SetPtEtaPhiM((*ele_pt)->at(elidx[0]),(*ele_eta)->at(elidx[0]),(*ele_phi)->at(elidx[0]),0.0005);
    subleadel.SetPtEtaPhiM((*ele_pt)->at(elidx[1]),(*ele_eta)->at(elidx[1]),(*ele_phi)->at(elidx[1]),0.0005);
    leadsubleaddielM->Fill((leadel+subleadel).M());
    leadsubleadqprod->Fill((*ele_charge)->at(elidx[0])*(*ele_charge)->at(elidx[1]));
    if(abs((*ele_eta)->at(elidx[0]))<1.479 && abs((*ele_eta)->at(elidx[1]))<1.479)  {
      leadbarsubleadbardielM->Fill((leadel+subleadel).M());
    }
    if(abs((*ele_eta)->at(elidx[0]))>1.479 && abs((*ele_eta)->at(elidx[1]))>1.479){
      leadecsubleadecdielM->Fill((leadel+subleadel).M());
    }
  }

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
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      bareleta->Fill((*ele_eta)->at(elidx[ctr]));
      barelphi->Fill((*ele_phi)->at(elidx[ctr]));
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
      ecelpt->Fill((*ele_pt)->at(elidx[ctr]));
      eceleta->Fill((*ele_eta)->at(elidx[ctr]));
      ecelphi->Fill((*ele_phi)->at(elidx[ctr]));
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

}  

// Add histograms for gen matching
void robustanalyzer::addgenmchhist(TString selection) {
  
  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";
    
  // Variables before gen match
  all1dhists.push_back(new TH1F(selection+"genelsctbar_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctbar_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  
  // Variables after gen match
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));

  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elpt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_eleta","gen e #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elphi","gen e #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_dielM","M(e,e)",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"genmchsct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"genmchsct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_dielM","all M(e,e)",1000,-10,990));
}

// Function to add a set of histograms for gen electrons
void robustanalyzer::addgenhist(TString selection) {

  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";
    
  all1dhists.push_back(new TH1F(selection+"gen_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_elphi","#phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_lead_elpt","e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_lead_eleta","e_{1} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_lead_elphi","e_{1} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_sublead_elpt","e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_sublead_eleta","e_{1} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_sublead_elphi","e_{1} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_diel_M","all M(e,e)",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_diel_deta","#Delta#eta(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dphi","#Delta#phi(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dR","#Delta R(e_{1}, e_{2})",1000,0,10));

}

// Function to add a set of histograms for scouting electrons
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_rho","rho",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_dielM","M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_chargeprod","q_{e1}*q_{e2}",3,-1,1));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbar_dielM","bar.e_{1} M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadec_dielM","ec. e_{1} M(e_{1},e_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctbar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctec_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctec_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elooemoop","1/E-1/p",20000,0,0.2));
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
void robustanalyzer::sort(int* idx, TTreeReaderValue<std::vector<float>> *factor, int n) {
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

// Function to do gen matching for to electrons
vector< pair<int,int> > robustanalyzer::diElecGenMatching(vector<int> genidx, vector<int> sctelidx) {

  if(genidx.size()!=2) {
    throw "Error! Analysis code only suitable for 2 gen electrons";
  }
  
  if(!isMC) {
    throw "Error! Cannot do gen matching. Not MC file.";
  }
  
  // Counter for the total no.of gen matches
  // and 1 history variable for the first gen matched position
  // Logic is not made to work with more than 2 gen el.
  vector< pair<int,int> > *gensctelmch = new vector< pair<int,int> >;
  gensctelmch->push_back(make_pair(genidx[0],-1));
  gensctelmch->push_back(make_pair(genidx[1],-1));
  
  // Find the sctelidx with the best angular match
  for(int elidx : sctelidx) {
    // Loop over the gen particles
    for(auto it=gensctelmch->begin(); it!=gensctelmch->end(); it++) {

      int genidx = (*it).first;
      int mchsctelidx = (*it).second;
      if(genidx==-1 || mchsctelidx!=-1) continue;

      double diffeta = abs((*ele_eta)->at(elidx)-(*gen_eta)->at(genidx)); 
      TLorentzVector vecsctel, vecgen;
      vecgen.SetPtEtaPhiM((*gen_pt)->at(genidx),(*gen_eta)->at(genidx),(*gen_phi)->at(genidx),0.0005);
      vecsctel.SetPtEtaPhiM((*ele_pt)->at(elidx),(*ele_eta)->at(elidx),(*ele_phi)->at(elidx),0.0005);
      double qdiffphi = ( (*gen_pdg)->at(genidx)/abs((*gen_pdg)->at(genidx)) )*( vecgen.DeltaPhi(vecsctel) );
      // Condition for gen matching
      if(abs((*ele_eta)->at(elidx))<1.479) {
	if(diffeta<0.2 && qdiffphi<0.05 && qdiffphi>-0.3) {
	  (*it).first = genidx;
	  (*it).second = elidx;
	}
      }
      else {
	if(diffeta<0.1 && qdiffphi<0.05 && qdiffphi>-0.2) {
	  (*it).first = genidx;
	  (*it).second = elidx;
	}
      }
    } // End of gen loop
    
  } // End of unseeded egamma object loop
  
  return (*gensctelmch);
}

// Return if selection in Z window
bool robustanalyzer::inZwind(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if(abs(lead.Eta())<1.479 || abs(sublead.Eta())<1.479) {
    if(invm>86 && invm<97) {
      return true;
    }
    else {
      return false;
    }
  }
  else {
    if(invm>79 && invm<98) {
      return true;
    }
    else {
      return false;
    }
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand1(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>40 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand2(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>45 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand3(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>50 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}
