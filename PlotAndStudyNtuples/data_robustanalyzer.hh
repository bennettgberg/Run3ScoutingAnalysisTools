#ifndef DATA_ROBUSTANALYZER_H
#define DATA_ROBUSTANALYZER_H

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

class data_robustanalyzer {

 public:
  data_robustanalyzer(TString, TString, bool);
  ~data_robustanalyzer();
  
  void analyzersinglefile(int);
  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
  pair<int,int> inZwindow(vector<int>);
  
 private:

  bool isDiEl;
  
  TTreeReader* tree;
  TTreeReaderArray<float> *bsx;
  TTreeReaderArray<float> *bsy;
  TTreeReaderArray<float> *bsz;
  TTreeReaderValue<UInt_t> *n_ele;
  TTreeReaderValue<vector<float>> *ele_pt;
  TTreeReaderValue<vector<float>> *ele_eta;
  TTreeReaderValue<vector<float>> *ele_phi;
  TTreeReaderValue<vector<float>> *ele_m;
  TTreeReaderValue<vector<float>> *ele_d0;
  TTreeReaderValue<vector<float>> *ele_dz;
  TTreeReaderValue<vector<float>> *ele_detain;
  TTreeReaderValue<vector<float>> *ele_dphiin;
  TTreeReaderValue<vector<float>> *ele_sigmaietaieta;
  TTreeReaderValue<vector<float>> *ele_hoe;
  TTreeReaderValue<vector<float>> *ele_ooemoop;
  TTreeReaderValue<vector<int>> *ele_mhits;
  TTreeReaderValue<vector<int>> *ele_charge;
  TTreeReaderValue<vector<float>> *ele_ecaliso;
  TTreeReaderValue<vector<float>> *ele_hcaliso;
  TTreeReaderValue<vector<float>> *ele_tkiso;
  TTreeReaderValue<vector<float>> *ele_r9;
  TTreeReaderValue<vector<float>> *ele_smin;
  TTreeReaderValue<vector<float>> *ele_smaj;
  TTreeReaderValue<vector<unsigned int>> *ele_seedid;
  TTreeReaderValue<vector<vector<float>>> *ele_enemat;
  TTreeReaderValue<vector<vector<float>>> *ele_timmat;
  
  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
