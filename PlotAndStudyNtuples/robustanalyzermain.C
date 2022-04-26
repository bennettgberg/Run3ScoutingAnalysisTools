#include "data_robustanalyzer.hh"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {

  stringstream ss;
  ss << argv[1];
  int cnt;
  ss >> cnt;
  int cnt2;
  stringstream ss3;
  ss3 << argv[3];
  ss3 >> cnt2;
  stringstream ss1;
  //ss1<<"hists_DYToLLM50_"<<cnt<<".root";
  ss1<<"hists_"<<argv[2]<<"_"<<cnt<<".root";
  //data_robustanalyzer drana_DYToLLM50("./data/DYToLLM50_ScoutingSkim220127.root",ss1.str(), true);
  //                                    input file path/name, output name, isMC
  stringstream ss2;
  //ss2 << "../Analysis/test/" << argv[2] << "_0_" << ss.str() << ".root";
  ss2 << "/eos/user/b/bgreenbe/scouting/ntuples/" << argv[2] << "/" << argv[2] << "_0_" << ss.str() << ".root";
  data_robustanalyzer drana_DYToLLM50(ss2.str(),ss1.str(), true);
  drana_DYToLLM50.analyzersinglefile(cnt2);

  return -1;
}
