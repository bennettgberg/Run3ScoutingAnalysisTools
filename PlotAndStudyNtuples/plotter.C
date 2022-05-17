#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 1.09/19; // Data rate - 1.1 Hz, 22 events from parent trigger selection
double sig3cmsf = 1.0/15;
double sig30cmsf = 1.0/28;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
//TFile* sigDYToLLfile = TFile::Open("hists_DYToLL.root","READ");
TFile* sigDYToLLM50file = TFile::Open("hists_DYToLLM50.root","READ");
//TFile* sigDoubleElectronPt1To300file = TFile::Open("hists_DoubleElectron_Pt1To300.root","READ");
TFile* sigQCDPt20To30file = TFile::Open("hists_QCDPt20To30EmEnriched.root","READ");
TFile* sigQCDPt30To50file = TFile::Open("hists_QCDPt30To50EmEnriched.root","READ");
//TFile* sigQCDPt50To80file = TFile::Open("hists_QCD_Pt50To80_EMEnriched.root","READ");
//TFile* sigQCDPt15To20bcToEfile = TFile::Open("hists_QCD_Pt15To20bcToE_EMEnriched.root","READ");
//TFile* sigQCDPt20To30bcToEfile = TFile::Open("hists_QCD_Pt20To30bcToE_EMEnriched.root","READ");
//TFile* sigQCDPt30To80bcToEfile = TFile::Open("hists_QCD_Pt30To80bcToE_EMEnriched.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis") {

  // Check if the size of file vector is same as the cutnames
  if(file.size()!=cutname.size()) {
    cout<<"Error! Mismatching size of vectors"<<endl;
    return -1;
  }

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_";

  std::vector<TH1F*> allhists;
  for(unsigned int histctr=0; histctr<cutname.size(); histctr++) {
    TString cutwithvar = cutname[histctr]+"_"+var;
    allhists.push_back((TH1F*)file[histctr]->Get(cutwithvar));
    cutname[histctr] += "_"+((TString)file[histctr]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","");
  }
  
  // Get the title from histogram title
  TString dummytitle = "xaxis";
  TString histtitle = allhists[0]->GetTitle();
  TString xtitle = xaxistitle.CompareTo(dummytitle)?xaxistitle:histtitle;

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?(underflow?0:1):xbinlow;
  xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    if(underflow) {
      allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
      allhists[histctr]->SetBinError(xbinlow,err);
    }
    err = 0.0;
    if(overflow) {
      allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
      allhists[histctr]->SetBinError(xbinhigh,err);
    }

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    if(normalize) allhists[histctr]->GetYaxis()->SetTitle("normalized number of events (a.u.)");
    else allhists[histctr]->GetYaxis()->SetTitle("number of events");

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");

  return -1;
}

int efficiency(std::vector<TFile*> file, std::vector<TString> cutnames, int nbins=0, double *rebin=0, TString xaxistitle="p_{T} / GeV") {

  if(cutnames.size()/2!=file.size()) {
    cout<<"Inconsistent entries for file size and cutnames: Suggested cutnameSize = 2* fileSize"<<endl;
    return -1;
  }
  
  std::vector<TEfficiency*> pEff;
  TH1F* demo;

  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    
    TH1F* nosel = (TH1F*) file[filenum]->Get(cutnames[filenum*2]);
    TH1F* sel = (TH1F*) file[filenum]->Get(cutnames[filenum*2+1]);
    
    nosel = (TH1F*) nosel->Rebin(nbins,"newx",rebin);
    sel = (TH1F*) sel->Rebin(nbins,"newx",rebin);
    if(filenum==0) demo = (TH1F*) sel->Clone();
    
    pEff.push_back(0);
    if(TEfficiency::CheckConsistency((*sel),(*nosel))) {
      pEff[filenum] = new TEfficiency((*sel),(*nosel));
    }
    
    sel->SetLineColor(kRed);
    gStyle->SetOptStat(0);
    pEff[filenum]->SetLineWidth(3);
    pEff[filenum]->SetLineColor(coloropt[filenum]);
  }
  
  std::vector<TH1F*> allhists;
  demo->SetTitle("");
  demo->GetXaxis()->SetTitle(xaxistitle);
  demo->GetYaxis()->SetTitle("RECO eff.");
  demo->SetLineColorAlpha(kWhite,1);
  demo->SetFillColorAlpha(kWhite,1);
  allhists.push_back(demo);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){0.7,0.65,0.9,0.95},false,(float []){0.0,1.1},false);
  auto pad = c1->GetPad(3);
  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    pEff[filenum]->Draw("same");
  }
  c1->SaveAs("./dirplots/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
  return -1;
}

int plotter() {

  std::vector<TFile*> file;
  std::vector<TString> cutname;
  std::vector<int> color;  
  std::vector<TString> legend;  

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselgen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 11, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){50,2e5}, (float []){0.65,0.75,0.9,0.95}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #phi");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen M(e, e) / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("ptetaminselgen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 11, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){50,2e5}, (float []){0.65,0.75,0.9,0.95}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #phi");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen M(e, e) / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselsct");
  color.push_back(8);
  legend.push_back("scouting electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "scouting electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "scouting M(e, e) / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselAptetaminselgenelsctbar");
  color.push_back(kBlue);
  legend.push_back("no cut");
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselAptetaminselgenelsctmchbar");
  color.push_back(kRed);
  legend.push_back("|#Delta#eta(sct,gen)|<0.06");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", 1900, 2100, 5, true, false, false, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "#Delta#eta (sct,gen)");
  //comparesamevariable(file, cutname, "qdPhi", 1900, 2200, 5, true, false, false, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "sct e charge #times #Delta#phi (sct,gen)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselAptetaminselgenelsctee");
  color.push_back(kBlue);
  legend.push_back("no cut");
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselAptetaminselgenelsctmchee");
  color.push_back(kRed);
  legend.push_back("|#Delta#eta(sct,gen)|<0.04");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", 1900, 2100, 5, true, false, false, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "#Delta#eta (sct,gen)");
  //comparesamevariable(file, cutname, "qdPhi", 1900, 2200, 5, true, false, false, (float []){1,1e7}, (float []){0.5,0.65,0.75,0.95}, false, "sct e charge #times #Delta#phi (sct,gen)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("ptetaminselgen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselZwindAptetaminselsctmchgenel");
  color.push_back(kRed);
  legend.push_back("sct matched");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 11, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){50,2e5}, (float []){0.65,0.75,0.9,0.95}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #phi");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen M(e, e) / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("narZwindselgen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselZwindAnarZwindselsctmchgenel");
  color.push_back(kRed);
  legend.push_back("sct matched");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 11, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){0.7,2e5}, (float []){0.65,0.75,0.9,0.95}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #phi");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen M(e, e) / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("ptetaminselgen_dielM");
  cutname.push_back("noselZwindAptetaminselsctmchgenel_dielM");
  legend.push_back("efficiency");
  legendEntries = legend;
  double binscutidM[] = {40,60,70,75,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,110,115,120,140};
  //efficiency(file, cutname, 33, binscutidM, "M(e, e) / GeV");
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("ptetaminselgen_eleta");
  cutname.push_back("noselZwindAptetaminselsctmchgenel_eleta");
  legend.push_back("efficiency");
  legendEntries = legend;
  double binseta[] = {-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  //efficiency(file, cutname, 50, binseta, "#eta");
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("ptetaminselgen_elpt");
  cutname.push_back("noselZwindAptetaminselsctmchgenel_elpt");
  legend.push_back("efficiency");
  legendEntries = legend;
  double binspt[] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
  //efficiency(file, cutname, 37, binspt, "p_{T} / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("narZwindselgen_elpt");
  cutname.push_back("noselZwindAnarZwindselsctmchgenel_elpt");
  legend.push_back("efficiency");
  legendEntries = legend;
  //double binspt[] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
  //efficiency(file, cutname, 37, binspt, "p_{T} / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM50file);
  cutname.push_back("noselZwindsctbar");
  color.push_back(8);
  legend.push_back("DY, Z window");
  file.push_back(sigQCDPt20To30file);
  cutname.push_back("noselsctbar");
  color.push_back(2);
  legend.push_back("QCD P_{T} 20To30 ");
  file.push_back(sigQCDPt30To50file);
  cutname.push_back("noselsctbar");
  color.push_back(46);
  legend.push_back("QCD P_{T} 30To50 ");
  coloropt = color;
  legendEntries = legend;
  comparesamevariable(file, cutname, "elsigmaietaieta", -1, -1, 10, true, true, true, (float []){1e-5,1}, (float []){-1,0.65,0.4,0.95}, true);
  //comparesamevariable(file, cutname, "elpt", 10, 80, 1, true, true, true, (float []){0.7,2e5}, (float []){0.65,0.75,0.9,0.95}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.5,0.65,0.75,0.95}, false, "gen electron #phi");
  //comparesamevariable(file, cutname, "dielM", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen M(e, e) / GeV");

  return -1;
}
