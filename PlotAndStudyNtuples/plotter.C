#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"
#include "RooMsgService.h"

double datasf = 1.09/19; // Data rate - 1.1 Hz, 22 events from parent trigger selection
double sig3cmsf = 1.0/15;
double sig30cmsf = 1.0/28;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
TFile* datahistfile = TFile::Open("hists_data_large.root","READ");
TFile* dym50histfile = TFile::Open("hists_DYToLLM50.root","READ");
//TFile* qcdpt30To50histfile = TFile::Open("hists_QCDPt30To50.root","READ");
TFile* qcdhistfile = TFile::Open("hists_QCD.root","READ");
TFile* tempfile = TFile::Open("hists_temp.root","RECREATE");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
std::vector<double> scale{-1, 1};

int autoplotter(std::vector<TFile*> file, std::vector<TString> cutname) {

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_autoplotted";

  TList* listOfKeys = file[0]->GetListOfKeys(); // Get the list of histograms
  for(unsigned int ctr=0; ctr<listOfKeys->GetEntries(); ctr++) {
    TString keyName = listOfKeys->At(ctr)->GetName(); // Get the histogram name
    bool isWithCut = keyName.BeginsWith(cutname[0]); // Check if the histogram begins with a certain cut

    if(!isWithCut) continue;

    int maxBins = 10000;
    TH1F* histVarOrig = (TH1F*) file[0]->Get(keyName);
    TH1F* histVar = (TH1F*) histVarOrig->Clone(histVarOrig->GetName());
    histVar->Sumw2();
    //if(histVar->GetNbinsX()>maxBins) histVar->Rebin(maxBins);
    /*
    // Go over the histogram and rebin it based on the contents
    std::vector<double> binlowedge;
    binlowedge.push_back(histVar->GetBinLowEdge(1));
    double maxvalue = histVar->GetMaximum();
    double binval = 0.0;
    double cutFracWMax = 0.05;
    for(unsigned int binctr=1; binctr<histVar->GetNbinsX(); binctr++) {
      binval += histVar->GetBinContent(binctr);
      if(binval>cutFracWMax*maxvalue) {
	binlowedge.push_back(histVar->GetBinLowEdge(binctr+1));
	binval = 0.0;
      }
    }
    if(binlowedge[binlowedge.size()-1]!=histVar->GetBinLowEdge(histVar->GetNbinsX()+1)) {
      binlowedge.push_back(histVar->GetBinLowEdge(histVar->GetNbinsX()+1));
    }

    histVar = histVar->Rebin(binlowedge.size()-1,histVar->GetTitle(),binlowedge)
    */
    TCanvas* c1;
    c1 = new TCanvas();
    gStyle->SetOptStat(0);
    //c1->SetLogx(true);
    histVar->Draw();
    c1->SaveAs("./dirplots/"+foldername+"/"+keyName+".png");
    
  }

  return -1;
}


int invmee_specialplot(TString selection, double res) {

  if(res>1) {
    cout<<"Resolution input as fraction. Cannot be greater than one. Invalid input to invmee_specialplot()."<<endl;
    return -1;
  }
  
  auto invmeehistorig = (TH1F*) datahistfile->Get(selection);
  TH1F* invmeehist = (TH1F*) invmeehistorig->Clone(invmeehistorig->GetName());
  
  vector<double> xbins;
  xbins.push_back(0.5);
  while(xbins[xbins.size()-1]<100) {
    double val = ceil(100.0*(1+res)*xbins[xbins.size()-1])/100.0;
    xbins.push_back(val);
  }

  invmeehist = (TH1F*) invmeehist->Rebin(xbins.size()-1, selection, &xbins[0]);
  invmeehist->SetTitle("");
  
  TCanvas* c1;
  c1 = enhance_plotter({invmeehist}, {"Run 3 Scouting"}, "M(e,e)", "number of events", (float []){0.15,0.7,0.4,0.9}, true, false, (float []){0,static_cast<float>(1.5*invmeehist->GetMaximum())}, false, {"hist e1"}, {20}, {2}, {"lep"});
  c1->SaveAs(selection+".png");
  
  return -1;
}

int fitinvmee(TString selection) {

  auto invmeehistorig = (TH1F*) datahistfile->Get(selection);
  TH1F* invmeehist = (TH1F*) invmeehistorig->Clone(invmeehistorig->GetName());

  invmeehist->SetTitle("");
  invmeehist->Rebin(100);
  invmeehist->GetXaxis()->SetRange(71,150);

  //TF1* fitbkg = new TF1("fitbkg","exp([0]+[1]*x)",40,65);
  TF1* fitbkg = new TF1("fitbkg","exp([0]+[1]*x)",60,65);
  fitbkg->SetParameters(5,-0.002);
  invmeehist->Fit("fitbkg","R");
  fitbkg->SetLineColor(kBlue);
  double b0 = fitbkg->GetParameter(0);
  double be0 = fitbkg->GetParError(0);
  double b1 = fitbkg->GetParameter(1);
  double be1 = fitbkg->GetParError(1);

  TF1* fitpeak = new TF1("fitpeak","[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",80,100);
  fitpeak->SetParameters(702,91,5);
  invmeehist->Fit("fitpeak","R");
  fitpeak->SetLineColor(kGreen+2);
  double p0 = fitpeak->GetParameter(0);
  double pe0 = fitpeak->GetParError(0);
  double p1 = fitpeak->GetParameter(1);
  double pe1 = fitpeak->GetParError(1);
  double p2 = fitpeak->GetParameter(2);
  double pe2 = fitpeak->GetParError(2);

  TF1* fitfunc = new TF1("fitZmass","exp([0]+[1]*x)+[2]*exp(-0.5*((x-[3])/[4])*((x-[3])/[4]))",60,140);
  fitfunc->SetParameters(b0, b1, p0, p1, p2);
  fitfunc->SetParLimits(0, b0-10*be0, b0+10*be0);
  fitfunc->SetParLimits(1, b1-10*be1, b1+10*be1);
  fitfunc->SetParLimits(2, p0-10*pe0, p0+10*pe0);
  fitfunc->SetParLimits(3, p1-10*pe1, p1+10*pe1);
  fitfunc->SetParLimits(4, p2-10*pe2, p2+10*pe2);
  
  TFitResultPtr fitres = invmeehist->Fit("fitZmass","RS");

  cout<<"Chi2/Ndf = "<<fitres->Chi2()<<"/"<<fitres->Ndf()<<endl;

  TF1* Zpeak = new TF1("Zpeak","[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",60,140);
  Zpeak->SetParameters(fitfunc->GetParameter(2), fitfunc->GetParameter(3), fitfunc->GetParameter(4));
  cout<<"Signal integral in the (83.73,92.53): "<<Zpeak->Integral(83.73,92.53)<<endl;
  cout<<"Bkg integral in the (83.73,92.53): "<<fitfunc->Integral(83.73,92.53)-Zpeak->Integral(83.73,92.53)<<endl;
  cout<<"Signal integral in the (80,100): "<<Zpeak->Integral(80,100)<<endl;
  cout<<"Bkg integral in the (80,100): "<<fitfunc->Integral(80,100)-Zpeak->Integral(80,100)<<endl;
  cout<<"Signal integral in the (70,100): "<<Zpeak->Integral(70,100)<<endl;
  cout<<"Bkg integral in the (70,100): "<<fitfunc->Integral(70,100)-Zpeak->Integral(70,100)<<endl;
  cout<<"Signal integral in the (50,110): "<<Zpeak->Integral(50,110)<<endl;
  cout<<"Bkg integral in the (50,110): "<<fitfunc->Integral(50,110)-Zpeak->Integral(50,110)<<endl;
  cout<<"Signal integral in the (45,115): "<<Zpeak->Integral(45,115)<<endl;
  cout<<"Bkg integral in the (45,115): "<<fitfunc->Integral(45,115)-Zpeak->Integral(45,115)<<endl;
  cout<<"Signal integral in the (40,120): "<<Zpeak->Integral(40,120)<<endl;
  cout<<"Bkg integral in the (40,120): "<<fitfunc->Integral(40,120)-Zpeak->Integral(40,120)<<endl;

  // Make signal significance plot
  double sigmafrommean[1000], signalsignificance[1000];
  int n = 1000, maxsignificancen=-1;
  for(unsigned int xctr = 0; xctr<n; xctr++) {
    double xval = 0.01+xctr*0.01;
    sigmafrommean[xctr] = xval;
    double s = Zpeak->Integral(fitfunc->GetParameter(3)-xval*fitfunc->GetParameter(4),fitfunc->GetParameter(3)+xval*fitfunc->GetParameter(4));
    double spb = fitfunc->Integral(fitfunc->GetParameter(3)-xval*fitfunc->GetParameter(4),fitfunc->GetParameter(3)+xval*fitfunc->GetParameter(4));
    double b = spb-s;
    signalsignificance[xctr] = (s/sqrt(b));
    if(xctr!=0) {
      if(signalsignificance[xctr]>signalsignificance[xctr-1]) {
	maxsignificancen = xctr;
      }
    }
  }
  TGraph *signalsig_graph = new TGraph(n, sigmafrommean, signalsignificance);
  cout<<"Maximum significance in the range ["<<(fitfunc->GetParameter(3)-sigmafrommean[maxsignificancen]*fitfunc->GetParameter(4))<<", "<<(fitfunc->GetParameter(3)+sigmafrommean[maxsignificancen]*fitfunc->GetParameter(4))<<"] with significance - "<<signalsignificance[maxsignificancen]<<endl;
  
  fitfunc->SetLineWidth(3);
  
  TCanvas* c1;
  //invmeehist->Draw();
  //fitbkg->Draw("SAME");
  //fitpeak->Draw("SAME");
  c1 = enhance_plotter({invmeehist}, {"Scouting"}, "M(e,e)", "number of events", (float []){0.15,0.7,0.4,0.9}, false, false, (float []){0,/*static_cast<float>(1.5*invmeehist->GetMaximum())*/2750}, false, {"hist"}, {20}, {2}, {"lep"});
  fitfunc->Draw("same");
  c1->SaveAs(selection+"_invmeefit.png");

  TCanvas* c2 = new TCanvas();
  signalsig_graph->Draw("AC*");
  c2->SaveAs(selection+"_signalsignificance.png");

  return -1;
}

int fitinvmee_roofit(TString selection, double fitrange[2], double *SpBinSpB, double *BinSpB, double SpBrange[2]) {
  using namespace RooFit;

  auto invmeehistorig = (TH1F*) datahistfile->Get(selection);
  TH1F* invmeehist = (TH1F*) invmeehistorig->Clone(invmeehistorig->GetName());
  invmeehist->SetTitle("");
  invmeehist->Rebin(100);
  invmeehist->GetXaxis()->SetRange(30,150);

  invmeehist->SetDirectory(0);

  // Set RooFit message service level
  RooMsgService::instance().setStreamStatus(1,false);
  RooMsgService::instance().setStreamStatus(0,false);
  
  // Declare observable x
  double invm_min = 0.0, invm_max=150.0;
  RooRealVar x("x", "x", invm_min, invm_max);
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable
  RooDataHist dh("dh", "dh", x, Import(*invmeehist));
  
  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot *frame1 = x.frame(Title("Signal(CB)+Bkg(exp)"));
  //RooPlot *frame2 = x.frame(Title("CrystalBall"));
  dh.plotOn(frame1, Name("data"), DataError(RooAbsData::SumW2));
  //dh.plotOn(frame2, Name("data"), DataError(RooAbsData::SumW2));
  /*
  // Build a Gaussian pdf
  RooGaussian Gaus("gaus", "Gaussian", x, meanG, sigma);
  
  // Build a skewed Gaussian pdf 
  RooRealVar meanSKG("meanSKG", "meanSKG", 91, 81, 101);
  RooRealVar sigmaSKG("sigmaSKG", "sigmaSKG", 10, 1, 21);
  RooRealVar tail("tail", "tail", 0, -4, 4);
  RooNovosibirsk SkGaus("SkGaus", "SkewedGaussian", x, meanSKG, sigmaSKG, tail);
  */
  
  // Build a crystall ball pdf
  RooRealVar meanG("meanG", "meanG", -2, -10, -0.1);
  RooRealVar sigma("sigma", "sigma", 2.34, 1, 11);
  RooRealVar a("a", "alpha", 1, 0.1, 10);
  RooRealVar n("n", "n", 1, 0, 20);
  RooCBShape CBShape("CBShape", "Crystal_Ball", x, meanG, sigma, a, n);

  // Build a breit-wigner pdf
  RooRealVar mean("mean", "mean", 91, 80, 100);
  RooRealVar width("width", "width", 5, 1, 10);
  RooBreitWigner bwMassPeak("bwMassPeak", "BreitWigner", x, mean, width);

  // Build a Voigtian pdf
  //RooVoigtian voigt("voigt", "voigt", x, mean, width, sigma);
  
  // Build convolution for signal
  x.setBins(10000, "cache");
  //RooFFTConvPdf sig("sig", "BW(x)G", x, bwMassPeak, Gaus);
  RooFFTConvPdf sig("SigFitModel", "BW(x)CB", x, bwMassPeak, CBShape);
  
  // Build Exponential pdf
  RooRealVar decay("decay", "decay", -0.05, -0.1, -1e-4);
  RooExponential bkg("BkgFitModel", "Exponential", x, decay);
  
  // Build composite pdf
  RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.19, 0.0, 1.0);
  //RooAddPdf model("model", "sig+bkg", RooArgList(bkg, voigt), bkgfrac);
  RooAddPdf model("SigPBkgFitModel", "sig+bkg", RooArgList(bkg, sig), bkgfrac);
  
  x.setRange("R1",fitrange[0],fitrange[1]) ;

  //RooCmdArg arg = RooCmdArg::PrintLevel(-1);
  model.fitTo(dh,Range("R1"), PrintLevel(-1));
  model.plotOn(frame1, Name("Background"), Components(bkg), LineStyle(kDashed), LineColor(kRed));
  //model.plotOn(frame1, Name("Signal"), Components(voigt), LineStyle(10), LineColor(kBlack));
  model.plotOn(frame1, Name("Signal"), Components(sig), LineStyle(10), LineColor(kBlack));
  model.plotOn(frame1, Name("SigPBkgFitModel"), Components("SigPBkgFitModel"), LineColor(kBlue));

  x.setRange("SpBRange", SpBrange[0], SpBrange[1]);
  RooAbsReal* spbcnt = model.createIntegral(x, NormSet(x), Range("SpBRange"));
  RooAbsReal* bcnt = bkg.createIntegral(x, NormSet(x), Range("SpBRange"));
  double evtcnt = invmeehist->Integral(invmeehist->FindBin(invm_min), invmeehist->FindBin(invm_max-0.0001));
  *SpBinSpB = spbcnt->getValV()*evtcnt;
  *BinSpB = bcnt->getValV()*bkgfrac.getValV()*evtcnt;

  //tempfile->WriteObject(&model,selection+"_"+model.GetName());

  // Optional plotting of the invariant mass distribution
  // And the signal significance plot when required. 
  TCanvas* c1;
  c1 = enhance_plotter({invmeehist}, {"Scouting"}, "M(e,e)", "number of events", (float []){0.15,0.7,0.4,0.9}, false, false, (float []){0,1500}, false, {"hist"}, {20}, {2}, {"lep"});
  frame1->Draw();
  c1->SaveAs(selection+"_roofit1_invmeefit.png");
  /*  
  // Make signal significance plot
  double sigmafrommean[1000], signalsignificance[1000];
  int maxcnt = 1000, maxsignificance=-1;
  for(unsigned int xctr = 0; xctr<maxcnt; xctr++) {
    double xval = 0.01+xctr*0.01;
    sigmafrommean[xctr] = xval;
    x.setRange("xvalsigma",mean.getValV()-xval*sigma.getValV(), mean.getValV()+xval*sigma.getValV());
    RooAbsReal* scnt = sig.createIntegral(x, NormSet(x), Range("xvalsigma"));
    RooAbsReal* bcnt = bkg.createIntegral(x, NormSet(x), Range("xvalsigma"));
    signalsignificance[xctr] = (scnt->getValV()/sqrt(bcnt->getValV()));
    if(xctr!=0) {
      if(signalsignificance[xctr]>signalsignificance[xctr-1]) {
	maxsignificance = xctr;
      }
    }
  }
  TGraph *signalsig_graph = new TGraph(maxcnt, sigmafrommean, signalsignificance);
  cout<<"Maximum significance in the range ["<<mean.getValV()-sigmafrommean[maxsignificance]*sigma.getValV()<<", "<<mean.getValV()+sigmafrommean[maxsignificance]*sigma.getValV()<<"] with significance - "<<signalsignificance[maxsignificance]<<endl;

  TCanvas* c2 = new TCanvas();
  signalsig_graph->Draw("AC*");
  c2->SaveAs(selection+"_roofit1_signalsignificance.png");
  */

  return -1;
}

int subtractsideband(TFile* histfile, TString SpBhistname, double scaleSpB, TString Bhistname, double scaleB, int xbinlow, int xbinhigh, int rebin, TString xaxistitle, TString yaxistitle, bool logY, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}) {

  TString foldername = "";
  foldername += ((TString)histfile->GetName()).ReplaceAll(".root","_noB").ReplaceAll("hists_","");

  TH1F* temphistholder = (TH1F*) histfile->Get(SpBhistname);
  TH1F* hist_SpB = (TH1F*)temphistholder->Clone(temphistholder->GetName());
  temphistholder = (TH1F*) histfile->Get(Bhistname);
  TH1F* hist_B = (TH1F*)temphistholder->Clone(temphistholder->GetName());

  cout<<"Normalisation factor for SpB histo. (should be close to 1): "<<scaleSpB/hist_SpB->Integral()<<endl;
  hist_SpB->Scale(scaleSpB/hist_SpB->Integral());
  hist_B->Scale(scaleB/hist_B->Integral());
  TH1F* hist_S = (TH1F*) hist_SpB->Clone(SpBhistname.ReplaceAll("Zwind","Zsignal"));
  hist_S->Add(hist_B, -1.0);

  tempfile->WriteObject(hist_SpB,hist_S->GetName());  

  coloropt.clear();
  legendEntries.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  
  std::vector<TH1F*> allhists;

  TString SpBCloneName = hist_SpB->GetName();
  SpBCloneName.Append("_clone");
  allhists.push_back((TH1F*) hist_SpB->Clone(SpBCloneName));
  coloropt.push_back(kMagenta);
  legendEntries.push_back("SpB");
  histtype.push_back("hist e1");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");
  scale.push_back(1);
  
  TString BCloneName = hist_B->GetName();
  BCloneName.Append("_clone");
  allhists.push_back((TH1F*) hist_B->Clone(BCloneName));
  coloropt.push_back(kRed);
  legendEntries.push_back("B");
  histtype.push_back("same hist e1");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");
  scale.push_back(1);
  
  TString SCloneName = hist_S->GetName();
  SCloneName.Append("_clone");
  allhists.push_back((TH1F*) hist_S->Clone(SCloneName));
  coloropt.push_back(kBlue);
  legendEntries.push_back("S");
  histtype.push_back("same hist e1");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");
  scale.push_back(1);
  
  if(rebin==-1) rebin = 1;
  //xbinlow = xbinlow==-1?(underflow?0:1):xbinlow;
  //xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow==-1?1:xbinlow;
  xbinhigh = xbinhigh==-1?allhists[0]->GetNbinsX():xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;
  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);
    allhists[histctr]->SetTitle("");
    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
    allhists[histctr]->SetMarkerColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, xaxistitle, yaxistitle, legPos, false, logY, yrange, false, histtype, markerstyle, markersize, legendmarkerstyle);
  c1->SaveAs("./dirplots/"+foldername+"/"+hist_S->GetName()+".png");

  return -1;
}

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
    TH1F* temphistholder = (TH1F*) file[histctr]->Get(cutwithvar);
    allhists.push_back((TH1F*)temphistholder->Clone(temphistholder->GetName()));
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

    if(scale[0]!=-1) allhists[histctr]->Scale(scale[histctr]);
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
    allhists[histctr]->SetMarkerColor(coloropt[histctr]);
    if(markerstyle[histctr]==0) allhists[histctr]->SetFillColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos, false, logY,yrange,normalize,histtype,markerstyle,markersize,legendmarkerstyle);
  if(!normalize) c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");
  else  c1->SaveAs("./dirplots/"+foldername+"/"+var+"_normed.png");

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
    
    TH1F* noselorig = (TH1F*) file[filenum]->Get(cutnames[filenum*2]);
    TH1F* nosel = (TH1F*) noselorig->Clone(noselorig->GetName());
    TH1F* selorig = (TH1F*) file[filenum]->Get(cutnames[filenum*2+1]);
    TH1F* sel = (TH1F*) selorig->Clone(selorig->GetName());
    
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
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){0.7,0.65,0.9,0.95},false,false,(float []){0.0,1.1},false);
  auto pad = c1->GetPad(3);
  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    pEff[filenum]->Draw("same");
  }
  c1->SaveAs("./dirplots/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
  return -1;
}

void group_plotter(std::vector<TFile*> file, std::vector<TString> cutname, bool isBar) {

  //all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  //all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  //all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  //all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  comparesamevariable(file, cutname, "elpt", 5, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron p_{T} [GeV]");
  comparesamevariable(file, cutname, "eld0", 8500, 11500, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron d_{0} [cm]");
  comparesamevariable(file, cutname, "eldetain", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta(SC seed, track)");
  if(isBar) {
    comparesamevariable(file, cutname, "eldphiin", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#phi(SC, track)");
    comparesamevariable(file, cutname, "elsigmaietaieta", -1, 400, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #sigmai#etai#eta");
    comparesamevariable(file, cutname, "elhoe", -1, -1, 50, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron H/E");
    comparesamevariable(file, cutname, "elooemoop", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron E^{-1}-p^{-1} [GeV^{-1}]");
    comparesamevariable(file, cutname, "elmhits", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron missing hits");
    comparesamevariable(file, cutname, "elecaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron ECAL iso. [GeV]");
    comparesamevariable(file, cutname, "elhcaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron HCAL iso. [GeV]");
    comparesamevariable(file, cutname, "eltkiso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron track iso. [GeV]");
    comparesamevariable(file, cutname, "elr9", -1, 150, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron r9");
    comparesamevariable(file, cutname, "elsmin", -1, 1100, 5, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smin");
    comparesamevariable(file, cutname, "elsmaj", -1, 700, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smaj");
  }
  else {
    comparesamevariable(file, cutname, "eldetain", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta(SC seed, track)");
    comparesamevariable(file, cutname, "eldphiin", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#phi(SC, track)");
    comparesamevariable(file, cutname, "elsigmaietaieta", -1, 800, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #sigmai#etai#eta");
    comparesamevariable(file, cutname, "elhoe", -1, -1, 50, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron H/E");
    comparesamevariable(file, cutname, "elooemoop", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron E^{-1}-p^{-1} [GeV^{-1}]");
    comparesamevariable(file, cutname, "elmhits", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron missing hits");
    comparesamevariable(file, cutname, "elecaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron ECAL iso. [GeV]");
    comparesamevariable(file, cutname, "elhcaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron HCAL iso. [GeV]");
    comparesamevariable(file, cutname, "eltkiso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron track iso. [GeV]");
    comparesamevariable(file, cutname, "elr9", -1, 300, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron r9");
    comparesamevariable(file, cutname, "elsmin", -1, 2200, 5, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smin");
    comparesamevariable(file, cutname, "elsmaj", -1, 1400, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smaj");
  }
    
}

int plotter() {

  std::vector<TFile*> file;
  std::vector<TString> cutname;
  std::vector<TString> legend;  

  //////////////////////// DELETE EVERYTHING BEFORE THIS \\\\\\\\\\\\\\\\\\\\\\\\\ \

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("noselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 20, 1, true, true, true, (float []){8e-1,2e8}, (float []){0.55,0.7,0.75,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", -1, 400, 1, true, true, true, (float []){8e-1,2e9}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 230, 770, 1, true, false, false, (float []){8e-1,2e8}, (float []){0.55,0.7,0.75,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, false, true, true, (float []){0,5e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron #phi");
  //comparesamevariable(file, cutname, "dielM", 500, 20000, 100, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");
  //comparesamevariable(file, cutname, "leadsublead_dielM", 3500, 20000, 100, false, true, true, (float []){8e-1,1e4}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");
  //invmee_specialplot("noselsct_dielM", 0.02);

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("leadepemselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Atleast (e_{12}^{+}, e_{12}^{-})");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("leadepem_only_selsct");
  coloropt.push_back(kRed);
  legend.push_back("(e_{12}^{+}, e_{12}^{-}) only");
  histtype.push_back("same p e1");
  markerstyle.push_back(20);
  markersize.push_back(1);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("leadepem_Zwind_selsct");
  coloropt.push_back(kGreen+2);
  legend.push_back("80<M(e_{12}^{+}, e_{12}^{-})<100");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "leadsublead_dielM", 3500, 20000, 100, false, true, true, (float []){8e-1,1e4}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("noselsct");
  coloropt.push_back(kBlack);
  legend.push_back("no sel. bkg. rich");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("leadepem_Zwind_selsct");
  coloropt.push_back(kGreen+2);
  legend.push_back("80<M(e_{12}^{+}, e_{12}^{-})<100");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;

  //comparesamevariable(file, cutname, "elmult", 5, 20, 1, true, true, true, (float []){2e-8,2}, (float []){0.55,0.7,0.75,0.95}, true, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", -1, 150, 1, true, true, true, (float []){2e-4,2}, (float []){0.55,0.7,0.75,0.95}, true, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 230, 770, 1, true, false, false, (float []){2e-8,2}, (float []){0.55,0.7,0.75,0.95}, true, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){5e-5,2}, (float []){0.55,0.7,0.75,0.95}, true, "electron #phi");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  /*
  file.push_back(datahistfile);
  cutname.push_back("leadepem_Zwindptgt20_selsct");
  coloropt.push_back(kBlack);
  legend.push_back("Z wind., p_{T}>20");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);
  
  file.push_back(dym50histfile);
  cutname.push_back("leadepem_ptgt20_selsct");
  coloropt.push_back(kRed);
  legend.push_back("(Z#rightarrowee, p_{T}>20");
  histtype.push_back("same hist");
  markerstyle.push_back(1);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");
  scale.push_back(1);
  */
  file.push_back(datahistfile);
  cutname.push_back("leadepem_ptgt20_selsct");
  coloropt.push_back(kGreen+2);
  legend.push_back("(e_{12}^{+}, e_{12}^{-}), p_{T}>20");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");
  scale.push_back(1);

  legendEntries = legend;

  //comparesamevariable(file, cutname, "leadsublead_dielM", 3500, 20000, 100, false, false, false, (float []){8e-1,1.5e3}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");
  //fitinvmee("tightselsct_leadbarsubleadbar_dielM");
  //fitinvmee("tightselsct_leadbarsubleadbar_dielM");
  //fitinvmee_roofit("tightselsct_leadbarsubleadbar_dielM");
  //fitinvmee("tightselsct_leadecsubleadec_dielM");
  //fitinvmee_roofit("tightselsct_leadecsubleadec_dielM");

  cutname.clear();
  scale.clear();
  scale.push_back(1);
  scale.push_back((2740.0/(2740.0+8365.0))*(15509.0/5147603.0));
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elsigmaietaieta","noselsctbar_elsigmaietaieta"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elhoe","noselsctbar_elhoe"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_eltkiso","noselsctbar_eltkiso"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elsmin","noselsctbar_elsmin"});

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  scale.push_back(-1);

  file.push_back(dym50histfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctbar");
  coloropt.push_back(kRed-6);
  legend.push_back("Z#rightarrowee, p_{T}>20");
  histtype.push_back("same hist");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");

  file.push_back(tempfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctbar");
  coloropt.push_back(kRed);
  legend.push_back("data, Z, p_{T}>20");
  histtype.push_back("same p");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("noselsctbar");
  coloropt.push_back(kGreen+2);
  legend.push_back("no selection, bkg.");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;

  //comparesamevariable(file, cutname, "elsigmaietaieta", -1, 250, 4, true, true, true, (float []){5e-4,1}, (float []){0.6,0.7,0.85,0.95}, true, "electron #sigmai#etai#eta");
  //comparesamevariable(file, cutname, "elhoe", -1, 20000, 200, true, true, true, (float []){5e-4,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron H/E");
  //comparesamevariable(file, cutname, "eltkiso", -1, -1, 200, true, true, true, (float []){5e-4,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron track iso. [GeV]");
  //comparesamevariable(file, cutname, "elsmin", -1, 1500, 5, true, true, true, (float []){5e-4,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron smin");

  cutname.clear();
  scale.clear();
  scale.push_back(1);
  scale.push_back((2740.0/(2740.0+8365.0))*(6775.0/5956977.0));
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elsigmaietaieta","noselsctec_elsigmaietaieta"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elhoe","noselsctec_elhoe"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_eltkiso","noselsctec_eltkiso"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elsmin","noselsctec_elsmin"});

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  scale.push_back(-1);

  file.push_back(dym50histfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctec");
  coloropt.push_back(kRed-6);
  legend.push_back("Z#rightarrowee, p_{T}>20");
  histtype.push_back("same hist");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");

  file.push_back(tempfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctec");
  coloropt.push_back(kRed);
  legend.push_back("data, Z, p_{T}>20");
  histtype.push_back("same p");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("noselsctec");
  coloropt.push_back(kGreen+2);
  legend.push_back("no selection, bkg.");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;

  //comparesamevariable(file, cutname, "elsigmaietaieta", 100, 600, 4, true, true, true, (float []){1e-3,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron #sigmai#etai#eta");
  //comparesamevariable(file, cutname, "elhoe", -1, 15000, 200, true, true, true, (float []){1e-3,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron H/E");
  //comparesamevariable(file, cutname, "eltkiso", -1, -1, 200, true, true, true, (float []){1e-3,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron track iso. [GeV]");
  //comparesamevariable(file, cutname, "elsmin", 100, 1000, 10, true, true, true, (float []){1e-3,8e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron smin");

  cutname.clear();
  scale.clear();
  scale.push_back(1);
  scale.push_back((2740.0/(2740.0+8365.0))*(15509.0/5147603.0));
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elsigmaietaieta","noselsctbar_elsigmaietaieta"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elhoe","noselsctbar_elhoe"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_eltkiso","noselsctbar_eltkiso"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctbar_elsmin","noselsctbar_elsmin"});

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  scale.push_back(-1);

  file.push_back(qcdhistfile);
  cutname.push_back("noselsctbar");
  coloropt.push_back(kGreen-10);
  legend.push_back("QCD MC");
  histtype.push_back("same hist");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");

  file.push_back(tempfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctbar");
  coloropt.push_back(kRed);
  legend.push_back("data, Z, p_{T}>20");
  histtype.push_back("same p");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("noselsctbar");
  coloropt.push_back(kGreen+2);
  legend.push_back("no selection, bkg.");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;

  //comparesamevariable(file, cutname, "elsigmaietaieta", -1, 250, 4, true, true, true, (float []){5e-4,1}, (float []){0.6,0.7,0.85,0.95}, true, "electron #sigmai#etai#eta");
  //comparesamevariable(file, cutname, "elhoe", -1, 20000, 200, true, true, true, (float []){5e-4,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron H/E");
  //comparesamevariable(file, cutname, "eltkiso", -1, -1, 200, true, true, true, (float []){5e-4,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron track iso. [GeV]");
  //comparesamevariable(file, cutname, "elsmin", -1, 1500, 5, true, true, true, (float []){5e-4,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron smin");

  cutname.clear();
  scale.clear();
  scale.push_back(1);
  scale.push_back((2740.0/(2740.0+8365.0))*(6775.0/5956977.0));
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elsigmaietaieta","noselsctec_elsigmaietaieta"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elhoe","noselsctec_elhoe"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_eltkiso","noselsctec_eltkiso"});
  //subtractsideband(datahistfile, {"leadepem_Zwindptgt20_selsctec_elsmin","noselsctec_elsmin"});

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();
  scale.clear();
  scale.push_back(-1);

  file.push_back(qcdhistfile);
  cutname.push_back("noselsctec");
  coloropt.push_back(kGreen-10);
  legend.push_back("QCD MC");
  histtype.push_back("same hist");
  markerstyle.push_back(0);
  markersize.push_back(0);
  legendmarkerstyle.push_back("f");

  file.push_back(tempfile);
  cutname.push_back("leadepem_Zwindptgt20_selsctec");
  coloropt.push_back(kRed);
  legend.push_back("data, Z, p_{T}>20");
  histtype.push_back("same p");
  markerstyle.push_back(24);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  file.push_back(datahistfile);
  cutname.push_back("noselsctec");
  coloropt.push_back(kGreen+2);
  legend.push_back("no selection, bkg.");
  histtype.push_back("same p");
  markerstyle.push_back(40);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;

  //comparesamevariable(file, cutname, "elsigmaietaieta", 100, 600, 4, true, true, true, (float []){1e-3,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron #sigmai#etai#eta");
  //comparesamevariable(file, cutname, "elhoe", -1, 15000, 200, true, true, true, (float []){1e-3,2e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron H/E");
  //comparesamevariable(file, cutname, "eltkiso", -1, -1, 200, true, true, true, (float []){1e-3,2}, (float []){0.6,0.7,0.85,0.95}, true, "electron track iso. [GeV]");
  //comparesamevariable(file, cutname, "elsmin", 100, 1000, 10, true, true, true, (float []){1e-3,8e-1}, (float []){0.6,0.7,0.85,0.95}, true, "electron smin");

  //invmee_specialplot("noselsct_leadsublead_dielM", 0.02);
  //invmee_specialplot("vetoselsct_leadsublead_dielM", 0.02);
  //invmee_specialplot("looseselsct_leadsublead_dielM", 0.02);
  //invmee_specialplot("mediumselsct_leadsublead_dielM", 0.02);
  //invmee_specialplot("tightselsct_leadsublead_dielM", 0.02);
  
  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("looseselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "leadsublead_dielM", 1250, 1470, 5, false, true, true, (float []){0,400}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("mediumselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "leadsublead_dielM", 1250, 1470, 5, false, true, true, (float []){0,400}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "leadsublead_dielM", 1250, 1470, 5, false, true, true, (float []){0,400}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightsel_Zwind_sct");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 20, 1, true, true, true, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", -1, 150, 1, true, true, true, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 230, 770, 1, true, false, false, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron #phi");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightsel");
  coloropt.push_back(kBlack);
  legend.push_back("Run2022B Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "Zwind_sctbar_elpt", -1, 150, 1, true, true, true, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "SideBand_sctbar_elpt", -1, 150, 1, true, true, true, (float []){1e-1,1e5}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");

  double SpBinSpB = 0.0, BinSpB = 0.0;
  fitinvmee_roofit("tightselsct_leadbarsubleadbar_dielM", (double []){50.0, 140.0}, &SpBinSpB, &BinSpB, (double []){86.0, 97.0});
  cout<<"Integral counts: "<<SpBinSpB<<"\t"<<BinSpB<<endl;
  subtractsideband(datahistfile, "tightsel_Zwind_sctbar_elpt", SpBinSpB, "tightsel_SideBand_sctbar_elpt", BinSpB, -1, 150, 1, "electron p_{T} [GeV]", "number of events", false, (float []){-1e2,1e3}, (float []){0.55,0.7,0.75,0.95});
  fitinvmee_roofit("tightselsct_leadecsubleadec_dielM", (double []){50.0, 140.0}, &SpBinSpB, &BinSpB, (double []){79.0, 98.0});
  cout<<"Integral counts: "<<SpBinSpB<<"\t"<<BinSpB<<endl;
  subtractsideband(datahistfile, "tightsel_Zwind_sctec_elpt", SpBinSpB, "tightsel_SideBand_sctec_elpt", BinSpB, -1, 150, 1, "electron p_{T} [GeV]", "number of events", false, (float []){-5e2,1e3}, (float []){0.55,0.7,0.75,0.95});
  
  tempfile->Close();
  return -1;
}
