import ROOT

#passedsig->Rebin(10)
#(TH1 *) 0x5597551c34d0
#root [3] oldpassedsig->Rebin(10)
#(TH1 *) 0x559754f8c660
#root [4] passedsig->Sumw2()
#root [5] oldpassedsig->Sumw2()
#root [6] allsig->Rebin(10)
#(TH1 *) 0x5597555d0f70
#root [7] oldpassedsig->SetMarkerColor(2)
#root [8] oldpassedsig->SetLineColor(2)
#root [9] oldpassedsig->Divide(allsig)
#(bool) true
#root [10] passedsig->Divide(allsig)
#(bool) true
#root [11] passedsig->Draw()
#Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
#root [12] oldpassedsig->Draw("same")
#root [13] passedsig->SetLineColor(2)
#root [14] passedsig->SetMarkerColor(2)
#root [15] oldpassedsig->SetLineColor(1)
#root [16] oldpassedsig->SetMarkerColor(1)
#root [17] oldpassedsig->Draw()
#Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
#root [18] passedsig->Draw("same")

rebin = 20 #10
bkg = False # True
#true to divide, false for just momentum spectra
eff = True
#fname = "hists_DYToLL_ScoutingSkim220510.root"
fname = "hists_DYToLL_cut16_ScoutingSkim220510.root"
#fname = "hists_QCD_ScoutingSkim220510.root"
#fname = "hists_DoubleElectron_ScoutingSkim220517.root"
#fname = "hists_DoubleElectron_wOld_ScoutingSkim220411.root"
#fname = "hists_DoubleElectron_barrel_ScoutingSkim220411.root"
#fname = "hists_DoubleElectron_endcap_ScoutingSkim220411.root"
f = ROOT.TFile.Open(fname)
if not bkg:
    oldpassed = f.Get("oldpassedsig")
    passed = f.Get("passedsig")
    alll = f.Get("allsig")
else:
    oldpassed = f.Get("oldpassedbkg")
    passed = f.Get("passedbkg")
    alll = f.Get("allbkg")
allgen = f.Get("allgen")
allgen.Rebin(rebin)
oldpassed.Rebin(rebin)
passed.Rebin(rebin)
alll.Rebin(rebin)
if eff:
    alll.Sumw2()
    oldpassed.Sumw2()
    passed.Sumw2()
passed.SetLineColor(2)
passed.SetMarkerColor(2)
oldpassed.SetLineColor(3)
oldpassed.SetMarkerColor(3)
if eff:
    #passed.Divide(alll)
    #oldpassed.Divide(alll)
    alll.Divide(allgen)
    passed.Divide(allgen)
    oldpassed.Divide(allgen)

alll.SetLineWidth(2)
alll.SetLineColor(1)
alll.SetMarkerColor(1)
passed.SetLineWidth(2)
oldpassed.SetLineWidth(2)

if eff:
    alll.GetYaxis().SetTitle("Efficiency")
    alll.GetXaxis().SetTitle("Electron pT (GeV)")
else:
    binwidth = alll.GetXaxis().GetBinWidth(1)
    allgen.GetYaxis().SetTitle("Events / %.2f GeV"%(binwidth))
    allgen.GetXaxis().SetTitle("Gen electron pT (GeV)")
    allgen.SetLineWidth(2)
    allgen.SetLineColor(4)
if not eff:
    allgen.Draw()
    alll.Draw("same")
else:
    alll.Draw()
oldpassed.Draw("same")

passed.Draw("same")
leg = ROOT.TLegend()
sigorbkg = "Background" if bkg else "Signal"
if not eff:
    leg.AddEntry(allgen, "Gen electrons")
leg.AddEntry(alll, "All " + sigorbkg)
#leg.AddEntry(passed, "H < 0.2E + 16")
leg.AddEntry(passed, "H < 0.2E + 16")
leg.AddEntry(oldpassed, "H/E < 0.2")
leg.Draw("same")
raw_input("h")
