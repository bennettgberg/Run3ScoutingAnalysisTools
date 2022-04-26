import ROOT

#fname = "hists_DoubleElectron_ScoutingSkim220404.root"
#fname = "hists_DoubleElectron_ScoutingSkim220411.root"
fname = "hists_new.root"
f = ROOT.TFile(fname)

#here's the 2D hist that we have to fit.
#hname = "hoevsptBkg"
hname = "hVe"

#fnew = ROOT.TFile("hists_new.root", "recreate") 
for fitnum in range(1, 6):

    h2d = f.Get(hname+str(fitnum))
    h2d.Print()
    hnew = h2d

#    #flip the axes
#    hnew = ROOT.TH2F("hnew"+str(fitnum), "H vs. E", 5000, 0, 500, 2000, 0, 200.0)
#    for i in range(h2d.GetNbinsX()):
#        for j in range(h2d.GetNbinsY()):
#            xcent = h2d.GetXaxis().GetBinCenter(i)
#            ycent = h2d.GetYaxis().GetBinCenter(j)
#            content = h2d.GetBinContent(i, j)
#            hnew.Fill(ycent, xcent, content)
#            
#    fnew.cd()
#    hnew.Write()
    #h2d.GetYaxis().SetLimits(0, 400)
    #now do the fit
    #pt doesn't go lower than 2 GeV
    xmin = 0 #2 
    #only interested in low pt (below 10 GeV)
    xmax = 500 #10
    #fitting function
    #f1 = ROOT.TF1("f1", "(gaus(0))", xmin, xmax)
    f1 = ROOT.TF1("f1", "([0]*x + [1])", xmin, xmax)
    #can add option Q for quiet fit, option N to now show.
    #h2d.FitSlicesX(f1, 0, -1, 0, "R")
    #slices = ROOT.TObjArray()
    #print("slices made: " + str(slices))
    #h2d.FitSlicesX(f1, 0, -1, 0, "G5 NR", slices)
    hnew.Fit(f1)
    #param 0 is the const multiplier of the gaussian
    #slices[0].Print()
    #slices[0].Draw("hist")
    hnew.GetXaxis().SetTitle("E (GeV)")
    hnew.GetYaxis().SetTitle("H (GeV)")
    ROOT.gStyle.SetOptFit(1)
    raw_input("h")

#h2d.Draw()
fnew.Close()
