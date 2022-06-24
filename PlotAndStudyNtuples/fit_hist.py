import ROOT

#can only read fast if hists_new.root exists and is filled in.
read_fast = False
rebin = 20 #what factor to rebin by
#fname = "hists_DoubleElectron_ScoutingSkim220404.root"
#fname = "hists_DoubleElectron_ScoutingSkim220411.root"
#fname = "hists_DoubleElectron_ScoutingSkim220506.root"
fname = "hists_DYToLL_ScoutingSkim220510.root"
#fname = "hists_QCD_ScoutingSkim220510.root"
#fname = "hists_new.root"
f = ROOT.TFile(fname)
f.Print()
#here's the 2D hist that we have to fit.
#hname = "hoevsptBkg"
hname = "hVe"
if read_fast:
    hname = "hnew"

if not read_fast:
    fnew = ROOT.TFile("hists_new.root", "recreate") 

#for fitnum in range(0, 7):
for fitnum in range(7, 8):

    h2d = f.Get(hname+str(fitnum))
    h2d.Print()
    #h2d.Draw()
    if read_fast:
        hnew = h2d
    #make 1 d histogram, x bin vs. mean of all y values in that bin
    h2d.RebinX(rebin)
    nbins = h2d.GetXaxis().GetNbins()
    #new_nbins = nbins / rebin
    if not read_fast:
        nybins = h2d.GetYaxis().GetNbins()
    xmin = h2d.GetXaxis().GetXmin()
    xmax = h2d.GetXaxis().GetXmax()
    print("xmin: %f, xmax: %f, nbins: %d "%(xmin, xmax, nbins))
    if not read_fast:
        #hnew = ROOT.TH1F("hnew"+str(fitnum), "hnew"+str(fitnum), new_nbins, xmin, xmax)
        hnew = ROOT.TH1F("hnew"+str(fitnum), "hnew"+str(fitnum), nbins, xmin, xmax)
    #for each bin, get the mean y value and std dev y value
    #  then set the point to be mean +/- stdev 
    if not read_fast:
        binmean = 0.
        binrms = 0.
        nentries = 0
        h2d.Draw()
        #raw_input("drawing h2d.")
        for xbin in range(nbins):
            #if xbin%rebin == 1:
            if xbin > 20: break
            binmean = 0.
            binrms = 0.
            nentries = 0
            #h2d.GetXaxis().SetRange(xbin,xbin)
            #find the number of entries that were used for this
            # along with bin mean and std
            #nentries = 0
            #binmean = 0. # h2d.GetMean(2)
            #binrms = 0. #h2d.GetRMS(2)
            projection = h2d.ProjectionY("hx"+str(xbin),xbin+1, xbin+1)
            #projection.GetYaxis().SetRangeUser(0, 10000)
            projection.Draw("hist")
            nentries = projection.GetEntries()
            raw_input("showing projection for xbin "+str(xbin))
            #for ybin in range(nybins):
            #    binmean += h2d.GetBinContent(xbin, ybin)*h2d.GetYaxis().GetBinCenter(ybin)
            #    binrms += h2d.GetBinContent(xbin, ybin) * h2d.GetYaxis().GetBinCenter(ybin)**2
            #    #nentries += h2d.GetBinContent(xbin, ybin)
                
            binmean = projection.GetMean()
            binrms = projection.GetRMS()
            #only need to calculate stuff for the last old (smaller) bin in the new bin
            #if xbin%rebin == 0:
            if nentries == 0: 
                print("error: nentries = 0.") 
                continue
            binmean /= nentries
            binrms /= nentries
            binrms = binrms**0.5
            #raw_input("continue?")
            #fix stupid roundoff error
            if binrms - binmean < .00001: 
                print("error: binrms = %f, binmean = %f"%(binrms, binmean)) 
                continue
            print("bin: %d, mean: %f, rms: %f, nentries: %d"%(xbin, binmean, binrms, nentries)) 
            binerr = ( (binrms**2 - binmean**2) / nentries )**0.5
            #if binerr == 0: continue
            print("bin: %d, error: %f"%(xbin, binerr)) 
            hnew.SetBinContent(xbin, binmean)
            hnew.SetBinError(xbin, binerr)

        ###flip the axes
        ##hnew = ROOT.TH2F("hnew"+str(fitnum), "H vs. E", 5000, 0, 500, 2000, 0, 200.0)
        ##for i in range(h2d.GetNbinsX()):
        ##    for j in range(h2d.GetNbinsY()):
        ##        xcent = h2d.GetXaxis().GetBinCenter(i)
        ##        ycent = h2d.GetYaxis().GetBinCenter(j)
        ##        content = h2d.GetBinContent(i, j)
        ##        hnew.Fill(ycent, xcent, content)
                
        fnew.cd()
        hnew.Write()
    #h2d.GetYaxis().SetLimits(0, 400)
    #now do the fit
    #pt doesn't go lower than 2 GeV
    #xmin = 0 #2 
    #only interested in low pt (below 10 GeV)
    #xmax = 500 #10
    #ymax = 40
    #fitting function
    #f1 = ROOT.TF1("f1", "(gaus(0))", xmin, xmax)
    #can add option Q for quiet fit, option N to now show.
    #h2d.FitSlicesX(f1, 0, -1, 0, "R")
    #slices = ROOT.TObjArray()
    #print("slices made: " + str(slices))
    #h2d.FitSlicesX(f1, 0, -1, 0, "G5 NR", slices)
    #hnew.Draw("colz")
    #hnew.SetMarkerStyle(2)
    cnew = ROOT.TCanvas("cnew"+str(fitnum), "cnew"+str(fitnum))
    cnew.cd()
    #hnew.Draw("hist pe")
    hnew.GetXaxis().SetTitle("E (GeV)")
    hnew.GetYaxis().SetTitle("H (GeV)")
    #raw_input("h")
    f1 = ROOT.TF1("f1", "([0]*x + [1])", xmin, 10.)
    fitresult = hnew.Fit(f1, 'S', "", 0, 10)
    #param 0 is the const multiplier of the gaussian
    #slices[0].Print()
    #slices[0].Draw("hist")
    ROOT.gStyle.SetOptFit(1)
#    hnew.Draw("hist pe")
#    fitresult.Draw("same")
    raw_input("h")

#h2d.Draw()
if not read_fast:
    fnew.Close()
