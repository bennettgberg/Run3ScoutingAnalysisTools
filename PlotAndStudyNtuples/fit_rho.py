import ROOT
import math

#can only read fast if hists_new.root exists and is filled in.
read_fast = False
rebin = 20 #what factor to rebin by
#fname = "hists_DoubleElectron_ScoutingSkim220404.root"
#fname = "hists_DoubleElectron_ScoutingSkim220411.root"
#fname = "hists_DoubleElectron_ScoutingSkim220506.root"
#fname = "hists_DYToLL_ScoutingSkim220510.root"
#fname = "hists_QCD_ScoutingSkim220510.root"
#fname = "hists_new.root"
fname = "hists_EtaTo2Mu2E_ScoutingSkim221208.root"
if read_fast:
    fname = "hists_new_Rho.root"
f = ROOT.TFile(fname)
f.Print()
#here's the 2D hist that we have to fit.
#hname = "hoevsptBkg"
#hname = "hVe"
hname = "RhovsH"
if read_fast:
    hname = "hnew"

if not read_fast:
    fnew = ROOT.TFile("hists_new_Rho.root", "recreate") 

#for fitnum in range(0, 7):
#for fitnum in range(7, 8):

#h2d = f.Get(hname+str(fitnum))
h2d = f.Get(hname)
h2d.Print()
#h2d.Draw()
if read_fast:
    hnew = h2d
nbins = h2d.GetXaxis().GetNbins()
print("old nbins: %d"%nbins) 
#make 1 d histogram, x bin vs. mean of all y values in that bin
h2d.RebinX(rebin)
#rebin in Y too???
#h2d.RebinY(rebin)
nbins = h2d.GetXaxis().GetNbins()
print("new nbins: %d"%nbins) 
print("nybins: %d"%h2d.GetYaxis().GetNbins()) 
#new_nbins = nbins / rebin
if not read_fast:
    nybins = h2d.GetYaxis().GetNbins()
xmin = h2d.GetXaxis().GetXmin()
xmax = h2d.GetXaxis().GetXmax()
ymin = h2d.GetYaxis().GetXmin()
ymax = h2d.GetYaxis().GetXmax()
binsize = (xmax - xmin) / nbins
ybinsize = (ymax - ymin) / nybins
print("xmin: %f, xmax: %f, nbins: %d "%(xmin, xmax, nbins))
if not read_fast:
    #hnew = ROOT.TH1F("hnew"+str(fitnum), "hnew"+str(fitnum), new_nbins, xmin, xmax)
    #hnew = ROOT.TH1F("hnew"+str(fitnum), "hnew"+str(fitnum), nbins, xmin, xmax)
    hnew = ROOT.TH1F("hnew", "hnew", nbins, xmin, xmax)
#for each bin, get the mean y value and std dev y value
#  then set the point to be mean +/- stdev 
if not read_fast:
    binmean = 0.
    binrms = 0.
    nentries = 0
    h2d.Draw()
    raw_input("drawing h2d.")
    #if xbin%rebin == 1:
    binmean = 0.
    binrms = 0.
    nentries = 0
        
    h_1 = ROOT.TH1F("H90", "H90", nbins, xmin, xmax)
    for xbin in range(nbins):
    #    ##so instead of doing this projection then get mean etc.,
    #    ## should use FitSlicesX() to get the 90% cutoff of each bin (instead of the 50% cutoff)
    #    ##  then fit those points with a line.
        cut90 = -1
        summie = 0
        projection = h2d.ProjectionY("hx"+str(xbin),xbin+1, xbin+1)
        myN = projection.GetEntries()
        #print("xbin %d, nentries: %d"%(xbin, myN)) 
        for ybin in range(nybins):
            yN = projection.GetBinContent(ybin)
            summie += yN
            #print("ybin: %d; yN: %d; new sum: %d"%(ybin, yN, summie)) 
            if summie >= 0.9*myN:
                cut90 = ybin
                break
        #print("cut90: %d, sum: %d"%(cut90, summie)) 
        if cut90 == -1:
            h1val = 0.0
        else:
            h1val = projection.GetBinCenter(cut90)
        print("h1val: %f"%h1val) 
        h1err = ybinsize
        h_1.SetBinContent( xbin, h1val ) 
        h_1.SetBinError( xbin, h1err )
        h_1.SetTitle("90% cutoff H value vs. rho")
        h_1.GetYaxis().SetTitle("H (GeV)")
        h_1.GetXaxis().SetTitle("#rho") 
    
        #projection.Draw()
        #raw_input("wait....") 
    if 0 == 0:
        #exit()
        #arr = []
        #print("boutta project")
        #ROOT.gDirectory.ls()
        ##why 10k bins instead of 30?????
        ##h2d.FitSlicesX(   )
        ##function to fit each slice to: exponential decay
        #f1 = ROOT.TF1( "f1", "[0]*exp(-x/[1])", 0, 30  )
        #f1.SetParameter(0, 10)
        #f1.SetParameter(1, 10)
        ####start at ybin 1 instead of 0 because huge spike in first bin??-- hmmm doesn't seem to change anything...
        ##ymin = 0.2
        ##ymax = 50
        ##f1.SetRange(ymin,ymax)
        #h2d.FitSlicesY( f1, 1, -1, 0, "R"  )
        ##h2d.FitSlicesY(   )
        #print("done with FitSlicesY()") 
        #ROOT.gDirectory.ls()
        #h_0 = ROOT.gDirectory.Get("RhovsH_0")
        ##param 1 * log(10) should give the 90% cutoff!
        #h_1 = ROOT.gDirectory.Get("RhovsH_1")
        ##h_sigma = ROOT.gDirectory.Get("RhovsH_2")
        #h_chi2 = ROOT.gDirectory.Get("RhovsH_chi2")
        #h2d.Draw()
        #raw_input("wait....") 
        ##print("h_mean: " + str(h_mean)) 
        #h_0.Draw()
        #raw_input("wait....") 
        #h_1.Draw()
        #raw_input("wait....") 
        #h_chi2.Draw()
        #raw_input("wait....") 
        #now fit the resulting hist to a quadratic!
            #NO must be linear!! only question is what range to use...
        ##Quadratic fit!
        #f2 = ROOT.TF1( "f2", "[0]+[1]*x+[2]*x*x" )
        #f2.SetParameters(4., -1, .25)
        h_1.SetMarkerStyle(8)
        #Linear fit!
        f2 = ROOT.TF1( "f2", "[0]+[1]*x" )
        #h_1.Fit(f2, "", "", 0., 18.)
        h_1.Fit(f2, "", "", 18., 40.)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(11)
        h_1.Draw('PE')
        raw_input("wait....")
        exit()
        #binmean = projection.GetMean()
        #binrms = projection.GetRMS()
        nentries = h_1.GetEntries()
        ##only need to calculate stuff for the last old (smaller) bin in the new bin
        ##if xbin%rebin == 0:
        if nentries == 0: 
            print("error: nentries = 0.") 
            #continue
        #binrms /= nentries
        #binrms = binrms**0.5
        #raw_input("continue?")
        #fix stupid roundoff error
        if binrms - binmean < .00001: 
            print("error: binrms = %f, binmean = %f"%(binrms, binmean)) 
            #continue
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
    #cnew = ROOT.TCanvas("cnew"+str(fitnum), "cnew"+str(fitnum))
    cnew = ROOT.TCanvas("cnew", "cnew")
    cnew.cd()
    #hnew.Draw("hist pe")
    #hnew.GetXaxis().SetTitle("E (GeV)")
    hnew.GetXaxis().SetTitle("pileup #rho")
    hnew.GetYaxis().SetTitle("H (GeV)")
    #raw_input("h")
    #f1 = ROOT.TF1("f1", "([0]*x + [1])", xmin, 10.)
    #quadratic instead of linear??
    f1 = ROOT.TF1("f1", "([0]*x + [1] + [2]*x*x)", xmin, 50.)
    #fitresult = hnew.Fit(f1, 'S', "", 0, 10)
    fitresult = hnew.Fit(f1, 'S', "", 0, 50)
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
