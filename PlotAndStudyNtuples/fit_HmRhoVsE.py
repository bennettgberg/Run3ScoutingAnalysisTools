import ROOT

fname = "hists_EtaTo2Mu2E_ScoutingSkim221208.root"
f = ROOT.TFile(fname)

hname = "HmRhoVsE"

h2d = f.Get(hname)

rebin = 5
h2d.RebinX(rebin)
h2d.RebinY(rebin)

nbins = h2d.GetXaxis().GetNbins()
xmin = h2d.GetXaxis().GetXmin()
xmax = h2d.GetXaxis().GetXmax()

h_1 = ROOT.TH1F("HmRhoMean", "HmRhoMean", nbins, xmin, xmax)
#param obtained from fit_rho fit of H vs. rho
Crho = 0.07

#first get the average within each E bin
for xbin in range(nbins):
    projection = h2d.ProjectionY("hx"+str(xbin),xbin+1, xbin+1)
    my_mean = projection.GetMean()
    my_rms = projection.GetRMS()
    entries = projection.GetEntries()
    if entries > 0:
        my_sdev = my_rms / entries**0.5
    else:
        my_sdev = my_mean
    print("xbin %d; mean %f, rms %f, entries %d, sdev %f"%(xbin, my_mean, my_rms, entries, my_sdev)) 
    h_1.SetBinContent(xbin, my_mean)
    h_1.SetBinError(xbin, my_sdev)

#define function to fit to
#f1 = ROOT.TF1("f1", "([0] + [1]/(x-[2]) + [3]/x)", xmin, 350.)
f1 = ROOT.TF1("f1", "([0] + [1]/(x-[2]) + [3]/x+[4]/(x*x)+[5]/(x-[6]))", xmin, 400.)
#f1.SetParameter(0, .05)
#f1.SetParameter(1, 350)
#f1.SetParameter(2, -1200)
#f1.SetParameter(-3, .05)

#f1 = ROOT.TF1("f1", "([0] + ([1]*x*x+[2]*x+[3])/(x*x))", xmin, 500.)
fitresult = h_1.Fit(f1, 'S', "", 0., 400.) 

h_1.GetXaxis().SetTitle("E_{SC} (GeV)")
h_1.GetYaxis().SetTitle("Mean (H-C_{#rho}*#rho)/E")
h_1.SetMarkerStyle(8)
h_1.Draw()
raw_input("....")
