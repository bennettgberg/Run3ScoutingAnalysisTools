
import ROOT
from array import array
import glob

newfilename = "electronsHE.root"

#make new output file
newfile = ROOT.TFile.Open(newfilename, "recreate")
newt = ROOT.TTree("eleTree", "eleTree")

#array for H
Harr = array('f', [0.])
#array for E
Earr = array('f', [0.])
#array for classification
# 0: QCD (background)
# 1: DYToLL (signal)
classArr = array('l', [0]) 

newt.Branch("H", Harr, "H/f")
newt.Branch("E", Earr, "E/f")
newt.Branch("class", classArr, "class/l")

#path to the ntuple names
path = "/eos/user/b/bgreenbe/scouting/ntuples/"
#QCD name, then DY name
classNames = ["QCD_ScoutingSkim220510", "DYToLL_ScoutingSkim220510"]

for num,cn in enumerate(classNames):
    print("CLASS NAME: " + cn) 
    fullpath = path + cn 
    fullfile = fullpath + "/" + cn + "_*.root"

    for fname in glob.glob(fullfile):
        print("filename: " + fname)
        f = ROOT.TFile.Open(fname)
        events = f.Get("mmtree/tree")

        for i,e in enumerate(events):
            for j in range(e.n_ele):
                #H/E
                hoe = e.Electron_hoe[j]
                
                #energy matrix
                energyMatrix = e.Electron_energymatrix[j]

                #get E from the energyMatrix
                E = 0.
                for en in energyMatrix:
                    E += en
                
                #H = H/E*E
                H = hoe * E
                
                #now fill the arrays
                Earr[0] = E
                Harr[0] = H
                classArr[0] = num
                
                #now fill the tree
                newt.Fill()

#now write the output file
newfile.cd()
newt.Write()
newfile.Close()
