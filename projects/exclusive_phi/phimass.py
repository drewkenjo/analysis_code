from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TMath, TFormula
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from array import array
import math

fin = TFile("../../analysis_out_skim4_rgaf18_final.root","read")

phim = fin.Get("phim")



gaus = "[0]*TMath::Gaus(x,[1],[2])"
hvyside = "(x>[4])*2/TMath::Pi()*TMath::ATan( (x-[4])/[5])"
bck = "[3]*TMath::Exp(-TMath::Power(TMath::Abs([6]*(x-[4])),[7]))"
#form = TFormula("form",gaus +"+"+ hvyside + "*" + bck)
form = gaus +"+"+ hvyside + "*" + bck

#print(" my function to fit phi is " + form)
phif = TF1("fitPhi",form, 0.98, 1.3)
print(phif)
fitparm = [900.2, 1.02, 0.022, 250, 0.987, 0.0045, 2.19, 3]
for pp in range(len(fitparm)):
    print("par number %d, value of %f" %(pp, fitparm[pp]) )
    phif.SetParameter(pp,fitparm[pp])
    phif.SetParLimits(pp, 0.5*fitparm[pp], 1.5*fitparm[pp])

phif.SetNpx(10000)    
phif.FixParameter(4,fitparm[4])
phim.Fit("fitPhi","BRN")

bck = TF1("bck",form,0.98,1.3)
bck.SetNpx(10000)
for pp in range(len(fitparm)):
    print("par number %d, value of %f" %(pp, fitparm[pp]) )
    bck.SetParameter(pp,phif.GetParameter(pp))

bck.FixParameter(0,0)


c1 = TCanvas("c1","c1",900,500)
c1.cd(1)
phim.SetTitle("Exclusive #phi from epK^{+}K^{-}; mass (GeV); counts")

phim.GetXaxis().CenterTitle()
phim.GetYaxis().CenterTitle()
phim.Draw("PE")
bck.SetLineColor(4)
phif.Draw("same")
bck.Draw("same")
c1.SaveAs("phiPeak.pdf")
