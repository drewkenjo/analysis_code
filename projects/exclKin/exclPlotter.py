from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from array import array
import math

#plot the final kinematic distributions for the particles in exclusive phi production
gStyle.SetOptStat(0000)
gStyle.SetPalette(55);
gStyle.SetStatStyle(1001); 
gStyle.SetStatFont(62);  
gStyle.SetStatBorderSize(2);
gStyle.SetTitleX(.20);
gStyle.SetTitleY(.99);
gStyle.SetPadBottomMargin(0.1);
gStyle.SetPadLeftMargin(0.12);   
gROOT.ForceStyle()

fin = TFile("../../analysis_out_skim4_rgaf18_final.root","read")

hw = fin.Get("W")
hmissemm2 = fin.Get("missemm2")
hphim = fin.Get("phim")
hmm2phim = fin.Get("missephim")
hmismm2phim = fin.Get("mismm2phim")
hmissmm2 = fin.Get("missmm2")

cmm = TCanvas("cmm","cmm",900,900)
cmm.Divide(2,2)
cmm.cd(1)
hmissemm2.Draw("colz")
cmm.cd(2)
hmm2phim.Draw("colz")
cmm.cd(3)
hmismm2phim.Draw("colz")
cmm.cd(4)
hmissmm2.Draw("colz")
cmm.SaveAs("final_kin_integrated_over_sect_phi.pdf")

hq2x = fin.Get("q2x")
htxb = fin.Get("txb")
hxb = fin.Get("xb")
htq2 = fin.Get("tq2")
ht = fin.Get("t")
hq2 = fin.Get("q2")
hcoskp = fin.Get("cos_kp")

ckin = TCanvas("ckin","ckin",900,900)
ckin.Divide(2,2)
ckin.cd(1)
hxb.SetTitle("X_{b}; X_{b}; counts")
hxb.Draw()
ckin.cd(2)
hq2.SetTitle("Q^{2} ; Q^{2} (GeV^{2}); counts")
hq2.Draw()
ckin.cd(3)
ht.SetTitle("-t; -t (GeV^{2}); counts")
ht.Draw()
ckin.SaveAs("final_kin_integrated_over_sect_phi.pdf")

ckin2 = TCanvas("ckin2","ckin2",900,900)
ckin2.Divide(2,2)
ckin2.cd(1)
hq2x.SetTitle("Q^{2} vs x_{b}; x_{b}; Q^{2} (GeV^{2})")
hq2x.Draw("colz")
ckin2.cd(2)
htxb.SetTitle(" x_{b} vs -t ; -t^{2} (GeV^{2}); x_{b}")
htxb.Draw("colz")
ckin2.cd(3)
htq2.SetTitle("Q^{2} vs -t; -t (GeV^{2}); Q^{2} (GeV^{2})")
htq2.Draw("colz")
ckin2.SaveAs("final_2dkin_integrated_over_sect_phi.pdf")




