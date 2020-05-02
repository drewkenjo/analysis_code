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

hel_ptheta = fin.Get("el_ptheta")
hpr_ptheta = fin.Get("pro_ptheta")
hkp_ptheta = fin.Get("kp_ptheta")
hkm_ptheta = fin.Get("km_ptheta")

hel_phitheta = fin.Get("el_phitheta")
hpr_phitheta = fin.Get("pro_phitheta")
hkp_phitheta = fin.Get("kp_phitheta")
hkm_phitheta = fin.Get("km_phitheta")


cel = TCanvas("cel","cel",900,900)
cel.Divide(2,1)
cel.cd(1)
hel_ptheta.SetTitle("Final State El: p vs #theta; p (GeV); #theta (deg)")
hel_ptheta.SetTitle("Final State El: p vs #theta; p (GeV); #theta (deg)")
hel_ptheta.Draw("colz")
cel.cd(2)
hel_phitheta.SetTitle("Final State El: #phi vs #theta; phi (deg); #theta (deg)")
hel_phitheta.SetTitle("Final State El: #phi vs #theta; phi (deg); #theta (deg)")
hel_phitheta.Draw("colz")
cel.SaveAs("final_el_kin.pdf")

cpr = TCanvas("cpr","cpr",900,900)
cpr.Divide(2,1)
cpr.cd(1)
hpr_ptheta.SetTitle("Final State Pr: p vs #theta; p (GeV); #theta (deg)")
hpr_ptheta.SetTitle("Final State Pr: p vs #theta; p (GeV); #theta (deg)")
hpr_ptheta.Draw("colz")
cpr.cd(2)
hpr_phitheta.SetTitle("Final State Pr: #phi vs #theta; phi (deg); #theta (deg)")
hpr_phitheta.SetTitle("Final State Pr: #phi vs #theta; phi (deg); #theta (deg)")
hpr_phitheta.Draw("colz")
cpr.SaveAs("final_pr_kin.pdf")

ckp = TCanvas("ckp","ckp",900,900)
ckp.Divide(2,1)
ckp.cd(1)
hkp_ptheta.SetTitle("Final State Kp: p vs #theta; p (GeV); #theta (deg)")
hkp_ptheta.SetTitle("Final State Kp: p vs #theta; p (GeV); #theta (deg)")
hkp_ptheta.Draw("colz")
ckp.cd(2)
hkp_phitheta.SetTitle("Final State Kp: #phi vs #theta; phi (deg); #theta (deg)")
hkp_phitheta.SetTitle("Final State Kp: #phi vs #theta; phi (deg); #theta (deg)")
hkp_phitheta.Draw("colz")
ckp.SaveAs("final_kp_kin.pdf")

ckm = TCanvas("ckm","ckm",900,900)
ckm.Divide(2,1)
ckm.cd(1)
hkm_ptheta.SetTitle("Final State Km: p vs #theta; p (GeV); #theta (deg)")
hkm_ptheta.SetTitle("Final State Km: p vs #theta; p (GeV); #theta (deg)")
hkm_ptheta.Draw("colz")
ckm.cd(2)
hkm_phitheta.SetTitle("Final State Km: #phi vs #theta; phi (deg); #theta (deg)")
hkm_phitheta.SetTitle("Final State Km: #phi vs #theta; phi (deg); #theta (deg)")
hkm_phitheta.Draw("colz")
ckm.SaveAs("final_km_kin.pdf")




