#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])

lumi = 1.0558*0.0001
xsec = 9.0285*100000
tot_rate = lumi * xsec
phi_acceptance = 4.5/360
ratio = tot_rate/ff.Get("spec_h_totalevent").GetBinContent(1)
ratio = ratio*phi_acceptance
ratio = ratio*24*3600

ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetMaxDigits(3);


c1 = ROOT.TCanvas('c1','c1',1100,800)

c1.cd()
h1 = ff.Get("epg_h_Q2_xB")
h1.SetTitle("Q^{2} - x_{B}")
h1.GetXaxis().SetTitle("x_{B}");
h1.GetYaxis().SetTitle("Q^{2} (GeV^{2})");
h1.GetYaxis().SetRangeUser(0,6);
h1.Draw("colz")
c1.Print("dvcs_desy.pdf(")

c1.cd()
h1 = ff.Get("rates_h_ele_rate")
h1.SetTitle("electron rates over 4.5 #circ")
h1.GetXaxis().SetTitle("#theta (#circ)")
h1.GetYaxis().SetTitle("Counts/ day")
h1.Scale(ratio)
h1.Draw("ehist")
c1.Print("dvcs_desy.pdf")

c1.cd()
h1 = ff.Get("rates_h_pro_rate")
h1.SetTitle("Proton rates over 4.5 #circ")
h1.GetXaxis().SetTitle("#theta (#circ)")
h1.GetYaxis().SetTitle("Counts/ day")
h1.Scale(ratio)
h1.Draw("ehist")
c1.Print("dvcs_desy.pdf")

c1.cd()
h1 = ff.Get("rates_h_gam_rate")
h1.SetTitle("Photon rates over 4.5 #circ")
h1.GetXaxis().SetTitle("#theta (#circ)")
h1.GetYaxis().SetTitle("Counds/ day")
h1.Scale(ratio)
h1.Draw("ehist")
c1.Print("dvcs_desy.pdf")

# c1.cd()
# h1 = ff.Get("angular_h_ep_azimuth")
# h1.SetTitle("e-p azimuthal angle")
# h1.GetXaxis().SetTitle("#phi_{p} (#circ)")
# h1.GetYaxis().SetTitle("#phi_{e} (#circ)")
# h1.Draw("colz")
# c1.Print("dvcs_desy.pdf")

c1.cd()
h1 = ff.Get("angular_h_ep_azimuth_diff")
h1.SetTitle("e-p azimuthal angle difference")
h1.GetXaxis().SetTitle("|#phi_{e}-#phi_{p}| (#circ)")
h1.Draw("colz")
c1.Print("dvcs_desy.pdf")

c1.cd()
h1 = ff.Get("angular_h_ep_polar")
h1.SetTitle("e-p polar angle")
h1.GetXaxis().SetTitle("#theta_{p} (#circ)")
h1.GetYaxis().SetTitle("#theta_{e} (#circ)")
h1.Draw("colz")

f1 = ROOT.TF1('f1', '180-2*180/TMath::Pi()*TMath::ATan(5.938/0.938*TMath::Tan(x*TMath::Pi()/180))', 0,90)
f1.Draw("same")


c1.Print("dvcs_desy.pdf)")
