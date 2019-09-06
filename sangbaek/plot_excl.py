#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])
f1 = ROOT.TF1('f1', 'gaus(0)+pol0(3)', 0,1)

c1 = ROOT.TCanvas('c1','c1',1100,800)

c1.cd()
h1 = ff.Get("epg_hmm2_ep")
h1.SetTitle("Missing Mass Squared, ep")
h1.GetXaxis().SetTitle("MM^2_{ep} (GeV^2)")
h1.Draw("colz")
c1.Print("dvcs.pdf(")


c1.cd()
h1 = ff.Get("epg_hmm2_eg")
h1.SetTitle("Missing Mass Squared, e#gamma")
h1.GetXaxis().SetTitle("MM^2_{e#gamma} (GeV^2)")
h1.Draw("colz")
c1.Print("dvcs.pdf")

c1.cd()
h1 = ff.Get("epg_hmm2_epg")
h1.SetTitle("Missing Mass Squared, ep#gamma")
h1.GetXaxis().SetTitle("MM^2_{ep#gamma} (GeV^2)")
h1.Draw("colz")
c1.Print("dvcs.pdf")

c1.cd()
h1 = ff.Get("epg_hangle_epg")
h1.SetTitle("Angle between #gamma and epX")
h1.GetXaxis().SetTitle("Angle between #gamma and epX (deg)")
h1.Draw("colz")
c1.Print("dvcs.pdf)")
