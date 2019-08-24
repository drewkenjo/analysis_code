#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])
f1 = ROOT.TF1('f1', 'gaus(0)+pol0(3)', 0,1)

c1 = ROOT.TCanvas('c1','c1',2200,2400)
c1.Divide(2,3)


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1)).ProjectionX()
	h1.SetTitleX("mom (GeV/c)")
	h1.SetTitleY("Sampling Fraction")
	h1.Draw()
c1.Print("elec.pdf(")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1)).ProjectionX()
	h1.SetTitleX("mom (GeV/c)")
	h1.SetTitleY("Sampling Fraction")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_S"+str(i+1))
	f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
	f1.SetParameter(2, h1.GetRMS())
	mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	h1.Fit(f1)
	mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	h1.Fit(f1, 'R')
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_S"+str(i+1))
	f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
	f1.SetParameter(2, h1.GetRMS())
	mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	h1.Fit(f1)
	mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	h1.Fit(f1, 'R')
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_nphe_S"+str(i+1))
	h1.Draw()
c1.Print("elec.pdf")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_nphe_S"+str(i+1))
	h1.Draw()
c1.Print("elec.pdf)")