#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])
f1 = ROOT.TF1('f1', 'gaus(0)+pol0(3)', 0,1)

c1 = ROOT.TCanvas('c1','c1',2200,2400)
c1.Divide(2,3)


# Samplic Fraction
for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1))
	h1.GetXaxis().SetTitle("mom (GeV/c)")
	h1.GetYaxis().SetTitle("Sampling Fraction")
	h1.Draw("colz")
c1.Print("elec.pdf(")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1))
	h1.GetXaxis().SetTitle("mom (GeV/c)")
	h1.GetYaxis().SetTitle("Sampling Fraction")
	h1.Draw("colz")
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Electron, Sampling Fraction, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Electron, Sampling Fraction, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Electron, Sampling Fraction, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_Sampl_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Electron, Sampling Fraction, momentum (7.95,10.6) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Negative, Sampling Fraction, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Negative, Sampling Fraction, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Negative, Sampling Fraction, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_Sampl_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Negative, Sampling Fraction, momentum (7.95,10.6) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")



for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_S"+str(i+1))
	# f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
	# f1.SetParameter(2, h1.GetRMS())
	# mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	# f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	# h1.Fit(f1)
	# mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	# f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	# h1.Fit(f1, 'R')
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_S"+str(i+1))
	# f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
	# f1.SetParameter(2, h1.GetRMS())
	# mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	# f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	# h1.Fit(f1)
	# mu,sig = [f1.GetParameter(i) for i in range(1,3)]

	# f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
	# h1.Fit(f1, 'R')
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Electron, vz, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Electron, vz, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Electron, vz, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_vz_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Electron, vz, momentum (7.95,10.6) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Negative, vz, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Negative, vz, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Negative, vz, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_vz_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Negative, vz, momentum (7.95,10.6) GeV/c")
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
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_nphe_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Electron, nphe, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_nphe_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Electron, nphe, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_nphe_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Electron, nphe, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_nphe_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Electron, nphe, momentum (7.95,10.6) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_nphe_mom_S"+str(i+1)).ProjectionY("",0,24,"[-cutg]")
	h1.SetTitle("Negative, nphe, momentum (0,2.65) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_nphe_mom_S"+str(i+1)).ProjectionY("",25,49,"[-cutg]")
	h1.SetTitle("Negative, nphe, momentum (2.65,5.3) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_nphe_mom_S"+str(i+1)).ProjectionY("",50,74,"[-cutg]")
	h1.SetTitle("Negative, nphe, momentum (5.3,7.95) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")

for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_nphe_mom_S"+str(i+1)).ProjectionY("",75,99,"[-cutg]")
	h1.SetTitle("Negative, nphe, momentum (7.95,10.6) GeV/c")
	h1.Draw()
c1.Print("elec.pdf")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_elec_PCALECAL_S"+str(i+1))
	h1.GetXaxis().SetTitle("E PCAL (GeV)")
	h1.GetYaxis().SetTitle("E ECAL (GeV)")
	h1.Draw("colz")
c1.Print("elec.pdf")


for i in range(0,6):
	c1.cd(i+1)
	h1 = ff.Get("elec_H_neg_PCALECAL_S"+str(i+1))
	h1.GetXaxis().SetTitle("E PCAL (GeV)")
	h1.GetYaxis().SetTitle("E ECAL (GeV)")
	h1.Draw("colz")

c1.Print("elec.pdf)")

