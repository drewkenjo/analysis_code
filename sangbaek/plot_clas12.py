#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])

# lumi = 1.0558*0.0001
# xsec = 9.0285*100000
# tot_rate = lumi * xsec
# phi_acceptance = 4.5/360
# ratio = tot_rate/ff.Get("spec_h_totalevent").GetBinContent(1)
# ratio = ratio*phi_acceptance
# ratio = ratio*24*3600

ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetMaxDigits(3);

c1 = ROOT.TCanvas('c1','c1',1024,592)
c1.Divide(3,3)

xB_array = [0.1, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35, 0.38, 0.42, 0.58]
t_array = [0.09, 0.13, 0.18, 0.23, 0.3, 0.39, 0.52, 0.72, 1.1, 2]
Integral = 0

h2 = ROOT.TH1F("test","test",24,0,360)

for i in range(1,190):
	c1.cd((i-1)%9+1)
	tbin = (i-1)%9
	xBQbin = (i-1)//9+1
	h1 = ff.Get("root_dvcs_h_phi_pro_CD_gam_FT_bin_"+str(21*tbin + xBQbin))
	if (not h1):
		h1 = h2.Clone(str(i))

	h1.SetTitle("bin"+str((i-1)//9+1)+"\t"+str(t_array[(i-1)%9])+"<-t<"+str(t_array[(i-1)%9+1])+" GeV ^{2}")
	h1.GetXaxis().SetTitle("#phi (#circ)");
	Integral = Integral + h1.Integral()
	if i%9 == 0:
		h1.Draw()
		c1.Update()
		if i == 9:
			c1.Print("clas12_dvcs_raw_yields.pdf(")
		elif i == 189:
			c1.Print("clas12_dvcs_raw_yields.pdf)")		
		else:
			c1.Print("clas12_dvcs_raw_yields.pdf")		
	else:
		h1.Draw()
		c1.Update()

print(Integral)

# c1 = ROOT.TCanvas('c1','c1',480,480)

# for i in range(1,190):
# 	h1 = ff.Get("root_dvcs_h_Q2_xB_bin_"+str(i))
# 	if not h1:
# 		 continue
# 	h1.SetTitle("Q^{2} - x_{B}")
# 	h1.GetXaxis().SetTitle("x_{B}");
# 	h1.GetYaxis().SetTitle("Q^{2} (GeV^{2})");
# 	h1.GetYaxis().SetRangeUser(0.5,5.5);
# 	h1.GetXaxis().SetLimits(.05,.7);
# 	if (i-1)%21==0:
# 		h1.Draw("same")
# 	elif ((i-1)%21)%4 == 0 or ((i-1)%21)%4 == 1:
# 		h1.SetMarkerColor(2)
# 		h1.Draw("same")
# 	elif ((i-1)%21)%4 == 2 or ((i-1)%21)%4 == 3:
# 		h1.SetMarkerColor(4)
# 		h1.Draw("same")

# c1.Print("dvcs_binning_clas12.pdf")
# c1.cd()
# h1 = ff.Get("epg_h_Q2_xB")
# h1.SetTitle("Q^{2} - x_{B}")
# h1.GetXaxis().SetTitle("x_{B}");
# h1.GetYaxis().SetTitle("Q^{2} (GeV^{2})");
# # h1.GetYaxis().SetRangeUser(0,6);
# h1.Draw("colz")
# c1.Print("dvcs_clas12.pdf(")

# c1.cd()
# h1 = ff.Get("epg_h_Q2_xB_W<2")
# h1.SetTitle("Q^{2} - x_{B}")
# h1.GetXaxis().SetTitle("x_{B}");
# h1.GetYaxis().SetTitle("Q^{2} (GeV^{2})");
# h1.GetYaxis().SetRangeUser(0.5,5.5);
# h1.GetXaxis().SetLimits(.05,.7);
# h2 = ff.Get("epg_h_Q2_xB_5_theta_10")
# h2.SetMarkerColor(2)
# h3 = ff.Get("epg_h_Q2_xB_10_theta_15")
# h3.SetMarkerColor(3)
# h4 = ff.Get("epg_h_Q2_xB_15_theta_20")
# h4.SetMarkerColor(4)
# h5 = ff.Get("epg_h_Q2_xB_20_theta_25")
# h5.SetMarkerColor(6)
# h6 = ff.Get("epg_h_Q2_xB_25_theta_30")
# h6.SetMarkerColor(7)
# h7 = ff.Get("epg_h_Q2_xB_30_theta_35")
# h7.SetMarkerColor(28)

# h1.Draw()
# h2.Draw("same")
# h3.Draw("same")
# h4.Draw("same")
# h5.Draw("same")
# h6.Draw("same")
# h7.Draw("same")
# h1.Draw("same")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("epg_h_t_xB")
# h1.SetTitle("t - x_{B}")
# h1.GetXaxis().SetTitle("x_{B}");
# h1.GetYaxis().SetTitle("-t (GeV^{2})");
# # h1.GetYaxis().SetRangeUser(0,6);
# h1.Draw("colz")
# c1.Print("dvcs_clas12.pdf)")

# c1.cd()
# h1 = ff.Get("rates_h_ele_rate")
# h1.SetTitle("electron rates over 4.5 #circ")
# h1.GetXaxis().SetTitle("#theta (#circ)")
# h1.GetYaxis().SetTitle("Counts/ day")
# h1.Scale(ratio)
# h1.Draw("ehist")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("rates_h_pro_rate")
# h1.SetTitle("Proton rates over 4.5 #circ")
# h1.GetXaxis().SetTitle("#theta (#circ)")
# h1.GetYaxis().SetTitle("Counts/ day")
# h1.Scale(ratio)
# h1.Draw("ehist")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("rates_h_gam_rate")
# h1.SetTitle("Photon rates over 4.5 #circ")
# h1.GetXaxis().SetTitle("#theta (#circ)")
# h1.GetYaxis().SetTitle("Counds/ day")
# h1.Scale(ratio)
# h1.Draw("ehist")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("angular_h_ep_azimuth")
# h1.SetTitle("e-p azimuthal angle")
# h1.GetXaxis().SetTitle("#phi_{p} (#circ)")
# h1.GetYaxis().SetTitle("#phi_{e} (#circ)")
# h1.Draw("colz")
# c1.Print("dvcs_clas12.pdf")


# c1.cd()
# h1 = ff.Get("epg_hmm2_epg")
# h1.SetTitle("Missing Mass Squared, ep#gamma")
# h1.GetXaxis().SetTitle("MM^{2}_{ep#gamma} (GeV^{2})")
# h1.Draw("colz")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("angular_h_ep_azimuth_diff")
# h1.SetTitle("e-p azimuthal angle difference")
# h1.GetXaxis().SetTitle("|#phi_{e}-#phi_{p}| (#circ)")
# h1.Draw("colz")
# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("angular_h_ep_polar")
# h1.SetTitle("e-p polar angle")
# h1.GetXaxis().SetTitle("#theta_{p} (#circ)")
# h1.GetYaxis().SetTitle("#theta_{e} (#circ)")
# h1.Draw("colz")

# f1 = ROOT.TF1('f1', '180-2*180/TMath::Pi()*TMath::ATan(5.938/0.938*TMath::Tan(x*TMath::Pi()/180))', 0,90)
# f1.Draw("same")


# c1.Print("dvcs_clas12.pdf")

# c1.cd()
# h1 = ff.Get("xsec_h_cross_section")
# h1.SetTitle("DVCS cross section (rg-a train)")
# h1.GetXaxis().SetTitle("#phi (#circ)")
# h1.GetYaxis().SetTitle("Number of Events")
# print(h1.Integral())
# c1.SetLogy()

# h1.Draw("ehist")
# c1.Print("dvcs_clas12.pdf)")
