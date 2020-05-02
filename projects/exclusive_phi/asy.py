from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from array import array
import math

fin = TFile("../../analysis_out_skim4_rgaf18_final.root","read")

h_pos_asym = fin.Get("pos_asy")
h_neg_asym = fin.Get("neg_asy")
h_pos_asym.RebinX(2)
h_neg_asym.RebinX(2)

h_p_copy = TH1F(h_pos_asym)
h_n_copy = TH1F(h_neg_asym)

h_p_copy2 = TH1F(h_pos_asym)
h_n_copy2 = TH1F(h_neg_asym)

h_p_copy.Add(h_n_copy,-1)
#h_p_copy2.Add(h_n_copy2,0.86)

#h_p_copy.Divide(h_p_copy2)

tren_phi = array('d')
asy = array('d')
tren_phi_er = array('d')
asy_er = array('d')

for i in range(1, h_pos_asym.GetNbinsX()+1):
    print(" bin number %d " % (i) )

    P = h_pos_asym.GetBinContent(i)
    N = h_neg_asym.GetBinContent(i)
    bc = h_pos_asym.GetXaxis().GetBinCenter(i)
    
    y = (P-N)/(P+N)
    y_err = math.sqrt( (1 - y**2)/(P+N))

    print(" number of positive values %d, number of negative values %d" % (P, N))        
    print(" bin center %d, asy %d, bin center error %d, asy error %d " % (bc, y, 0, y_err))

    tren_phi.append(bc)
    asy.append(y*(1/0.85355))
    tren_phi_er.append(0)
    asy_er.append(y_err)


# create fit function for asy
fasy = TF1("fasy","[0]*sin(x*3.1415926/180.0)",0.0,361.0)
fasy.SetParameter(0,1)

c1 = TCanvas("c1", "canvas", 800, 800)
gStyle.SetOptStat(0000)
c1.Divide(2,2)

c1.cd(1)
h_pos_asym.SetTitle("Counts with Pos. Helicity; #phi_{trento}; Counts")
h_pos_asym.Draw()

c1.cd(2)
h_neg_asym.SetTitle("Counts with Neg. Helicity; #phi_{trento}; Counts")
h_neg_asym.Draw()

c1.cd(3)
h_p_copy.SetTitle("Pos.-Neg. Helicity; #phi_{trento}; Counts")
h_p_copy.Draw()

c1.cd(4)
g_asy = TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
g_asy.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy.Fit("fasy")
g_asy.Draw("AP")
fasy.Draw("same")
c1.SaveAs("phi_asy.pdf")

