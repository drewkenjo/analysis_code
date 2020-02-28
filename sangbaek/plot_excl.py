#!/usr/bin/python

import sys
sys.argv.append('-b')

import rootplot.root2matplotlib as r2m
import ROOT
from matplotlib import pyplot as plt
import numpy as np

import matplotlib
from matplotlib import rc
import matplotlib.patches as patches
from scipy import stats
from scipy.odr import *
import math
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
# plt.style.use('seaborn-paper')
rc('text',usetex=True)

pgf_with_latex = {
    "pgf.texsystem": "xelatex",     # Use xetex for processing
    "text.usetex": True,            # use LaTeX to write all text
    "font.family": "Helvetica",         # use serif rather than sans-serif
    "font.serif": "Helvetica",
    "font.sans-serif": [],          # Unset sans-serif
    "font.monospace": "Helvet",
    "axes.labelsize": 12,
    "font.size": 12,
    "legend.fontsize": 14,
    "axes.titlesize": 13,           # Title size when one figure
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 12,         # Overall figure title
    "pgf.rcfonts": False,           # Ignore Matplotlibrc
    "text.latex.unicode": True,
    "pgf.preamble": [
        # r'\usepackage{xcolor}',     # xcolor for colours
        # r'\usepackage{fontspec}',
        # r'\usepackage{color}',     # xcolor for colours
        #r'\usepackage{sans-serif}',
        # r'\setmainfont{Helvetica}',
        # r'\setmonofont{Helvetica}',
        #r'\usepackage{unicode-math}',

        # r'\setmathfont{Helvetica}'
    ]
}
matplotlib.rcParams.update(pgf_with_latex)


ff = ROOT.TFile(sys.argv[1])
f1 = ROOT.TF1('f1', 'gaus(0)+pol0(3)', 0,1)

# c1 = ROOT.TCanvas('c1','c1',1100,800)

# c1.cd()
# h1 = ff.Get("epg_hmm2_ep")
# h1.SetTitle("Missing Mass Squared, ep")
# h1.GetXaxis().SetTitle("MM^2_{ep} (GeV^2)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf(")


# c1.cd()
# h1 = ff.Get("epg_hmm2_eg")
# h1.SetTitle("Missing Mass Squared, e#gamma")
# h1.GetXaxis().SetTitle("MM^2_{e#gamma} (GeV^2)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_hmm2_epg")
# h1.SetTitle("Missing Mass Squared, ep#gamma")
# h1.GetXaxis().SetTitle("MM^2_{ep#gamma} (GeV^2)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_hangle_epg")
# h1.SetTitle("Angle between #gamma and epX")
# h1.GetXaxis().SetTitle("Angle between #gamma and epX (deg)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_hangle_ep_eg")
# h1.SetTitle("Angle between two planes, ep and e#gamma")
# h1.GetXaxis().SetTitle("Angle between two planes, ep and e#gamma (deg)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_h_kine_ele")
# h1.GetXaxis().SetTitle("#theta (deg)")
# h1.GetYaxis().SetTitle("momentum (GeV/c)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_h_kine_pro")
# h1.GetXaxis().SetTitle("#theta (deg)")
# h1.GetYaxis().SetTitle("momentum (GeV/c)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_h_kine_gam")
# h1.GetXaxis().SetTitle("#theta (deg)")
# h1.GetYaxis().SetTitle("momentum (GeV/c)")
# h1.Draw("colz")
# c1.Print("dvcs.pdf")

# c1.cd()
# h1 = ff.Get("epg_h_Q2_xB")
# h1.SetTitle("")
# h1.GetXaxis().SetTitle("x_{B}");
# h1.GetYaxis().SetTitle("Q^{2} (GeV^2)");
# a = np.array([0], dtype=np.int32)
# h1.Draw("colz")
# ROOT.gStyle.SetPalette(1)
# ROOT.gPad.Modified()
# ROOT.gPad.Update()
# c1.Print("dvcs.pdf)")

h1 = ff.Get("spec_h_totalevent")
tot= h1.GetBinContent(1)

h1 = ff.Get("epg_hmm2_epg")
h1.SetTitle("Missing Mass Squared, ep"+"$\gamma$")
h1.GetXaxis().SetTitle("MM"+"${}^{2}_{ep}$"+"${}_{\gamma}$"+" (GeV"+"${}^{2}$"+")")

f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
f1.SetParameter(2, h1.GetRMS())
mu,sig = [f1.GetParameter(i) for i in range(1,3)]

f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
h1.Fit(f1)
mu,sig = [f1.GetParameter(i) for i in range(1,3)]

f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
h1.Fit(f1, 'R')

# Make a figure with width 6 inches and height 4 inches
plt.figure(1, (8, 6))
# Create an axes instance
ax1 = plt.axes()
hist = r2m.Hist(h1)

hist.hist(color='white', histtype='stepfilled')
hist.show_titles()

x= np.linspace(f1.GetParameter(1)-2.5*f1.GetParameter(2),f1.GetParameter(1)+2.5*f1.GetParameter(2),1000)
y=f1.GetParameter(0)*np.exp(-0.5*((x-[f1.GetParameter(1)])/f1.GetParameter(2))**2)+f1.GetParameter(3)

ax1.plot(x,y, color='r',label=''.format())
the_table = plt.table(cellText=[['{:1.0f}'.format(tot)],['{:-1.2f}'.format(f1.GetParameter(0))],['{:-1.5f}'.format(f1.GetParameter(1))],['{:-1.5f}'.format(f1.GetParameter(2))]],
                  colWidths = [0.1],
                  rowLabels=['events','amp',' $\mu$','$\sigma$'],
                  colLoc='center', rowLoc='center', cellLoc='center',
                  bbox=[0.8,0.37,0.2,0.2])
the_table.auto_set_font_size(False)
the_table.set_fontsize(18)
the_table.scale(2, 2)

# plt.tight_layout()

plt.savefig('dvcs.pdf')
# plt.show()
