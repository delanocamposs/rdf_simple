import ROOT
from plotting.style import ggHcmsstyle
writeExtraText=True

c= ROOT.TCanvas("c", "CMS Style Polynomial", 800, 700)

leg=ROOT.TLegend(0.3, 0.75, 0.5, 0.85)
poly = ROOT.TF1("poly", "[0]*x*x", 1e-3, 1)
polyy = ROOT.TF1("polyy", "[0]*(2*x*(1-x)+x*x)", 1e-3, 1)

#from xsec_{ggH+VBF} and also included efficiency of the trigger in signal MC
poly.SetParameters(52.143*0.9768)
poly.SetLineWidth(3)
poly.SetLineColor(ROOT.kAzure+5)
polyy.SetLineColor(ROOT.kRed)
poly.GetXaxis().SetTitle("BR(#phi #rightarrow #gamma#gamma)")
poly.GetXaxis().SetLimits(1e-3,1)
poly.GetYaxis().SetTitle("#frac{Total}{BR(H #rightarrow #phi#phi)}")
poly.GetYaxis().SetTitleSize(0.01)

#from EXO-23-012 using xsec_{VH}*BR(V-->e/mu)
#xsec_{ZH}*[BE(Z-->ee)+BR(Z-->mm)] + xsec_{W+/-,H}*[BR(W+/- --> ev)+BR(W+/- --> mv)]=0.35858
polyy.SetParameters(0.358583457)

c.Update()
e1 = leg.AddEntry("poly", "This work (ggH)", "l")
e2=leg.AddEntry("polyy", "VH", "l")
e1.SetLineWidth(3)
e1.SetLineColor(ROOT.kAzure+5)
e2.SetLineWidth(3)
e2.SetLineColor(ROOT.kRed)
leg.SetTextSize(0.1)
poly.Draw()
polyy.Draw("same")
leg.Draw("Same")
ggHcmsstyle.CMSstyle(c, leg, [""])

c.SetLeftMargin(0.2)
c.SetLogy()
c.SetLogx()
c.Update()
c.SaveAs("polynomial_CMS.pdf")
