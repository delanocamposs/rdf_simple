import ROOT 
from plotting.style.ggHcmsstyle import CMSstyle
import numpy

def plot_2D():
    c=ROOT.TCanvas("", "", 800, 800)
    right=0.05
    left=0
    up=0.05
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    temps=[]
    xmax=110
    ymax=110
    xmin=-10
    ymin=-10
    nbins=11

    with open('asymptotic_limits_avg.txt') as f:
        for line in f:
            temps.append(float(line.strip()))

    histo=ROOT.TH2F("hist", "", nbins, xmin, xmax, nbins, ymin, ymax)
    for ix in range(nbins):
        for iy in range(nbins):
            index=ix*nbins+iy
            x=ix*((xmax-abs(xmin))/(nbins-1))
            y=iy*((ymax-abs(ymin))/(nbins-1))
            histo.Fill(x, y, 1e-4*temps[index])
            

    histo.GetXaxis().SetTitle("#phi_{1} L_{xy} (cm)        ")
    histo.GetYaxis().SetTitle("   #phi_{2} L_{xy} (cm)")
    histo.GetZaxis().SetTitle("95% CL Median Upper Limit")
    ROOT.gStyle.SetPalette(ROOT.kBird)


    c.Update()
    histo.Draw("COLZ")
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.25)
    CMSstyle(c, l, ["H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
    histo.SetStats(0)
    histo.GetXaxis().SetTitleOffset(2)
    histo.GetXaxis().SetLabelSize(0.02)
    histo.GetYaxis().SetLabelSize(0.02)
    histo.GetZaxis().SetLabelSize(0.02)
    histo.GetZaxis().SetTitleSize(0.03)
    histo.GetZaxis().SetTitleOffset(2)
    histo.GetYaxis().SetTitleOffset(2)
    c.Update()
    c.SaveAs("limit_heatmap_2D.png")
    c.SaveAs("limit_heatmap_2D.pdf")
    input()



def plot_1D():
    c=ROOT.TCanvas("", "", 800, 800)
    right=0.05
    left=0
    up=0.05
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    temps=[]
    xmax=105
    xmin=-25
    nbins=13

    with open('asymptotic_limits_avg.txt') as f:
        for line in f:
            temps.append(float(line.strip()))

    histo=ROOT.TH1F("hist", "", nbins, xmin, xmax)
    xs=[]
    CLs=[]
    for ix in range(nbins):
        x=10*ix-20
        xs.append(x)
        CLs.append(temps[ix])

    for i in range(len(xs)):
        histo.Fill(xs[i], CLs[i])

    histo.GetXaxis().SetTitle("L_{xy} threshold (cm)")
    histo.GetYaxis().SetTitle("95% CL Median Upper Limit (#times 10^{-4})")
    ROOT.gStyle.SetPalette(ROOT.kBird)


    c.Update()

    c.SetLeftMargin(0.08)
    c.SetRightMargin(0.2)


    histo.Draw("HIST")
    histo.SetFillColor(ROOT.kAzure-9)
    histo.SetLineColor(ROOT.kAzure+2)
    histo.SetLineWidth(2)
    histo.SetFillStyle(1001)
    #c.SetLeftMargin(0.2)
    #c.SetRightMargin(0.25)
    CMSstyle(c, l, ["H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
    histo.SetStats(0)
    #histo.GetXaxis().SetTitleOffset(2)
    #histo.GetXaxis().SetLabelSize(0.03)
    #histo.GetXaxis().SetLabelSize(0.03)
    #histo.GetXaxis().SetNdivisions(810)
    #histo.GetYaxis().SetLabelSize(0.02)
    #histo.GetZaxis().SetLabelSize(0.02)
    #histo.GetZaxis().SetTitleSize(0.03)
    #histo.GetZaxis().SetTitleOffset(2)
    #histo.GetYaxis().SetTitleOffset(2)
    c.Update()
    c.SaveAs("limit_heatmap_1D.png")
    c.SaveAs("limit_heatmap_1D.pdf")
    input()

if __name__=="__main__":
    plot_1D()
