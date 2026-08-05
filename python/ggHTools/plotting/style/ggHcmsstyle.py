from plotting.style import CMS_lumi
from plotting.style import tdrstyle
import ROOT 

class CMSstyle:
    def __init__(self, c, l, year, lumi, tex_entries=[]):
        self.c=c
        self.l=l
        self.iPeriod=0
        self.iPos=0
        self.tex=ROOT.TLatex()
        self.tex_entries=tex_entries
        self.year=year
        self.lumi=lumi
        tdrstyle.setTDRStyle()
        CMS_lumi.CMS_lumi(self.c, self.iPeriod, self.iPos, self.year, self.lumi)

        self.c.SetTicks(1, 1)
        self.c.SetLeftMargin(0.12)
        self.c.SetRightMargin(0.2)
        self.c.SetBottomMargin(0.18)    
        self.c.GetFrame().SetLineWidth(3)

        self.l.SetBorderSize(0)
        self.l.SetFillStyle(0)
        self.l.SetTextSize(0.03)

        self.tex.SetNDC()
        self.tex.SetTextSize(0.03)
        self.tex.SetTextFont(42)
        i=0
        for entry in self.tex_entries:
            self.tex.DrawLatex(0.2, 0.85-0.05*i, f"{entry}")
            i+=1


        #self.tex.DrawLatex(0.4, 0.85, "LH")
        #self.tex.DrawLatex(0.39, 0.47, "LL")
        #self.tex.DrawLatex(0.73, 0.47, "HL")
        #self.tex.DrawLatex(0.73, 0.85, "HH")

        self.c.Update()
