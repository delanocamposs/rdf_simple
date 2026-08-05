import ROOT
import subprocess
from plotting.style import tdrstyle
from plotting.style import CMS_lumi
from datacard.ggHdatacardmaker import main
from ggHparameters import lumi, order_fit, signal_window
from plotting.plottingtools import fetchError, getPoisson, getPoisson2


def _load_prefit_band(finalstate, physics, mass, lifetime, cat, year):
    path = f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_fit.root"
    cf = ROOT.TFile(path)
    if not cf or cf.IsZombie():
        return None
    g1 = cf.Get("band_1sigma")
    g2 = cf.Get("band_2sigma")
    if not g1 or not g2:
        cf.Close()
        return None
    g1c = g1.Clone("pf_1sigma")
    g2c = g2.Clone("pf_2sigma")
    cf.Close()
    return g2c, g1c

def plot(MultiDimFit, fitDiagnosticsTest, cat, year, bins, finalstate, physics, mass, lifetime, order=order_fit):
    loop=0
    for j in range(2):
        tdrstyle.setTDRStyle()
        CMS_lumi.writeExtraText = True
        CMS_lumi.extraText="Work in Progress"
        CMS_lumi.lumi_13TeV=f"{year}, {lumi[year]/1000} fb^{-1}"
            
        ROOT.gStyle.SetEndErrorSize(2)
        f1 = ROOT.TFile(MultiDimFit)
        f2 = ROOT.TFile(fitDiagnosticsTest)
        r1 = f2.Get("fit_s")
        r2 = f2.Get("fit_b")
        w = f1.Get("w")
        
        #set the proportions and make it beautiful
        can = ROOT.TCanvas("c", "c", 700, 700)
        pad1 = ROOT.TPad("pad1", "pad1", 0,   0.3, 1, 1.0)
        pad1.SetTopMargin(0.08506945)
        pad1.SetBottomMargin(0.00)  
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, 0.25)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()
        pad1.cd()
        x = w.var("mass")
        x.setBins(bins[0])
        plot = x.frame()

        #create data_obs as a TH1 obj so i can look through the contents and throw away the error bars (roofit fucks them up. i compute them manually later for the case of N=0 bins):
        data_obs = w.data("data_obs")
        data_obs_binned = ROOT.RooDataHist("data_obs_binned", "binned data", ROOT.RooArgSet(x), data_obs)
        data_obs_TH1 = data_obs_binned.createHistogram("mass")

        cs=[]
        cloned_data = data_obs_TH1.Clone("int_hist")
        for i in range(1, cloned_data.GetNbinsX()+1):
            count = cloned_data.GetBinContent(i)
            cs.append(cloned_data.GetBinContent(i))
            cloned_data.SetBinError(i, 0)
        cloned_data_binned = ROOT.RooDataHist("data_obs_binned", "binned hist", ROOT.RooArgList(x), cloned_data)

        #create points at the locations of the data bins with correct error bars from fetcError earlier in the code:
        gres1=ROOT.TGraphAsymmErrors()
        q = (1-0.6827)/2.0
        bin_centers = [cloned_data.GetBinCenter(i) for i in range(1, cloned_data.GetNbinsX()+1)]
        for n in range(cloned_data.GetNbinsX()):
            b_n = bin_centers[n]
            c_n = cs[n]
            error=fetchError(q, c_n)
            gres1.SetPoint(n, b_n, c_n)
            gres1.SetPointError(n, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))

        gres1.SetLineColor(ROOT.kBlack)
        gres1.SetMarkerColor(ROOT.kBlack)
        gres1.SetMarkerStyle(20)

        #grab the s+b and b models from the RooWorkspace and get the normalizations of the models after the fits
        sb_model  = w.pdf("model_s").getPdf(f"{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}")
        b_model   = w.pdf("model_b").getPdf(f"{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}")
        if loop==0:
            w.loadSnapshot("MultiDimFit")
        norm_name = "norm_fit_s" if loop == 0 else "norm_prefit"
        norm_fit = f2.Get(norm_name)
        if not isinstance(norm_fit, ROOT.RooArgSet):
            print(f"skipping postfit plot for {cat} ({year}): combine fit failed, no {norm_name} normalizations")
            f1.Close()
            f2.Close()
            return
        bkg_norm = norm_fit.find(f"{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}/background").getVal()
        sig_norm = norm_fit.find(f"{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}/signal").getVal()
        fit_label = "postfit" if loop==0 else "prefit"
        print(f"[{fit_label}] SIGNAL YIELD in [{signal_window[0]},{signal_window[1]}]: {sig_norm:.4f}")
        print(f"[{fit_label}] BKG YIELD in [{signal_window[0]},{signal_window[1]}]: {bkg_norm:.4f}")

        #RooFormulaVar's that scale the pdf's from the RooWorkspace based on the norms of the fits
        b_func = ROOT.RooFormulaVar("bkg_func","N_bkg * B(x)",f"1*@0",ROOT.RooArgList(b_model))
        sb_func = ROOT.RooFormulaVar("sb_func","N_tot*sb_model",f"1*@0",ROOT.RooArgList(sb_model))

        #plotting. the oreder here matters alot in order to get the brazil plot colors to show up correctly and to get the data points on top of everything
        cloned_data_binned.plotOn(plot,ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"), XErrorSize=0, DataError=None)

        prefit_band = _load_prefit_band(finalstate, physics, mass, lifetime, cat, year) if loop==1 else None
        if loop==0 and data_obs_TH1.Integral()!=0:
            try:
                b_model.plotOn(plot,ROOT.RooFit.VisualizeError(r2, 2, ROOT.kFALSE),ROOT.RooFit.Normalization(bkg_norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.FillColor(ROOT.kYellow),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_2sigma"),ROOT.RooFit.DrawOption("F"))
                b_model.plotOn(plot,ROOT.RooFit.VisualizeError(r2, 1, ROOT.kFALSE),ROOT.RooFit.Normalization(bkg_norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.FillColor(ROOT.kGreen),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_1sigma"),ROOT.RooFit.DrawOption("F"))
                show_bands=True
            except Exception:
                print("skipping uncertainty bands: fit covariance unusable (likely a degenerate fit)")
                show_bands=False
        elif loop==1 and data_obs_TH1.Integral()!=0 and prefit_band is not None:
            show_bands=True
        else:
            show_bands=False
            if loop==0:
                print("data is 0. ignoring uncertainty bands because uncertainties on fit parameters are unstable")

        b_model.plotOn(plot,ROOT.RooFit.Normalization(bkg_norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Name("bkg_curve"))
        sb_model.plotOn(plot,ROOT.RooFit.Normalization(sig_norm + bkg_norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sb_curve"))
        cloned_data_binned.plotOn(plot, ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"), XErrorSize=0, DataError=None)

        integral_sb = sb_model.getVal(ROOT.RooArgSet(x))
        integral_b = b_model.getVal(ROOT.RooArgSet(x))

        y_ax_val=round((bins[2]-bins[1])/(bins[0]),2)
        plot.GetYaxis().SetTitle(f"Events/{y_ax_val} GeV")
        plot.GetXaxis().SetLabelSize(0)    
        plot.GetXaxis().SetTitleSize(0)
        plot.GetYaxis().SetTitleOffset(0.9)

        if cs:
            max_y = max(cs[n]+(fetchError(q,cs[n])[1]-cs[n]) for n in range(len(cs)))
            plot.SetMaximum(2.0*max_y)
        else:
            plot.SetMaximum(10)
        plot.Draw()

        g1_top = None
        g2_top = None
        if loop==1 and show_bands:
            g2_top, g1_top = prefit_band
            g2_top.SetFillColor(ROOT.kYellow)
            g2_top.SetLineColor(ROOT.kBlack)
            g2_top.SetLineWidth(1)
            g1_top.SetFillColor(ROOT.kGreen)
            g1_top.SetLineColor(ROOT.kBlack)
            g1_top.SetLineWidth(1)
            g2_top.Draw("F")
            g1_top.Draw("F")
            plot.getCurve("bkg_curve").Draw("L same")
            plot.getCurve("sb_curve").Draw("L same")
            ROOT.gPad.RedrawAxis()

        #legend and contents to draw on the canvas
        leg = ROOT.TLegend(0.17, 0.5, 0.68, 0.88)
        cate = ROOT.TLatex()
        cate.SetTextSize(0.043)
        cate.DrawLatexNDC(0.59, 0.84, f"c#tau = {lifetime} mm, m_{{#phi}} = {mass} GeV")
        cat_label = {"asym": "asymmetric", "none": "combined", "all": "combined"}.get(cat, cat)
        cate.DrawLatexNDC(0.59, 0.78, f"category: {cat_label}")
        #cate.DrawLatexNDC(0.59, 0.72, f"order = {order}")
        cate.DrawLatexNDC(0.59, 0.72, fit_label)
        leg.AddEntry(gres1,"Background Data", "pe")
        leg.AddEntry("sb_curve","S+B fit sum", "L")
        leg.AddEntry("bkg_curve", "B component", "L")
        if show_bands and loop==0:
            leg.AddEntry("bkg_1sigma",r"\pm 1 \sigma", "f")
            leg.AddEntry("bkg_2sigma",r"\pm 2 \sigma", "f")
        elif show_bands:
            leg.AddEntry(g1_top,r"\pm 1 \sigma", "f")
            leg.AddEntry(g2_top,r"\pm 2 \sigma", "f")
        else:
            blank_entries = [ROOT.TLine(), ROOT.TLine()]
            for b in blank_entries:
                leg.AddEntry(b, " ", "")
        leg.SetHeader("H #rightarrow #phi#phi #rightarrow 4#gamma")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.Draw("SAME")
        cate.Draw("SAME")
        gres1.Draw("p, same")

        sb_rc = plot.getCurve("sb_curve")
        bkg_rc = plot.getCurve("bkg_curve")
        g_smb = ROOT.TGraph(sb_rc.GetN())
        for i in range(sb_rc.GetN()):
            xi = sb_rc.GetX()[i]
            g_smb.SetPoint(i, xi, sb_rc.GetY()[i] - bkg_rc.Eval(xi))
        g_smb.SetLineColor(2)
        g_smb.SetLineWidth(2)

        band_graphs = []
        if show_bands:
            if loop==0:
                band_curves = [(plot.getCurve("bkg_2sigma"), ROOT.kYellow), (plot.getCurve("bkg_1sigma"), ROOT.kGreen)]
            else:
                band_curves = [(g2_top, ROOT.kYellow), (g1_top, ROOT.kGreen)]
            for bc, color in band_curves:
                gb = ROOT.TGraph(bc.GetN())
                for i in range(bc.GetN()):
                    bxi = bc.GetX()[i]
                    gb.SetPoint(i, bxi, bc.GetY()[i] - bkg_rc.Eval(bxi))
                gb.SetFillColor(color)
                gb.SetLineColor(color)
                band_graphs.append(gb)

        #cd into bottom panel and begin populating it with s-b curve and datapoints
        pad2.cd()
        resid_hist = plot.residHist("data_points", "bkg_curve")
        nres  = resid_hist.GetN()
        xs = resid_hist.GetX()
        ys = resid_hist.GetY()
        x_vals = [xs[i] for i in range(nres)]
        y_vals = [ys[i] for i in range(nres)]

        #define the residual data points and set errors/locations manually
        g_res=ROOT.TGraphAsymmErrors()
        q = (1-0.6827)/2.0
        for n in range(nres):
            x_n=xs[n]
            y_n=ys[n]
            c_n = cs[n]
            g_res.SetPoint(n, x_n, y_n)
            error=fetchError(q, c_n)
            g_res.SetPointError(n, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        g_res.SetMarkerStyle(20)
        g_res.SetMarkerColor(ROOT.kBlack)
        g_res.SetLineColor(ROOT.kBlack)

        #fill axis titles on the bottom panel and write on the canvas
        x.setRange(bins[1], bins[2])
        x.setBins(bins[0])
        lower_plot = x.frame(ROOT.RooFit.Range(bins[1], bins[2]))
        lower_plot.GetXaxis().SetTitle("m_{#gamma#gamma#gamma#gamma} (GeV)")
        lower_plot.GetYaxis().SetTitle("")
        lower_plot.GetXaxis().SetLabelSize(0.11)
        lower_plot.GetYaxis().SetLabelSize(0.11)
        lower_plot.GetXaxis().SetTitleSize(0.14)
        lower_plot.GetYaxis().SetTitleSize(0.14)
        lower_plot.GetXaxis().SetTitleOffset(1.1)
        lower_plot.GetYaxis().SetTitleOffset(0)

        #dynamic y axis scaling
        n_points = g_res.GetN()
        if n_points > 0:
            ymin = float("inf")
            ymax = float("-inf")
            for i in range(n_points):
                y = g_res.GetY()[i]
                y_err_low = g_res.GetErrorYlow(i)
                y_err_high = g_res.GetErrorYhigh(i)
                ymin = min(ymin, y-y_err_low)
                ymax = max(ymax, y+y_err_high)

            for gb in band_graphs:
                for i in range(gb.GetN()):
                    gy = gb.GetY()[i]
                    ymin = min(ymin, gy)
                    ymax = max(ymax, gy)

            y_range = ymax-ymin
            margin = 0.3*y_range if y_range>0 else 1.0
            lower_plot.SetMinimum(ymin-margin)
            lower_plot.SetMaximum(ymax+margin)
        else:
            lower_plot.SetMinimum(-2)
            lower_plot.SetMaximum(2)

        n_bins =bins[0]
        x_min=bins[1]
        x_max=bins[2]

        dx=(x_max-x_min)/n_bins
        g_bzb = ROOT.TGraph(n_bins)

        for i in range(n_bins):
            x_center=x_min+(i+0.5)*dx
            g_bzb.SetPoint(i, x_center, 0.0)

        lower_plot.Draw("axis")
        for gb in band_graphs:
            gb.Draw("F")
        g_bzb.SetLineColor(ROOT.kRed)
        g_bzb.SetLineStyle(2)
        g_bzb.SetLineWidth(2)
        g_bzb.Draw("L")
        g_smb.Draw("L")
        g_res.Draw("P")
        ROOT.gPad.RedrawAxis()

        leg2 = ROOT.TLegend(0.15, 0.8, 0.68, 0.95)
        leg2.SetHeader("B component subtracted")
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.Draw("Same")
        pad1.cd()
        iPeriod = 4 if year in {"2017", "2018", "Run2"} else 5
        CMS_lumi.cmsTextSize = 0.85
        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.lumiTextRightOffset = 0.0
        CMS_lumi.CMS_lumi(pad1, iPeriod, 0, year, lumi[year], lumi_13TeV=f"{lumi[year]/1000:.1f}", extraText="Work in Progress")
        can.Update()
        can.cd()
        can.Update()
        if loop==0:
            can.SaveAs(f"postfit_{cat}_{year}_m{mass}_ct{lifetime}.png")
            can.SaveAs(f"postfit_{cat}_{year}_m{mass}_ct{lifetime}.pdf")
        else:
            can.SaveAs(f"prefit_{cat}_{year}_m{mass}_ct{lifetime}.png")
            can.SaveAs(f"prefit_{cat}_{year}_m{mass}_ct{lifetime}.pdf")
        loop+=1

