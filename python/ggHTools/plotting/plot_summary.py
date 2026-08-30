import ROOT
from datacard import ggHfitter
from datacard.ggHdatacardmaker import recommended_photon_id
from plotting.style import tdrstyle
from plotting.style import CMS_lumi
from ggHparameters import (lumi, signal_path, bkg_path, signal_window, order_fit, lower_sb, upper_sb,summary_bin_width, year_config)
import ggHcuts as cuts
from plotting.plottingtools import fetchError, save_histos
import os
ROOT.gROOT.SetBatch(True)

sb_bins=[int(round((upper_sb[1]-lower_sb[0])/summary_bin_width)), lower_sb[0], upper_sb[1]]
residual_scale=1.3
pad2_top=0.25*residual_scale
pad1_bottom=pad2_top+0.05
data_cache={}


def sum_gen_weights(path):
    with ROOT.TFile.Open(path) as root_file:
        return sum(float(entry.genEventSumw) for entry in root_file.Get("Runs"))


def data_histogram(year, mass, config, bins, photon_id):
    key=(year, str(mass), tuple(bins), photon_id)
    if key not in data_cache:
        data_cut=cuts.background_selection(mass, photon_id)
        total=None
        for y in config["data_years"]:
            data_df=ROOT.RDataFrame("ggH4g", bkg_path(y)).Filter(data_cut)
            proxy=data_df.Histo1D((f"data_tmp_{y}_m{mass}", f"data_hist;4#gamma mass;Events",
                                   bins[0], bins[1], bins[2]), f"best_4g_corr_mass_m{mass}")
            histogram=proxy.GetValue().Clone(f"data_{y}_m{mass}")
            histogram.SetDirectory(0)
            if total is None:
                total=histogram.Clone("data_hist")
                total.SetDirectory(0)
            else:
                total.Add(histogram)
        data_cache[key]=total
    clone=data_cache[key].Clone("data_hist")
    clone.SetDirectory(0)
    return clone


def run(mass, ctau, year, output_dir, work_dir=None, formats=("png","pdf")):
    config = year_config[year]
    photon_id = recommended_photon_id(year)
    if work_dir is None:
        work_dir = os.path.join(output_dir, "fit_workspaces")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText="Work in Progress"
    CMS_lumi.lumi_13TeV=f"{year}, {lumi[year]/1000} fb^{-1}"
    ROOT.gStyle.SetEndErrorSize(2)

    
    bins_sig = list(sb_bins)
    bins_data = list(sb_bins)
    bins=bins_data
    var=f"best_4g_corr_mass_m{mass}"

    years_to_process = config["signal_years"]

    h_data=data_histogram(year, mass, config, bins_data, photon_id)

    #signal histogram built separately for each year.
    h_sig_total = None
    sig_cut=cuts.signal_selection(mass, photon_id)
    for y in years_to_process:
        sig_y = signal_path(mass, ctau, y)
        sumw_y = sum_gen_weights(sig_y)
        signal_df_y = ROOT.RDataFrame("ggH4g", sig_y)
        weight_y = cuts.mc_weight(sumw_y)
        signal_df_y = signal_df_y.Define("event_weight", weight_y)
        signal_df_y = signal_df_y.Filter(sig_cut)
        h_y = signal_df_y.Histo1D((f"sig_tmp_{y}", f"sig_tmp_{y}", bins_sig[0], bins_sig[1], bins_sig[2]), var, "event_weight")
        h_y_clone = h_y.GetValue().Clone(f"sig_scaled_{y}")
        h_y_clone.Scale(lumi[y])
        if h_sig_total is None:
            h_sig_total = h_y_clone.Clone("signal_hist")
        else:
            h_sig_total.Add(h_y_clone)
    h_sig = h_sig_total

    id_suffix = f"_{photon_id}"
    hist_file = save_histos(mass, year, ctau, h_data, h_sig, tag=photon_id, directory=work_dir)
    fit_file = os.path.join(work_dir, f"SB_fit_result_m{mass}_ct{ctau}_year{year}{id_suffix}.root")
    sresult, bresult=ggHfitter.fitSIGBKG(hist_file, "signal_hist", "data_hist", fit_file, order=order_fit)
    SB_file=ROOT.TFile.Open(fit_file)
    w = SB_file.Get("w")
    x= w.var("mass")
    s_model = w.pdf("model_s")
    b_model = w.pdf("model_b")

    s_model.removeStringAttribute("fitrange")
    b_model.removeStringAttribute("fitrange")

    bkg_norm = h_data.Integral()
    sig_norm = h_sig.Integral()

    x.setRange("sig_window", signal_window[0], signal_window[1])
    sig_window_frac = s_model.createIntegral(ROOT.RooArgSet(x), ROOT.RooFit.NormSet(ROOT.RooArgSet(x)), ROOT.RooFit.Range("sig_window")).getVal()
    print(f"SIGNAL YIELD total: {sig_norm:.4f}")
    print(f"SIGNAL YIELD in [{signal_window[0]},{signal_window[1]}]: {sig_norm * sig_window_frac:.4f}")

    x.setRange("sb_low", lower_sb[0], lower_sb[1])
    x.setRange("sb_high", upper_sb[0], upper_sb[1])
    nset = ROOT.RooArgSet(x)
    b_sr = b_model.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sig_window")).getVal()
    b_sb = (b_model.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sb_low")).getVal()
            + b_model.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sb_high")).getVal())
    bkg_sr_ratio = b_sr / b_sb
    bkg_norm_full = bkg_norm / b_sb
    print(f"BKG YIELD sideband data: {bkg_norm:.4f}")
    print(f"BKG YIELD in [{signal_window[0]},{signal_window[1]}]: {bkg_norm * bkg_sr_ratio:.4f}")

    n_sig = ROOT.RooRealVar("n_sig", "n_sig", sig_norm)
    n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", bkg_norm_full)
    sb_model = ROOT.RooAddPdf("model_sb", "S+B model",
        ROOT.RooArgList(s_model, b_model),
        ROOT.RooArgList(n_sig, n_bkg))

    #set the proportions and make it beautiful
    can = ROOT.TCanvas("c")
    pad1 = ROOT.TPad("pad1", "pad1", 0,   pad1_bottom, 1, 1.0)
    pad1.SetTopMargin(0.08506945)
    pad1.SetBottomMargin(0.00)  
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, pad2_top)
    pad2.SetTopMargin(0.05/residual_scale)
    pad2.SetBottomMargin(0.35/residual_scale)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()
    pad1.cd()
    x.setBins(bins[0])
    plot=x.frame()

    #create data_obs as a TH1 obj so i can look through the contents and throw away the error bars (roofit fucks them up. i compute them manually later for the case of N=0 bins):
    data_obs = w.data("data_b")
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
    n_point = 0
    for n in range(cloned_data.GetNbinsX()):
        b_n = bin_centers[n]
        if signal_window[0] <= b_n <= signal_window[1]:
            continue
        c_n = cs[n]
        error=fetchError(q, c_n)
        gres1.SetPoint(n_point, b_n, c_n)
        gres1.SetPointError(n_point, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        n_point += 1

    gres1.SetLineColor(ROOT.kBlack)
    gres1.SetMarkerColor(ROOT.kBlack)
    gres1.SetMarkerStyle(20)

    #plotting. the order here matters a lot in order to get the brazil plot colors to show up correctly and to get the data points on top of everything
    cloned_data_binned.plotOn(plot,ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"),ROOT.RooFit.XErrorSize(0),ROOT.RooFit.Invisible())

    show_bands = data_obs_TH1.Integral()!=0
    if show_bands:
        print("data obs integral: ", data_obs_TH1.Integral())
        b_model.plotOn(plot,ROOT.RooFit.VisualizeError(bresult, 2, ROOT.kFALSE),ROOT.RooFit.Normalization(bkg_norm_full, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.FillColor(ROOT.kYellow),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_2sigma"),ROOT.RooFit.DrawOption("F"))
        b_model.plotOn(plot,ROOT.RooFit.VisualizeError(bresult, 1, ROOT.kFALSE),ROOT.RooFit.Normalization(bkg_norm_full, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.FillColor(ROOT.kGreen),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_1sigma"),ROOT.RooFit.DrawOption("F"))
        sb_norm_to_use = bkg_norm_full + sig_norm
    else:
        print("data is 0. ignoring uncertainty bands because uncertainties on fit parameters are unstable")
        sb_norm_to_use = bkg_norm_full + sig_norm

    b_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Name("bkg_curve"),ROOT.RooFit.Normalization(bkg_norm_full, ROOT.RooAbsReal.NumEvent))
    sb_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sb_curve"),ROOT.RooFit.Normalization(sb_norm_to_use, ROOT.RooAbsReal.NumEvent))

    integral_sb = sb_model.getVal(ROOT.RooArgSet(x))
    print("norm after scaling sb: ", integral_sb)
    integral_b = b_model.getVal(ROOT.RooArgSet(x))
    print("norm after scaling b: ", integral_b)

    y_ax_val=round((bins[2]-bins[1])/(bins[0]),2)
    plot.GetYaxis().SetTitle(f"Events/{y_ax_val} GeV")
    plot.GetXaxis().SetLabelSize(0)
    plot.GetXaxis().SetTitleSize(0)
    plot.GetYaxis().SetTitleOffset(0.7)

    #scale the frame to whichever is taller: the highest data error bar or the S+B curve peak,
    #so a narrow signal peak cannot run into the legend
    data_max = max((fetchError(q, c_n)[1] for c_n in cs), default=0.0)
    range_curve = plot.getCurve("sb_curve")
    curve_max = max(range_curve.GetY()[i] for i in range(range_curve.GetN()))
    max_y = max(data_max, curve_max)
    plot.SetMaximum(2.0*max_y if max_y > 0 else 10)
    plot.Draw()

    #legend and contents to draw on the canvas
    leg = ROOT.TLegend(0.2, 0.5, 0.68, 0.88)
    cate = ROOT.TLatex()
    cate.SetTextSize(0.06)
    cate.DrawLatexNDC(0.55, 0.84, rf"c#tau = {ctau} mm, m_{{#phi}} = {mass} GeV")
    leg.AddEntry(gres1,"Blinded Data", "pe")
    leg.AddEntry("sb_curve","S+B fit sum", "L")
    leg.AddEntry("bkg_curve", "B component", "L")
    if show_bands:
        leg.AddEntry("bkg_1sigma",r"\pm 1 \sigma", "f")
        leg.AddEntry("bkg_2sigma",r"\pm 2 \sigma", "f")
    else:
        blank_entries = [ROOT.TLine(), ROOT.TLine()]
        for b in blank_entries:
            leg.AddEntry(b, " ", "")
    leg.SetHeader("H #rightarrow #phi#phi #rightarrow 4#gamma")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("SAME")
    gres1.Draw("p, same")

    #build the signal-only curve for the bottom panel by subtracting the B curve from the S+B curve on the upper panel
    #this guarantees the bottom panel signal shape matches what's shown on top
    sb_curve_obj = plot.getCurve("sb_curve")
    b_curve_obj = plot.getCurve("bkg_curve")

    n_curve_points = sb_curve_obj.GetN()
    g_sig_only = ROOT.TGraph(n_curve_points)
    for i in range(n_curve_points):
        x_i = sb_curve_obj.GetX()[i]
        y_sb = sb_curve_obj.GetY()[i]
        y_b = b_curve_obj.interpolate(x_i)
        g_sig_only.SetPoint(i, x_i, y_sb - y_b)
    g_sig_only.SetLineColor(ROOT.kRed)
    g_sig_only.SetLineWidth(2)

    band_graphs = []
    if show_bands:
        for cname, color in (("bkg_2sigma", ROOT.kYellow), ("bkg_1sigma", ROOT.kGreen)):
            bc = plot.getCurve(cname)
            gb = ROOT.TGraph(bc.GetN())
            for i in range(bc.GetN()):
                bxi = bc.GetX()[i]
                gb.SetPoint(i, bxi, bc.GetY()[i] - b_curve_obj.interpolate(bxi))
            gb.SetFillColor(color)
            gb.SetLineColor(color)
            band_graphs.append(gb)

    #cd into bottom panel and begin populating it with s-b curve and datapoints
    pad2.cd()
    resid_hist = plot.residHist("data_points", "bkg_curve")
    nres  = resid_hist.GetN()
    xs = resid_hist.GetX()
    ys = resid_hist.GetY()

    #define the residual data points and set errors/locations manually
    g_res=ROOT.TGraphAsymmErrors()
    q = (1-0.6827)/2.0
    n_res_point = 0
    for n in range(nres):
        x_n=xs[n]
        y_n=ys[n]
        if signal_window[0] <= x_n <= signal_window[1]:
            continue
        c_n = cs[n]
        g_res.SetPoint(n_res_point, x_n, y_n)
        error=fetchError(q, c_n)
        g_res.SetPointError(n_res_point, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        n_res_point += 1
    g_res.SetMarkerStyle(20)
    g_res.SetMarkerColor(ROOT.kBlack)
    g_res.SetLineColor(ROOT.kBlack)

    #fill axis titles on the bottom panel and write on the canvas
    x.setRange(bins[1], bins[2])
    x.setBins(bins[0])
    lower_plot = x.frame(ROOT.RooFit.Range(bins[1], bins[2]))
    lower_plot.GetXaxis().SetTitle("m_{4#gamma} (GeV)")
    lower_plot.GetYaxis().SetTitle("")
    lower_plot.GetXaxis().SetLabelSize(0.09/residual_scale)
    lower_plot.GetYaxis().SetLabelSize(0.07/residual_scale)
    lower_plot.GetXaxis().SetTitleSize(0.15/residual_scale)
    lower_plot.GetYaxis().SetTitleSize(0.10/residual_scale)
    lower_plot.GetXaxis().SetTitleOffset(0.9)
    lower_plot.GetYaxis().SetTitleOffset(0.5)
    lower_plot.GetYaxis().SetNdivisions(505)

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

        #the b-subtracted signal curve is drawn on this pad, so it has to set the range too
        for i in range(g_sig_only.GetN()):
            gy = g_sig_only.GetY()[i]
            ymin = min(ymin, gy)
            ymax = max(ymax, gy)

        span = max(abs(ymin), abs(ymax), 1.0)
        lower_plot.SetMaximum(1.2*ymax if ymax > 0 else 0.2*span)
        lower_plot.SetMinimum(1.2*ymin if ymin < 0 else -0.2*span)
    else:
        lower_plot.SetMinimum(-2)
        lower_plot.SetMaximum(2)

    lower_plot.Draw("axis")
    lower_plot.GetYaxis().SetTitle("")
    lower_plot.SetTitle("")

    n_bins =bins[0]
    x_min=bins[1]
    x_max=bins[2]

    dx=(x_max-x_min)/n_bins
    g_bzb = ROOT.TGraph(n_bins)

    for i in range(n_bins):
        x_center=x_min+(i+0.5)*dx
        g_bzb.SetPoint(i, x_center, 0.0)

    for gb in band_graphs:
        gb.Draw("F")
    g_bzb.SetLineColor(ROOT.kRed)
    g_bzb.SetLineStyle(2)
    g_bzb.SetLineWidth(2)
    g_bzb.Draw("L")
    g_sig_only.Draw("L same")
    g_res.Draw("P")
    ROOT.gPad.RedrawAxis()

    leg2 = ROOT.TLegend(0.63, 0.8, 1.03, 0.95)
    leg2.SetHeader("B component subtracted")
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.Draw("Same")
    pad1.cd()
    iPeriod = config["period"]
    CMS_lumi.cmsTextSize = 0.85
    CMS_lumi.lumiTextSize = 0.6
    CMS_lumi.lumiTextRightOffset = 0.0
    CMS_lumi.CMS_lumi(pad1, iPeriod, 0, year, lumi[year], lumi_13TeV=f"{lumi[year]/1000:.2f}", extraText="Work in Progress")
    can.Update()
    can.cd()
    can.Update()
    basename = f"summary_plot_m{mass}_ct{ctau}_{year}{id_suffix}"
    outputs = []
    for extension in formats:
        output = os.path.join(output_dir, f"{basename}.{extension}")
        can.SaveAs(output)
        outputs.append(output)

    can.Close()
    SB_file.Close()
    del can, pad1, pad2, plot, lower_plot
    del w, x, s_model, b_model, sb_model
    del bresult, sresult
    ROOT.gROOT.GetListOfCanvases().Clear()
    ROOT.gROOT.GetListOfFiles().Clear()
    return outputs
