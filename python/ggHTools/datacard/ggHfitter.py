import ROOT
from ggHparameters import signal_window, fit_window, lower_sb, upper_sb, dcb_mean, dcb_sigma, dcb_alpha1, dcb_n1, dcb_alpha2, dcb_n2, bernstein_coeff, n_bins
ROOT.gROOT.SetBatch(False)


def _decasteljau_split(ctrl, t):
    pts = list(ctrl)
    left = [pts[0]]
    right = [pts[-1]]
    while len(pts) > 1:
        pts = [pts[i] * (1.0 - t) + pts[i + 1] * t for i in range(len(pts) - 1)]
        left.append(pts[0])
        right.insert(0, pts[-1])
    return left, right


def bernstein_subinterval(ctrl, u, v):
    _, right = _decasteljau_split(ctrl, u)
    s = (v - u) / (1.0 - u)
    left, _ = _decasteljau_split(right, s)
    return left


class Fitter(object):
    def __init__(self,poi = ['x']):
        self.cache=ROOT.TFile("cache.root","RECREATE")
        self.cache.cd()
        self.w=ROOT.RooWorkspace("w","w")
        self.dimensions = len(poi)
        self.poi=poi
        for v in poi:
            self.w.factory(f"{v}[{signal_window[0]},{signal_window[1]}]")

    def setRange(self,name,poi,low,high):
        self.w.var(poi).setRange(name,low,high)

    def factory(self,expr):
        self.w.factory(expr)

    def bernstein(self,name = 'model',poi='x',*,order):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        cList = ROOT.RooArgList()
        for i in range(0,order):
            self.w.factory(f"c_{i}[{bernstein_coeff[0]},{bernstein_coeff[1]}]")
            cList.add(self.w.var("c_"+str(i)))
        bernsteinPDF = ROOT.RooBernsteinFast(order)(name,name,self.w.var(poi),cList)
        getattr(self.w,'import')(bernsteinPDF)


    def doubleCB(self,name = 'model',poi='x'):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory(f"mean[{dcb_mean[0]},{dcb_mean[1]},{dcb_mean[2]}]")
        self.w.factory(f"sigma[{dcb_sigma[0]},{dcb_sigma[1]},{dcb_sigma[2]}]")
        self.w.factory(f"alpha1[{dcb_alpha1[0]},{dcb_alpha1[1]},{dcb_alpha1[2]}]")
        self.w.factory(f"n1[{dcb_n1[0]},{dcb_n1[1]},{dcb_n1[2]}]")
        self.w.factory(f"alpha2[{dcb_alpha2[0]},{dcb_alpha2[1]},{dcb_alpha2[2]}]")
        self.w.factory(f"n2[{dcb_n2[0]},{dcb_n2[1]},{dcb_n2[2]}]")
        doubleCB = ROOT.RooDoubleCB(name,name,self.w.var(poi),self.w.var("mean"),self.w.var("sigma"),self.w.var("alpha1"),self.w.var("n1"),self.w.var("alpha2"),self.w.var("n2"))
        getattr(self.w,'import')(doubleCB)

    def importBinnedData(self,histogram,poi = ["x"],name = "data"):
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print( 'Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead')
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))

    def fit(self, model="model", data="data", options=[], fitRange=None, sumW2=False):
        opts = list(options)
        opts.append(ROOT.RooFit.Save(True))
        opts.append(ROOT.RooFit.PrintLevel(-1))
        opts.append(ROOT.RooFit.Verbose(False))
        if fitRange:
            opts.append(ROOT.RooFit.Range(fitRange))
        if sumW2:
            opts.append(ROOT.RooFit.SumW2Error(True))
        return self.w.pdf(model).fitTo(self.w.data(data), *opts)

    def projection(self,model = "model",data="data",poi="x",filename="fit.root"):
        self.frame=self.w.var(poi).frame()
        self.w.data(data).plotOn(self.frame)
        self.w.pdf(model).plotOn(self.frame)
        self.c=ROOT.TCanvas("c","c")
        self.c.cd()
        self.frame.Draw()
        self.c.SaveAs(filename)
        return self.frame.chiSquare()

    def DCBandBernstein(self, order, dcbname="s_model", bname="b_model", poi="x", model_name="sb_model"):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("mean[125,120,130]")
        self.w.factory("sigma[2,0.1,20]")
        self.w.factory("alpha1[2,1,20]")
        self.w.factory("n1[2,1,50]")
        self.w.factory("alpha2[2,1,20]")
        self.w.factory("n2[2,1,50]")
        cList = ROOT.RooArgList()
        for i in range(0,order):
            self.w.factory(f"c_{i}[{bernstein_coeff[0]},{bernstein_coeff[1]}]")
            cList.add(self.w.var("c_"+str(i)))
        bernsteinPDF = ROOT.RooBernsteinFast(order)(bname,bname,self.w.var(poi),cList)
        doubleCB = ROOT.RooDoubleCB(dcbname,dcbname,self.w.var(poi),self.w.var("mean"),self.w.var("sigma"),self.w.var("alpha1"),self.w.var("n1"),self.w.var("alpha2"),self.w.var("n2"))
        getattr(self.w,'import')(bernsteinPDF)
        getattr(self.w,'import')(doubleCB)
        self.w.factory(f"SUM::{model_name}(nsig[0,0,10000]*{dcbname}, nbkg[100,0,100000]*{bname})")
        return


def fitBKG(file, hist, output_name, *, order, POI="mass", verbose=False):
    f = ROOT.TFile(file)
    bkg = f.Get(hist)

    fitter = Fitter([POI])
    fitter.importBinnedData(bkg, [POI])
    fitter.bernstein('model', POI, order=order)
    fitter.setRange("lower", POI, lower_sb[0], lower_sb[1])
    fitter.setRange("upper", POI, upper_sb[0], upper_sb[1])
    bkg_fit_result = fitter.fit("model", "data", fitRange="lower,upper")
    chi2 = fitter.projection("model", "data", POI, filename=output_name)

    x = fitter.w.var(POI)
    x.setRange("sr", signal_window[0], signal_window[1])
    x.setRange("sb_low", lower_sb[0], lower_sb[1])
    x.setRange("sb_high", upper_sb[0], upper_sb[1])
    nset = ROOT.RooArgSet(x)
    pdf = fitter.w.pdf("model")
    I_sr = pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sr")).getVal()
    I_sb = (pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sb_low")).getVal()
            + pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("sb_high")).getVal())
    ratio = I_sr / I_sb if I_sb > 0 else 0.0

    a_dom = x.getMin()
    b_dom = x.getMax()
    full_ctrl = [1.0] + [fitter.w.var(f"c_{i}").getVal() for i in range(order)]
    u = (signal_window[0] - a_dom) / (b_dom - a_dom)
    v = (signal_window[1] - a_dom) / (b_dom - a_dom)
    sr_ctrl = bernstein_subinterval(full_ctrl, u, v)
    sr_coeffs = [sr_ctrl[i + 1] / sr_ctrl[0] for i in range(order)]

    fitter_sr = Fitter([POI])
    fitter_sr.bernstein('model', POI, order=order)
    for i in range(order):
        fitter_sr.w.var(f"c_{i}").setVal(sr_coeffs[i])
    fitter_sr.w.factory(f"sr_sb_ratio[{ratio}]")

    have_band = False
    try:
        if bkg_fit_result is None or bkg_fit_result.covQual() < 2:
            raise RuntimeError("sideband fit covariance unusable (covQual < 2)")
        pdf.removeStringAttribute("fitrange")
        M_band = (bkg.Integral() * ratio) / I_sr if I_sr > 0 else 0.0
        band_frame = x.frame(ROOT.RooFit.Bins(n_bins), ROOT.RooFit.Range(signal_window[0], signal_window[1]))
        pdf.plotOn(band_frame, ROOT.RooFit.VisualizeError(bkg_fit_result, 1, ROOT.kFALSE), ROOT.RooFit.Normalization(M_band, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Name("e1"), ROOT.RooFit.DrawOption("F"))
        pdf.plotOn(band_frame, ROOT.RooFit.VisualizeError(bkg_fit_result, 2, ROOT.kFALSE), ROOT.RooFit.Normalization(M_band, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Name("e2"), ROOT.RooFit.DrawOption("F"))
        c1 = band_frame.getCurve("e1")
        c2 = band_frame.getCurve("e2")
        g_band1 = ROOT.TGraph(c1.GetN())
        g_band2 = ROOT.TGraph(c2.GetN())
        for i in range(c1.GetN()):
            g_band1.SetPoint(i, c1.GetX()[i], c1.GetY()[i])
        for i in range(c2.GetN()):
            g_band2.SetPoint(i, c2.GetX()[i], c2.GetY()[i])
        have_band = True
    except Exception:
        print("skipping prefit band: sideband fit covariance unusable (likely a degenerate fit)")

    output = ROOT.TFile(output_name, "UPDATE")
    output.cd()
    fitter_sr.w.Write("w")
    if have_band:
        g_band1.Write("band_1sigma")
        g_band2.Write("band_2sigma")
    output.Close()
    if verbose:
        print("bkg chi-squared={} sr_sb_ratio={}".format(chi2, ratio))
    return ratio

def fitSIG(file, hist, output_name, POI="mass", verbose=False):
    f = ROOT.TFile(file)
    signal = f.Get(hist)
    fitterS = Fitter([POI])
    fitterS.doubleCB('model',POI)
    fitterS.importBinnedData(signal,[POI],name = "data")
    #fit twice
    fitterS.setRange("signal_window", "mass", *signal_window)
    fitterS.fit("model","data", sumW2=True, fitRange="signal_window")
    chi2 = fitterS.projection("model","data",POI,filename=output_name)
    if verbose:
        print("signal chi-squared={}".format(chi2))
    output = ROOT.TFile(output_name, "UPDATE")
    output.cd()
    fitterS.w.Write("w")
    output.Close()
    return

def fitSIGBKG(file, sighist, bkghist,output_name, order, POI="mass"):
    f = ROOT.TFile(file)
    signal = f.Get(sighist)
    bkg = f.Get(bkghist)
    fitter = Fitter([POI])
    fitter.DCBandBernstein(order,dcbname="model_s",bname="model_b",poi=POI,model_name="model_sb")
    fitter.importBinnedData(signal,[POI], "data_s")
    fitter.importBinnedData(bkg,[POI], "data_b")
    fitter.setRange("signal_window", "mass", signal_window[0], signal_window[1])
    s_fit_result=fitter.fit("model_s","data_s", sumW2=True, fitRange="signal_window")
    print(f"[SIG fit] status={s_fit_result.status()} covQual={s_fit_result.covQual()} "
      f"sigma={fitter.w.var('sigma').getVal():.3f} +/- {fitter.w.var('sigma').getError():.3f} "
      f"mean={fitter.w.var('mean').getVal():.3f}")
    fitter.setRange("lower", "mass", lower_sb[0], lower_sb[1])
    fitter.setRange("upper", "mass", upper_sb[0], upper_sb[1])
    b_fit_result=fitter.fit("model_b","data_b",fitRange="lower,upper", sumW2=False)
    #chi2S = fitter.projection("model_s","data_s",POI,filename=output_name)
    #chi2B = fitter.projection("model_b","data_b",POI,filename=output_name)
    output = ROOT.TFile(output_name, "UPDATE")
    output.cd()
    fitter.w.Write("w")
    output.Close()
    return s_fit_result, b_fit_result
