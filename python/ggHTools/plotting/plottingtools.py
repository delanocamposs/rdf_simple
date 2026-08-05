import ROOT


def fetchError(q, n):
    l=0
    if n!=0:
        l = ROOT.Math.chisquared_quantile_c(1.0-q, 2.0*n)/2.0
    u=ROOT.Math.chisquared_quantile_c(q, 2.0*n+2)/2.0
    return [l,u]


def getPoisson(h):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        gRate.SetPoint(n, thresh, N)
        error = fetchError(q, N)
        gRate.SetPointError(n, 0, 0, (N-error[0]), (error[1]-N))
        n+=1
    return gRate


def getPoisson2(h, scale):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        #if N == 0:
        #    continue
        gRate.SetPoint(n, thresh, scale*N)
        error = fetchError(q, N)
        gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, scale*(N-error[0]), scale*(error[1]-N))
        n+=1
    return gRate


def save_histos(mass, year, ctau, hist1, hist2):
    outfile=ROOT.TFile(f"sig_bkg_summary_histos_m{mass}_ct{ctau}_year{year}.root", "RECREATE")
    hist1.GetValue().Write()
    hist2.GetValue().Write()
    outfile.Close()
    return
