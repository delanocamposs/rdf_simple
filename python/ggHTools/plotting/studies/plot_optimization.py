import ROOT 
import subprocess 
from datacard.ggHdatacardmaker import main
from datacard.ggHdatacardmaker import cleanup

def two_dim(xmin, xmax,ymin,  ymax, xbins, ybins, order, year):
    ##generates the data points for a beautiful heat map to optimize categroeis

    ROOT.gROOT.SetBatch(True)
    limit_array=[]
    xbinw=int(xmax-xmin/xbins)
    ybinw=int(ymax-ymin/ybins)
    for i in range(xmin, xbins+1):
        for j in range(ymin, ybins+1):
            main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowlow"],period=f"{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}", bins=[90, 110, 140], order=order, lxy1=i*xbinw, lxy2=j*ybinw)
            main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowhigh"],period=f"{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}", bins=[90,110,140], order=order, lxy1=i*xbinw, lxy2=j*ybinw)
            main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highlow"],period=f"{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}", bins=[90,110,140], order=order, lxy1=i*xbinw, lxy2=j*ybinw)
            main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highhigh"],period=f"{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}", bins=[90, 110, 140], order=order, lxy1=i*xbinw, lxy2=j*ybinw)

            with open(f"datacard_combined_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt", "w") as f:
                subprocess.run(["combineCards.py", f"datacard_ggH_4g_lowlow_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt", f"datacard_ggH_4g_lowhigh_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt", f"datacard_ggH_4g_highlow_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt", f"datacard_ggH_4g_highhigh_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt"], stdout=f)

            subprocess.run(["text2workspace.py", f"datacard_combined_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.txt", "-o", f"datacard_combined_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.root"])
            subprocess.run(["combine", "-M", "AsymptoticLimits", f"datacard_combined_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.root", "-m", "125"])
            subprocess.run(["mv", f"higgsCombineTest.AsymptoticLimits.mH125.root",f"higgsCombineTest_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.AsymptoticLimits.mH125.root" ])
            file=ROOT.TFile.Open(f"higgsCombineTest_{year}_lxy1_{i*xbinw}_lxy2_{j*ybinw}.AsymptoticLimits.mH125.root")
            limit_tree=file.Get("limit")
            arr=[]
            for entry in limit_tree:
                upper_lim = entry.limit
                arr.append(upper_lim)
            print(arr[2])
            limit_array.append(arr[2])

    with open("asymptotic_limits_vals.txt", "w") as f:
        for l in limit_array:
            f.write(f"{l} \n")

def one_dim(xmin, xmax, xinterval, order, year):
    #generates 1D data for cat optimization
    #xinterval is the step between datapoints in the x value. ie, datapoints every 10 cm step xinterval=10

    ROOT.gROOT.SetBatch(True)
    limit_array=[]
    #xbinw=int((xmax-xmin)/xbins)
    npoints=int(((xmax-xmin)/(xinterval)))+1 ##+1 because we need data at the endpoints
    for i in range(0, npoints):
        main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowlow"],period=f"{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}", bins=[90, 110, 140], order=order, lxy1=i*xinterval+xmin, lxy2=i*xinterval+xmin)
        main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowhigh"],period=f"{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}",  bins=[90, 110, 140],order=order, lxy1=i*xinterval+xmin, lxy2=i*xinterval+xmin)
        main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highlow"],period=f"{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}",  bins=[90, 110, 140],order=order, lxy1=i*xinterval+xmin, lxy2=i*xinterval+xmin)
        main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highhigh"],period=f"{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}", bins=[90, 110, 140], order=order, lxy1=i*xinterval+xmin, lxy2=i*xinterval+xmin)

        with open(f"datacard_combined_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt", "w") as f:
            subprocess.run(["combineCards.py", f"datacard_ggH_4g_lowlow_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt", f"datacard_ggH_4g_lowhigh_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt", f"datacard_ggH_4g_highlow_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt", f"datacard_ggH_4g_highhigh_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt"], stdout=f)

        subprocess.run(["text2workspace.py", f"datacard_combined_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.txt", "-o", f"datacard_combined_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.root"])
        subprocess.run(["combine", "-M", "AsymptoticLimits", f"datacard_combined_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.root", "-m", "125"])
        subprocess.run(["mv", f"higgsCombineTest.AsymptoticLimits.mH125.root",f"higgsCombineTest_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.AsymptoticLimits.mH125.root" ])
        file=ROOT.TFile.Open(f"higgsCombineTest_{year}_lxy1_{i*xinterval+xmin}_lxy2_{i*xinterval+xmin}.AsymptoticLimits.mH125.root")
        limit_tree=file.Get("limit")
        arr=[]
        for entry in limit_tree:
            upper_lim = entry.limit
            arr.append(upper_lim)
        print(arr[2])
        limit_array.append(arr[2])

    with open("asymptotic_limits_vals.txt", "w") as f:
        for l in limit_array:
            f.write(f"{l} \n")



if __name__=="__main__":
    one_dim(-20, 100,10, 3, "Run-2")

