import ROOT
import sys
from array import array
import numpy as np
import pickle
import ctypes
import math
import common.tdrstyle as tdrstyle
import common.CMS_lumi as CMS_lumi
from python.VHTools.py_helpers import *
sty = tdrstyle.setTDRStyle()
from scipy.stats import chi2,poisson
import matplotlib.pyplot as plt
import mplhep as mh

drawprelim = True

@contextmanager
def open_root_file(filename, mode="READ"):
    """
    Context manager for safely opening and closing a ROOT TFile.

    This ensures that the ROOT file is automatically closed when the
    context exits, even in the event of an exception.

    Parameters:
    ----------
    filename : str
        Path to the ROOT file to open.

    mode : str, optional
        Mode in which to open the file (default is "READ").
        Typical options include:
        - "READ"   : Open existing file for reading
        - "UPDATE" : Open existing file for reading/writing
        - "RECREATE": Create new file, overwriting if exists
        - "NEW"    : Create new file, fail if exists

    Yields:
    ------
    ROOT.TFile
        The opened ROOT file object.

    Example:
    --------
    >>> with open_root_file("data.root") as f:
    ...     h = f.Get("myHistogram")
    ...     h.Draw()
    """
    file = ROOT.TFile.Open(filename, mode)
    try:
        yield file
    finally:
        if file and file.IsOpen():
            file.Close()

class plotter_base(object):

    def __init__(self):
        self.fillstyle=1001
        self.linestyle=1
        self.linecolor=1
        self.linewidth=2
        self.fillcolor=ROOT.kOrange-3
        self.markerstyle=20
        self.corrFactors=[]
        self.n_bootstraps=10000

    def addCorrectionFactor(self,value,model):
        corr=dict()
        corr['value']=value
        corr['model']=model
        self.corrFactors.append(corr)

    def setLineProperties(self,linestyle,linecolor,linewidth):
        self.linestyle=linestyle
        self.linecolor=linecolor
        self.linewidth=linewidth 

    def setFillProperties(self,fillstyle,fillcolor):
        self.fillstyle=fillstyle
        self.fillcolor=fillcolor

    def setMarkerProperties(self,markerstyle):
        self.markerstyle=markerstyle
        
    def poisson_confidence_interval(self,k, confidence=0.6827):
        """
        Computes the Garwood (exact) confidence interval for a Poisson rate.
        k: observed number of events (can be 0)
        confidence: confidence level (e.g., 0.95)
        """
        k_array = k
        alpha = 1 - confidence
    
        # Lower bound: use 0 where k is 0, otherwise calculate chi2
        # We use np.maximum(2*k, 0) to ensure we don't pass negatives
        lower = np.where(
            k_array > 0, 
            chi2.ppf(alpha / 2, 2 * k_array) / 2, 
            0.0
        )
    
        # Upper bound: always calculated using 2*k + 2
        upper = chi2.ppf(1 - alpha / 2, 2 * k_array + 2) / 2
        return lower, upper


    def poisson_bootstrap_ci(self,weights, n_bootstraps=1000, ci_level=0.6827):
        """
        Calculates the CI for a sum of weighted Poisson events.
        """
        # Convert to array for vectorization
        n_events = len(weights)
    
        # Generate a matrix of Poisson(1) multipliers
        # Shape: (n_bootstraps, n_events)
        poisson_counts = np.random.poisson(1, size=(n_bootstraps, n_events))
    
        # Calculate the weighted sum for every bootstrap iteration
        # Matrix multiplication: (n_bootstraps, n_events) @ (n_events, 1)
        bootstrap_sums = poisson_counts @ weights
    
        # Calculate percentiles for the confidence interval
        lower_bound = (1 - ci_level) / 2
        upper_bound = 1 - lower_bound
    
        ci_low, ci_high = np.percentile(bootstrap_sums, [lower_bound * 100, upper_bound * 100])
        del bootstrap_sums
        del poisson_counts
        return ci_low, ci_high
    
    def array1d(self,var,cuts,model,include_overflow=False,error_mode='w2'):
        hist=self.hist1d(var,cuts,model,titlex = "",units = "")
        axis=hist.GetXaxis()
        num_bins = hist.GetNbinsX()+2
        data = np.ndarray(num_bins, dtype=np.float64, buffer=hist.GetArray())
        edges = np.array([axis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+1)])
        edges=np.append(edges,axis.GetBinUpEdge(num_bins-2))
        if error_mode=='w2':
            w2   = np.ndarray(num_bins, dtype=np.float64, buffer=hist.GetSumw2().GetArray())
            w2=np.array([w2,w2])            
        elif error_mode=='poisson':
            low, high = self.poisson_confidence_interval(data,0.6827)
            w2=np.array([np.square(data-low),np.square(high-data)])
        elif error_mode=='poisson_bootstrap':
            low = data.copy()
            high = data.copy()
            for i in range(0,len(data)):
                if (data[i]==0.0):
                    continue;
                if i==0:
                    low_edge=edges[0]
                    c = "*".join([cuts,f"({var}<{low_edge})"])
                elif i==(len(data)-1):
                    up_edge=edges[-1]
                    c = "*".join([cuts,f"({var}>={up_edge})"])
                else:
                    lo=edges[i-1]
                    hi=edges[i]
                    c = "*".join([cuts,f"({var}>={lo}&&{var}<{hi})"])
                weights = self.event_weights(c)
                if len(weights)>0:
#                    print("Debug, poisson bootstrap weights",weights)
                    l,h = self.poisson_bootstrap_ci(weights, n_bootstraps=self.n_bootstraps,ci_level=0.6827)
                    weights=None
                    low[i]=l
                    high[i]=h
                    d=data[i]
#                    print(f"Poisson Bootstrap: Observation = {d} Interval E [{l},{h}]")
            w2=np.array([np.square(data-low),np.square(high-data)])            
        del hist    
        if include_overflow==False:
            return edges.copy(),data[1:-1].copy(),w2[:,1:-1].copy()
        else:
            return edges.copy(),data.copy(),w2.copy()

    def array2d(self,var1,var2,cuts,model,include_overflow=False,error_mode='w2'):
        h=self.hist2d(var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy="")        
        nx = h.GetNbinsX() + 2
        ny = h.GetNbinsY() + 2
        data1d = np.ndarray(shape=(nx*ny,), buffer=h.GetArray(), dtype=np.float64)
        data=data1d.reshape(ny,nx)
        xaxis=h.GetXaxis()
        yaxis=h.GetYaxis()
        xedges = np.array([xaxis.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 1)])
        xedges=np.append(xedges,xaxis.GetBinUpEdge(h.GetNbinsX()))                     
        yedges = np.array([yaxis.GetBinLowEdge(i) for i in range(1, h.GetNbinsY() + 1)])
        yedges=np.append(yedges,yaxis.GetBinUpEdge(h.GetNbinsY()))

        if error_mode=='w2':
            w22d   = np.ndarray((nx*ny,), buffer=h.GetSumw2().GetArray(), dtype=np.float64)            
            w2=w22d.reshape(ny,nx)
            w2=np.array([w2,w2])
        elif error_mode=='poisson':
            low, high = self.poisson_confidence_interval(data,0.6827)
            w2=np.array([np.square(data-low),np.square(high-data)])
        elif error_mode=='poisson_bootstrap':
            cx=None
            cy=None            
            low = data.copy()
            high = data.copy()
            for i in range(0,len(data[0,:])):
                if i==0:
                    low_edge=xedges[0]
                    cx = f"({var1}<{low_edge})"
                elif i==(len(data[0,:])-1):
                    up_edge=xedges[-1]
                    cx = f"({var1}>={up_edge})"
                else:
                    lo=xedges[i-1]
                    hi=xedges[i]
                    cx = f"({var1}>={lo}&&{var1}<{hi})"
                for j in range(0,len(data[:,0])):
                    if data[j,i]==0.0:
                        continue
                    if j==0:
                        low_edge=yedges[0]
                        cy = f"({var2}<{low_edge})"
                    elif j==(len(data[:,0])-1):
                        up_edge=yedges[-1]
                        cy = f"({var2}>={up_edge})"
                    else:
                        lo=yedges[j-1]
                        hi=yedges[j]
                        cy = f"({var2}>={lo}&&{var2}<{hi})"
#                    print(f"Poisson bootstrap for ({i},{j})","*".join([cx,cy]))
                    weights = self.event_weights("*".join([cuts,cx,cy]))
                    if len(weights)>0:
                        l,h = self.poisson_bootstrap_ci(weights, n_bootstraps=self.n_bootstraps,ci_level=0.6827)
                        low[j,i]=l
                        high[j,i]=h
                        d=data[j,i]
#                        print(f"Observation = {d} Interval E [{l},{h}]")
                    weights=None
            w2=np.array([np.square(data-low),np.square(high-data)])
        del h    
        if include_overflow==False:
            return xedges.copy(),yedges.copy(),data[1:-1,1:-1].copy(),w2[:,1:-1,1:-1].copy()
        else:
            return xedges.copy(),yedges.copy(),data.copy(),w2.copy()        




class rdf_plotter(plotter_base):
    def __init__(self,file,isMC=False,weight = "1.0",tree='Events',defaultCuts = "1.0",report='report'):
        self.rdf = ROOT.RDataFrame(tree,file)
        self.weight=weight
        self.isMC=isMC
        self.defaultCuts = defaultCuts
        self.file = file
        self.report = report
        super(rdf_plotter,self).__init__()
        #if MC read the sum of weights and weigh the event
        if isMC:
            with open_root_file(file) as f:
                t = f.Get('Runs')
                sumw = 0.0
                for event in t:
                    sumw += event.genEventSumw
                self.weight = self.weight+'*(genWeight/{})*(sample_sigma)*(Pileup_weight)'.format(sumw if sumw>0 else 1.0)

    def readReport(self) -> dict[str]:
        out = {}
        with open_root_file(self.file) as f:
            report = f.Get(self.report)
            if not report:
                raise RuntimeError(f"Report {self.report} not found in {self.file}")
            
            cuts = [b.GetName() for b in report.GetListOfBranches()]
            out = {c: 0 for c in cuts}
            for event in report:
                for c in cuts:
                    out[c] += getattr(event, c, 0)

            return out
        

    def define(self, var, definition):
        if var in self.rdf.GetColumnNames():
            self.rdf = self.rdf.Redefine(var, definition)           
        else:    
            self.rdf = self.rdf.Define(var, definition)

    def redefine(self, var, definition):
        self.rdf = self.rdf.Redefine(var, definition)
    
    def display(self,var):
        self.rdf.Display(var).Print()

    def filter(self,condition):
        '''
        Inputs
        ======
        condition(string): condition to be used for filtering the dataframe that evaluates to True (1) or False (0)
        new_var(string): name of the new column  
        
        Returns
        =======
        RDataFrame with filter applied
        '''
        self.rdf = self.rdf.Filter(condition)

    def event_weights(self,cuts):
        c="1.0"
        for corr in self.corrFactors:
            c = c+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+c+")"
        rdf=self.rdf.Define('plot_weight',c)
        filtered = rdf.Filter("plot_weight!=0")
        arr=filtered.AsNumpy(["plot_weight"])
        del rdf
        del filtered        
        return arr['plot_weight'].copy()
        
        
    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo1D(model,var,'plot_weight')
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        for n in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(n) < 0:
                h.SetBinContent(n, 0)
        return h.GetPtr()



    
    def profile1d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley = "",unitsy = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Profile1D(model,var1,var2,'plot_weight')
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h.GetPtr()


    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')

        #h.Sumw2(0)
        #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h.GetPtr()

    def unrolled2d(self,var1,var2,cuts,model):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')

        hu=ROOT.TH1D(h.GetName(),h.GetTitle(),h.GetNbinsX()*h.GetNbinsY(),0,h.GetNbinsX()*h.GetNbinsY())
        for i in range(1,h.GetNbinsX()+1):
            for j in range(1,h.GetNbinsY()+1):
                b=h.GetBin(i,j)
                bin1d=(i-1)*h.GetNbinsX()+j
                hu.SetBinContent(bin1d,h.GetBinContent(b))
                hu.SetBinError(bin1d,h.GetBinError(b))
            
        hu.SetLineStyle(self.linestyle)
        hu.SetLineColor(self.linecolor)
        hu.SetLineWidth(self.linewidth)
        hu.SetFillStyle(self.fillstyle)
        hu.SetFillColor(self.fillcolor)
        hu.SetMarkerStyle(self.markerstyle)
        return hu


    

class merged_plotter(plotter_base):

    def __init__(self, plotters):
        super(merged_plotter,self).__init__()
        self.plotters = plotters

    def define(self, var, definition):
        for plotter in self.plotters:
            plotter.define(var, definition)

    def filter(self, condition):
        for plotter in self.plotters:
            plotter.filter(condition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter.redefine(var, definition)


    def addCorrectionFactor(self,value,model):
        for plotter in self.plotters:
            plotter.addCorrectionFactor(value,model)
            
    def event_weights(self,cuts):
        arr =None
        for plotter in self.plotters:
            if arr is None:
                arr = plotter.event_weights(cuts)
            else:
                extra=plotter.event_weights(cuts)
                arr=np.concatenate([arr,extra])
        return arr        
        

    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist1d(var, cuts, model, titlex, units)
            else:
                h.Add(plotter.hist1d(var, cuts, model, titlex, units))
        if h is None:
            tmprdf = ROOT.RDataFrame(1)
            tmprdf = tmprdf.Define("x", "0")
            h = tmprdf.Define("goodX", "x!=0").Histo1D(model, "goodX")
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        return h


    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist2d(var1,var2, cuts, model, titlex, unitsx, titley, unitsy)
            else:
                h.Add(plotter.hist2d(var1,var2, cuts, model, titlex, unitsx, titley, unitsy))
        if h is None:
            return h
        else:
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
            if unitsx=="":
                h.GetXaxis().SetTitle(titlex)
            else:
                h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
            if unitsy=="":
                h.GetYaxis().SetTitle(titley)
            else:
                h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
            return h

    def unrolled2d(self,var1,var2,cuts,model):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.unrolled2d(var1,var2, cuts, model)
            else:
                h.Add(plotter.unrolled2d(var1,var2, cuts, model))
        if h is None:
            return h
        else:
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
        return h
                               
    def getArray(self, var, cuts):
        """
        Applies cuts to RDF and extracts numpy array of values for a specific variable.

        Parameters:
            vars (str): variable
            cuts (str): Selection cuts

            range (tuple) : range over

        Returns:
            list of values corresponding to cuts
        """        
        # Initialize empty lists to store the data
        combined_var = []

        # Iterate over the plotters and collect data
        for plotter in self.plotters:
            # Define weight expression
            corrString = "1.0"
            for corr in plotter.corrFactors:
                corrString += f"*({corr['value']})"

            # Apply cuts and filters based on the provided cuts and lumi
            weight_expr = f"({plotter.defaultCuts})*({cuts})*{plotter.weight}*({corrString})"
            rdf = plotter.rdf.Define("plot_weight", weight_expr)
            rdf = rdf.Filter("plot_weight")

            # Get data from the current plotter
            numpy_var = rdf.AsNumpy([var])[var]

            # Collect the data
            combined_var.extend(numpy_var)

        # If no data points were collected, return None
        if len(combined_var) == 0:
            return None

        

        return combined_var
    
    def display(self,var):
        for plotter in self.plotters:
            plotter.display(var)

    def readReport(self):
        out = {}
        for plotter in self.plotters:
            r = plotter.readReport()
            for cut in r:
                if cut in out:
                    out[cut] += r[cut]
                else:
                    out[cut] = r[cut]
        return out



#Tricky, use with care!
class abcd_plotter(merged_plotter):
    def __init__(self,cutsSR,cutsCR,cutsSSB,cutsCSB,plotters):
        self.cutsSR=cutsSR
        self.cutsCR=cutsCR
        self.cutsSSB=cutsSSB
        self.cutsCSB=cutsCSB
        super(abcd_plotter,self).__init__(plotters)
        super(abcd_plotter,self).define('scaleVar','1.0')

    def getScale(self):       
        numeratorCuts = self.cutsSSB
        denominatorCuts = self.cutsCSB
#        print(numeratorCuts)
#        print(denominatorCuts)
        
        hNum=super(abcd_plotter,self).hist1d('scaleVar',numeratorCuts,('a','a',10,-5,5))
        #replace cuts with cuts from control region
        hDenom=super(abcd_plotter,self).hist1d('scaleVar',denominatorCuts,('a','a',10,-5,5))
        if hDenom.Integral()==0:
            print("Error , denominator has zero events, returning scale=1")
            return 1.0
        else:
            return hNum.Integral()/hDenom.Integral()        
        
    def hist1d(self,var,cuts,model,titlex = "",units = ""):
        h=super(abcd_plotter,self).hist1d(var,cuts.replace(self.cutsSR,self.cutsCR),model,titlex,units)
        h.Scale(self.getScale())
        return h
    
    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h=super(abcd_plotter,self).hist2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),model,titlex,unitsx,titley,unitsy)
        h.Scale(self.getScale())
        return h

    def unrolled2d(self,var1,var2,cuts,model):
        h=super(abcd_plotter,self).unrolled2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),model)      
        h.Scale(self.getScale())
        return h



#Tricky, use with care!
class fakerate_plotter(merged_plotter):
    def __init__(self,cutsSR,cutsCR,plotters,fakeRateVar,definition):
        self.cutsSR=cutsSR
        self.cutsCR=cutsCR
        super(fakerate_plotter,self).__init__(plotters)
        self.fakeRateVar=fakeRateVar
        self.define(fakeRateVar,definition)

    def hist1d(self,var,cuts,model,titlex = "",units = "",variation=0):
        h=super(fakerate_plotter,self).hist1d(var,'('+cuts.replace(self.cutsSR,self.cutsCR)+f")*({self.fakeRateVar}[{variation}])",model,titlex,units)
        return h
    
    def hist2d(self,var1,var2,cuts,model,titlex = "",unitsx = "",titley="",unitsy="",variation=0):
        h=super(fakerate_plotter,self).hist2d(var1,var2,'('+cuts.replace(self.cutsSR,self.cutsCR)+f")*({self.fakeRateVar}[{variation}])",model,titlex,unitsx,titley,unitsy)
        return h

    def unrolled2d(self,var1,var2,cuts,model,variation=0):
        h=super(fakerate_plotter,self).unrolled2d(var1,var2,'('+cuts.replace(self.cutsSR,self.cutsCR)+f")*({self.fakeRateVar}[{variation}])",model)      
        return h

    def event_weights(self,cuts,variation=0):
        arr =None
        for plotter in self.plotters:
            if arr is None:
                arr = plotter.event_weights('('+cuts.replace(self.cutsSR,self.cutsCR)+f")*({self.fakeRateVar}[{variation}])")
            else:
                extra=plotter.event_weights('('+cuts.replace(self.cutsSR,self.cutsCR)+f")*({self.fakeRateVar}[{variation}])")
                arr=np.concatenate([arr,extra])
        return arr        

    
class stack_plotter(object):
    def __init__(self,mode='stack',defaultCut="1"):
        self.plotters = []        
        self.log=False
        self.defaultCut=defaultCut
        self.mode=mode
    def setLog(self,doLog):
        self.log=doLog

    def add_plotter(self,plotter,name="",label = "label",typeP = "background"):
        packet = {'plotter':plotter,
                  'type':typeP,
                  'label':label,
                  'name':name}
        self.plotters.append(packet)
        
    def define(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].define(var, definition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].redefine(var, definition)

            
    def __getattr__(self, name):

        def fill_histograms(*args,**kwargs):    
            bkg_histograms=[]
            sig_histograms=[]
            data_histograms=[]
            
        
            #extract the histograms and build the stack
            for p in self.plotters:
                if hasattr(p['plotter'], name):
                    method = getattr(p['plotter'], name)
                    if p['type']=='background':
                        bkg_histograms.append(method(*args,**kwargs))
                        bkg_histograms[-1].SetName(p['name'])
                    elif p['type']=='signal':
                        sig_histograms.append(method(*args,**kwargs))
                        sig_histograms[-1].SetName(p['name'])
                                                   
                    elif p['type']=='data':
                        data_histograms.append(method(*args,**kwargs))
                        data_histograms[-1].SetName(p['name'])

            maxYields = [s.GetMaximum() for s in sig_histograms+bkg_histograms+data_histograms]
            integrals = [s.GetMaximum()/s.Integral() for s in sig_histograms+bkg_histograms+data_histograms]
            maxYield = max(maxYields)*1.2
            maxNormYield = max(integrals)*1.5
            #ok now stack all background hists if stack is the mode
            if self.mode=='stack':
                for i in range(0,len(bkg_histograms)-1):
                    for j in range(i,len(bkg_histograms)):
                        bkg_histograms[i].Add(bkg_histograms[j])
            #also add all signals separately to the background stack
            if len(bkg_histograms)>0:
                for h in sig_histograms:
                    h.Add(bkg_histograms[0])
                
                
        
            canvas = ROOT.TCanvas("canvas", "", 800, 600)
            canvas.SetTopMargin(0.055) # testing tdr style
            canvas.SetBottomMargin(0.12) # testing
            canvas.SetLeftMargin(0.12) # testing

            canvas.cd()
            #first the background
            for i in range(0,len(bkg_histograms)):
                if i==0:
                    if self.mode=='stack':
                        bkg_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        bkg_histograms[i].Draw("HIST")
                        
                    else:
                        bkg_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)                        
                        bkg_histograms[i].DrawNormalized("HIST")                    
                else:    
                    if self.mode=='stack':
                        bkg_histograms[i].Draw("HIST,SAME")
                    else:
                        bkg_histograms[i].DrawNormalized("HIST,SAME")                    
            #then the signals
            for i in range(0,len(sig_histograms)):
                if i==0 and len(bkg_histograms)==0:
                    if self.mode=='stack':
                        sig_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        sig_histograms[i].Draw("HIST")
                    else:
                        sig_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)
                        sig_histograms[i].DrawNormalized("HIST")                        
                else:
                    if self.mode=='stack':
                        sig_histograms[i].Draw("HIST,SAME")
                    else:
                        sig_histograms[i].DrawNormalized("HIST,SAME")                    
            
            #then the data
            for i in range(0,len(data_histograms)):
                if i==0 and len(bkg_histograms)==0 and len(sig_histograms)==0:
                    if self.mode=='stack':
                        data_histograms[i].GetYaxis().SetRangeUser(0,maxYield)
                        data_histograms[i].Draw("")
                    else:
                        data_histograms[i].GetYaxis().SetRangeUser(0,maxNormYield)
                        data_histograms[i].DrawNormalized("")
                else:        
                    if self.mode=='stack':
                        data_histograms[i].Draw("SAME")
                    else:
                        data_histograms[i].DrawNormalized("SAME")

            #Then the legend            
            legend = ROOT.TLegend(0.61,0.64,0.91,0.94,"","brNDC") # if using tdrstyle
            legend.SetBorderSize(0)
            legend.SetLineColor(1)
            legend.SetLineStyle(1)
            legend.SetLineWidth(1)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)

            legend.SetFillColor(ROOT.kWhite)
            i=0
            for p in self.plotters:
                if p['type']=='data':
                    legend.AddEntry(data_histograms[i], p['label'], 'p')
                    i=i+1
            i=0
            for p in self.plotters:
                if p['type']=='signal':
                    legend.AddEntry(sig_histograms[i], p['label'], 'l')
                    i=i+1

            i=0
            for p in self.plotters:
                if p['type']=='background':
                    legend.AddEntry(bkg_histograms[i], p['label'], 'lf')
                    i=i+1
                
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            ROOT.SetOwnership(legend, False)
            legend.Draw()
            canvas.Update()
            
            return {'canvas': canvas, 'legend': legend,'hists': sig_histograms+bkg_histograms+data_histograms}
    
        return fill_histograms



#MPLHEP specific plotter
class mplhep_plotter(object):
    def __init__(self,label=None,lumi=137.62,data=True,com=13,defaultCut="1",stack=True,capsize=5):
        self.plotters = []        
        self.defaultCut=defaultCut
        mh.style.use('CMS')
        self.label=label
        self.lumi=lumi
        self.com=com
        self.data=data
        self.stack=stack
        self.capsize=capsize
        
    def add_plotter(self,plotter,name='name',label = "label",typeP = "background",error_mode='w2',color='black'):
        packet = {'plotter':plotter,
                  'type':typeP,
                  'error_mode':error_mode,
                  'label':label,
                  'name':name,
                  'color':color}
        self.plotters.append(packet)
        
    def define(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].define(var, definition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter['plotter'].redefine(var, definition)

    def hist1d(self,var,cuts,model,alpha=1.0,xlabel="",xunits="",legend_loc='upper right',show=True,logscale=False):
        background_hists=[]
        background_edges=[]
        background_w2=[]
        background_labels=[]
        background_colors=[]
        data_hists=[]
        data_edges=[]
        data_w2=[]
        data_labels=[]
        data_colors=[]
        signal_hists=[]
        signal_edges=[]
        signal_w2=[]
        signal_labels=[]
        signal_colors=[]

        bkgExists=False        
        for p in self.plotters:
            if p['type']=='data':
                edges,data,w2=p['plotter'].array1d(var,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)
                data_hists.append(data)
                data_edges.append(edges)
                data_w2.append(w2)
                data_labels.append(p['label'])
                data_colors.append(p['color'])                                    
            elif p['type']=='background':
                edges,data,w2=p['plotter'].array1d(var,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)
                
                if bkgExists==False:
                    background_sum=data
                    background_sumw2=w2
                    bkgExists=True
                else:
                    background_sum=background_sum+data
                    background_sumw2=background_sumw2+w2

                background_hists.append(data)
                background_edges.append(edges)
                background_w2.append(w2)
                background_labels.append(p['label'])
                background_colors.append(p['color'])                                                    
        for p in self.plotters:                   
            if p['type']=='signal':
                edges,data,w2=p['plotter'].array1d(var,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)
                
                if self.stack==True:
                    signal_hists.append(data+background_sum)
                    signal_w2.append(w2+background_sumw2)
                else:
                    signal_hists.append(data)
                    signal_w2.append(w2)                    
                signal_edges.append(edges)
                signal_labels.append(p['label'])
                signal_colors.append(p['color'])                                    

        fig,ax = plt.subplots()
        if len(signal_hists)>0:
            mh.histplot(signal_hists,signal_edges[0],
                        histtype='step',
                        stack=False,
                        label=signal_labels,
                        color=signal_colors,
                        yerr=None,
                        ax=ax,
                        density=(True if self.stack==False else False)                        
                        )                   
        if len(background_hists)>0:
            #plot background stack            
            mh.histplot(background_hists,background_edges[0],
                        histtype='fill',
                        stack=self.stack,
                        label=background_labels,
                        ax=ax,
                        alpha=(alpha if self.stack==False else 1.0), 
                        density=(True if self.stack==False else False)
                        )
            #plot background error band in a custom way (since we use old version of mplhep)
            if self.stack:
                mh.histplot(background_hists,background_edges[0],
                            histtype='step',
                            color=(['black']*len(background_hists)),
                            stack=self.stack,
                            ax=ax
                            )
                
                ax.fill_between(background_edges[0],
                                np.append(background_sum-np.sqrt(background_sumw2[0]),0),
                                np.append(background_sum+np.sqrt(background_sumw2[1]),0),
                                step='post', color='lightgray', alpha=0.5, hatch='////')            
        if len(data_hists)>0:           
            mh.histplot(data_hists,data_edges[0],
                        histtype='errorbar',
                        stack=False,
                        label=data_labels,
                        color=data_colors,
                        yerr = [np.sqrt(a) for a in data_w2],
                        capsize=self.capsize,
                        ax=ax,
                        density=(True if self.stack==False else False)                        
                        )

        #then stack backgrounds and then draw band
        ax.legend(loc=legend_loc)
        if self.data==True:
            mh.cms.label(self.label, data=self.data, lumi=self.lumi, com=self.com,ax=ax, loc=0)
        else:
            mh.cms.label(self.label, data=self.data, lumi=None, com=None,rlabel=f'{self.com} TeV',ax=ax, loc=0)
            
        #fix the lower limit
        lims=plt.ylim()
        plt.ylim(0.0, lims[1])
        if xlabel!="":
            if xunits!='':
                ax.set_xlabel(f"{xlabel} ({xunits})")
            else:
                ax.set_xlabel(f"{xlabel}")                
        if self.stack:
            ax.set_ylabel("Events")
        else:
           ax.set_ylabel("a.u")

        if logscale:            
            if self.stack==False:
                ax.set_ylim(0.0001,100)
            ax.set_yscale('log')    
        if show:
            plt.show()
            

    def unrolled2d(self,var1,var2,cuts,model,alpha=1.0,xlabel="",xunits="",legend_loc='upper right',show=True):
        background_hists=[]
        background_edges=[]
        background_w2=[]
        background_labels=[]
        background_colors=[]
        data_hists=[]
        data_edges=[]
        data_w2=[]
        data_labels=[]
        data_colors=[]
        signal_hists=[]
        signal_edges=[]
        signal_w2=[]
        signal_labels=[]
        signal_colors=[]

        bkgExists=False        
        for p in self.plotters:
            if p['type']=='data':
                xedges,yedges,data,w2=p['plotter'].array2d(var1,var2,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)

                data_hists.append(data)
                data_edges.append((xedges,yedges))
                data_w2.append(w2)
                data_labels.append(p['label'])
                data_colors.append(p['color'])                                    
            elif p['type']=='background':
                xedges,yedges,data,w2=p['plotter'].array2d(var1,var2,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)
                
                if bkgExists==False:
                    background_sum=data
                    background_sumw2=w2
                    bkgExists=True
                else:
                    background_sum=background_sum+data
                    background_sumw2=background_sumw2+w2

                background_hists.append(data)
                background_edges.append((xedges,yedges))
                background_w2.append(w2)
                background_labels.append(p['label'])
                background_colors.append(p['color'])                                                    
        for p in self.plotters:                   
            if p['type']=='signal':
                xedges,yedges,data,w2=p['plotter'].array2d(var1,var2,cuts,model,error_mode=p['error_mode'])
                if self.stack==False:
                    w2=w2/np.sum(data)
                    data=data/np.sum(data)
                
                if self.stack==True:
                    signal_hists.append(data+background_sum)
                    signal_w2.append(w2+background_sumw2)
                else:
                    signal_hists.append(data)
                    signal_w2.append(w2)                    
                signal_edges.append((xedges,yedges))
                signal_labels.append(p['label'])
                signal_colors.append(p['color'])                                    

        
        if len(background_hists)>0:
            xedges=background_edges[0][0]
            yedges=background_edges[0][1]
        elif len(signal_hists)>0:
            xedges=signal_edges[0][0]
            yedges=signal_edges[0][1]
        elif len(data_hists)>0:
            xedges=data_edges[0][0]
            yedges=data_edges[0][1]
        fig,ax = plt.subplots(1,len(yedges)-1,sharey=True,figsize=(25,10))
        plt.subplots_adjust(wspace=0)        
        for i in range(0,len(yedges)-1):
            if len(signal_hists)>0:
                mh.histplot([arr[i, :] for arr in signal_hists],xedges,
                            histtype=('step' if self.stack==True else 'fill'),
                            stack=False,
                            label=signal_labels,
                            color=signal_colors,
                            yerr=None,
                            ax=ax[i],
                            density=(True if self.stack==False else False)                        
                            )                   
            if len(background_hists)>0:
                #plot background stack            
                mh.histplot([arr[i, :] for arr in background_hists],xedges,
                            histtype=('fill' if self.stack==True else 'step'),
                            stack=self.stack,
                            label=background_labels,
                            ax=ax[i],
                            density=(True if self.stack==False else False)
                            )
                #plot background error band in a custom way (since we use old version of mplhep)
                if self.stack:
                    mh.histplot([arr[i, :] for arr in background_hists],xedges,
                            histtype='step',
                            color=(['black']*len(background_hists)),
                            stack=self.stack,
                            ax=ax[i]
                            )                    
                    ax[i].fill_between(xedges,
                                    np.append(background_sum[i,:]-np.sqrt(background_sumw2[0])[i,:],0),
                                    np.append(background_sum[i,:]+np.sqrt(background_sumw2[1])[i,:],0),
                                    step='post', color='lightgray', alpha=0.5, hatch='////')            
                if len(data_hists)>0:           
                    mh.histplot([arr[i, :] for arr in data_hists],xedges,
                                histtype='errorbar',
                                stack=False,
                                label=data_labels,
                                color=data_colors,
                                yerr = [np.array([np.sqrt(a[0,i,:]),np.sqrt(a[1,i,:])]) for a in data_w2],
                                capsize=self.capsize,                                
                                ax=ax[i],
                                density=(True if self.stack==False else False)                        
                                )

        #then stack backgrounds and then draw band
        ax[-1].legend(loc=legend_loc)

        #Labels
        mh.cms.label(self.label,data=self.data,rlabel="", ax=ax[0], loc=0)
        if self.data==False:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel=f'{self.com} TeV')
        else:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=self.lumi,com=self.com)

        
        #fix the lower limit
        if xlabel!="":
            ax[-1].set_xlabel(f"{xlabel} ({xunits})")
        if self.stack:
            ax[0].set_ylabel("Events")
        else:
           ax[0].set_ylabel("a.u")
             
        if show:
            plt.show()
        
    def unrolledCustom(self,var1,var2,cuts,custom_binning,alpha=1.0,xlabel="",xunits="",ylabel="",yunits="",legend_loc='upper right',legend_ax=-1,show=True,textx=0.5,texty=0.5):

        fig,ax = plt.subplots(1,len(custom_binning),sharey=True,figsize=(25,10))
        plt.subplots_adjust(wspace=0)        
        yields={}
        
        for i,binning in enumerate(custom_binning):
            xedges=binning[0]
            yedges=binning[1]
            
            background_hists=[]
            background_w2=[]
            background_labels=[]
            background_colors=[]
            data_hists=[]
            data_w2=[]
            data_labels=[]
            data_colors=[]
            signal_hists=[]
            signal_w2=[]
            signal_labels=[]
            signal_colors=[]

            bkgExists=False        
            for p in self.plotters:
            
                if p['type']=='data':
                    edg,data,w2=p['plotter'].array1d(var2,cuts+f"*({var1}>={xedges[0]}&&{var1}<{xedges[1]})",(p['name'],p['name'],len(yedges)-1,np.array(yedges)),error_mode=p['error_mode'])
                    if self.stack==False:
                        w2=w2/np.sum(data)
                        data=data/np.sum(data)
                    
                    if p['label'] in yields.keys():
                        yields[p['label']]=yields[p['label']]+np.sum(data)
                    else:
                        yields[p['label']]=np.sum(data)
                        
                    data_hists.append(data)
                    data_w2.append(w2)
                    data_labels.append(p['label'])
                    data_colors.append(p['color'])                                    
                elif p['type']=='background':
                    edg,data,w2=p['plotter'].array1d(var2,cuts+f"*({var1}>={xedges[0]}&&{var1}<{xedges[1]})",(p['name'],p['name'],len(yedges)-1,np.array(yedges)),error_mode=p['error_mode'])
                    if self.stack==False:
                        w2=w2/np.sum(data)
                        data=data/np.sum(data)
                    
                    if p['label'] in yields.keys():
                        yields[p['label']]=yields[p['label']]+np.sum(data)
                    else:
                        yields[p['label']]=np.sum(data)
                    
                    if bkgExists==False:
                        background_sum=data
                        background_sumw2=w2
                        bkgExists=True
                    else:
                        background_sum=background_sum+data
                        background_sumw2=background_sumw2+w2
                    background_hists.append(data)
                    background_w2.append(w2)
                    background_labels.append(p['label'])
                    background_colors.append(p['color'])                                                    
            for p in self.plotters:                
                if p['type']=='signal':
                    edg,data,w2=p['plotter'].array1d(var2,cuts+f"*({var1}>={xedges[0]}&&{var1}<{xedges[1]})",(p['name'],p['name'],len(yedges)-1,np.array(yedges)),error_mode=p['error_mode'])
                    if self.stack==False:
                        w2=w2/np.sum(data)
                        data=data/np.sum(data)
                    
                    if p['label'] in yields.keys():
                        yields[p['label']]=yields[p['label']]+np.sum(data)
                    else:
                        yields[p['label']]=np.sum(data)
                    
                    if self.stack==True:
                        if bkgExists:
                            signal_hists.append(data+background_sum)
                            signal_w2.append(w2+background_sumw2)                        
                        else:
                            signal_hists.append(data)                            
                            signal_w2.append(w2)
                    else:
                        signal_hists.append(data)
                        signal_w2.append(w2)                    
                    signal_labels.append(p['label'])
                    signal_colors.append(p['color'])                                    

        

            if len(background_hists)>0:
                #plot background stack
                if self.stack:
                    mh.histplot(background_hists,yedges,
                                histtype='fill',
                                stack=True,
                                label=background_labels,
                                ax=ax[i]
                                )
                    ax[i].fill_between(yedges,
                                       np.append(background_sum-np.sqrt(background_sumw2[0]),0),
                                       np.append(background_sum+np.sqrt(background_sumw2[1]),0),
                                       step='post', color='lightgray', alpha=0.5, hatch='////')            
                    
                else:
                    mh.histplot(background_hists,yedges,
                                histtype='fill',
                                stack=False,
                                label=background_labels,
                                alpha=alpha,
                                ax=ax[i],
                                density=True
                                )
            #then signal stack        
            if len(signal_hists)>0:
                mh.histplot(signal_hists,yedges,
                            histtype='step',
                            stack=False,
                            label=signal_labels,
                            color=signal_colors,
                            yerr=None,
                            ax=ax[i],
                            density=(True if self.stack==False else False)                        
                            )                   

            #Now redraw the background outlines if we are doing stack
            if len(background_hists)>0 and self.stack:
                mh.histplot(background_hists,yedges,
                            histtype='step',
                            color=['black']*len(background_hists),
                            stack=True,
                            label=None,
                            ax=ax[i]
                            )
                    
                    
            
            if len(data_hists)>0:           
                    mh.histplot(data_hists,yedges,
                                histtype='errorbar',
                                stack=False,
                                label=data_labels,
                                color=data_colors,
                                yerr = [np.sqrt(a) for a in data_w2],
                                capsize=self.capsize,                                
                                ax=ax[i],
                                density=(True if self.stack==False else False)                        
                                )
            ax[i].text(textx, texty, f'{xedges[0]}'+r"$\leq$"+f"{ylabel} <{xedges[1]} {yunits}" , transform=ax[i].transAxes,
                       fontsize=20, ha='center', va='center', 
                       bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            #then stack backgrounds and then draw band
        ax[legend_ax].legend(loc=legend_loc)

        #Labels
        mh.cms.label(self.label,data=self.data,rlabel="", ax=ax[0], loc=0)
        if self.data==False:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel=f'{self.com} TeV')
        else:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=self.lumi,com=self.com)
        #fix the lower limit
        if xlabel!="":
            ax[-1].set_xlabel(f"{xlabel} ({xunits})")
        if self.stack:
            ax[0].set_ylabel("Events")
        else:
           ax[0].set_ylabel("a.u")


        #print yields
        for label,y in yields.items():
            print(label,y)


           
        if show:
            plt.show()
        


#class for limits
class limit_plotter(object):
    def __init__(self,filename,label=None,lumi=137.62,data=True,com=13):
        mh.style.use('CMS')
        self.label=label
        self.lumi=lumi
        self.com=com
        self.data=data
        self.rdf = ROOT.RDataFrame('limit',filename)

    def brazilian_flag(self,band2sigma_color='lime',band1sigma_color='yellow',ax=None,show=False,quiet=False,legend_loc='upper right',xlabel='',xunits='',ylabel='95% CL Upper Limits'):
        expected025=self.rdf.Filter('TMath::Abs(quantileExpected-0.025)<0.001').AsNumpy(['mh','limit'])
        expected16=self.rdf.Filter('TMath::Abs(quantileExpected-0.16)<0.001').AsNumpy(['mh','limit'])
        expected=self.rdf.Filter('TMath::Abs(quantileExpected-0.5)<0.001').AsNumpy(['mh','limit'])
        expected84=self.rdf.Filter('TMath::Abs(quantileExpected-0.84)<0.001').AsNumpy(['mh','limit'])
        expected975=self.rdf.Filter('TMath::Abs(quantileExpected-0.975)<0.001').AsNumpy(['mh','limit'])
        observed=self.rdf.Filter('quantileExpected==-1').AsNumpy(['mh','limit'])
#        import pdb;pdb.set_trace()
        if ax is None:
            fig,ax = plt.subplots()
        
            
        ax.fill_between(expected025['mh'],
                        expected025['limit'],
                        expected975['limit'],
                        color=band2sigma_color, alpha=1,label=r'$\pm$ 2 std. deviations')            
        ax.fill_between(expected16['mh'],
                        expected16['limit'],
                        expected84['limit'],
                        color=band1sigma_color, alpha=1,label=r'$\pm$ 1 std. deviation')            
        
        ax.plot(expected['mh'],expected['limit'],linestyle='--',color='k',label=r'Asymptotic $CL_s$ expected')
        ax.plot(observed['mh'],observed['limit'],color='k',label='Observed')
        ax.margins(x=0)
                     
        if quiet==False:
            ax.legend(loc=legend_loc)
            if xlabel!="":
                ax.set_xlabel(f"{xlabel} ({xunits})")
            if ylabel!="":
                ax.set_ylabel(ylabel)
            if self.data==True:
                mh.cms.label(self.label, data=self.data, lumi=self.lumi, com=self.com,ax=ax, loc=0)
            else:
                mh.cms.label(self.label, data=self.data, lumi=None, com=None,rlabel=f'{self.com} TeV',ax=ax, loc=0)
        if show:
            plt.show()
        
    def expected_vs_observed(self,col='lime',ax=None,show=False,quiet=False,legend_loc='upper right',xlabel='',xunits='',ylabel='95% CL Upper Limits'):
        expected=self.rdf.Filter('TMath::Abs(quantileExpected-0.5)<0.001').AsNumpy(['mh','limit'])
        observed=self.rdf.Filter('quantileExpected==-1').AsNumpy(['mh','limit'])

        if ax is None:
            fig,ax = plt.subplots()
        
            
        ax.plot(expected['mh'],expected['limit'],linestyle='--',color=col,label=r'Asymptotic $CL_s$ expected')
        ax.plot(observed['mh'],observed['limit'],color=col,label='Observed')
        ax.margins(x=0)
                     
        if quiet==False:
            ax.legend(loc=legend_loc)
            if xlabel!="":
                ax.set_xlabel(f"{xlabel} ({xunits})")
            if ylabel!="":
                ax.set_ylabel(ylabel)
            if self.data==True:
                mh.cms.label(self.label, data=self.data, lumi=self.lumi, com=self.com,ax=ax, loc=0)
            else:
                mh.cms.label(self.label, data=self.data, lumi=None, com=None,rlabel=f'{self.com} TeV',ax=ax, loc=0)
        if show:
            plt.show()


class multilimit_plotter(object):
    def __init__(self,plotters,label=None,lumi=137.62,data=True,com=13):
        mh.style.use('CMS')
        self.label=label
        self.lumi=lumi
        self.com=com
        self.data=data

        self.plotters=plotters

    def brazilian_flag(self,band2sigma_color='lime',band1sigma_color='yellow',show=False,legend_loc='upper right',xlabel='',xunits='',ylabel='95% CL Upper Limits',log=True):
        
        fig,ax = plt.subplots(1,len(self.plotters),sharey=True,figsize=(7.5*len(self.plotters),20))
        plt.subplots_adjust(wspace=0)
        for i , plotter in enumerate(self.plotters):
            xdim = int(i%3)
            ydim = int(i/3)

            plotter['plotter'].brazilian_flag(band2sigma_color=band2sigma_color,band1sigma_color=band1sigma_color,ax=ax[i],show=False,quiet=True)
            
            
        ax[-1].legend(loc=legend_loc)
        if xlabel!="":
            ax[-1].set_xlabel(f"{xlabel} ({xunits})")
        if ylabel!="":
            ax[0].set_ylabel(ylabel)
        mh.cms.label(self.label,data=self.data,rlabel="", ax=ax[0], loc=0)
        if self.data==False:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel=f'{self.com} TeV')
        else:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=self.lumi,com=self.com)
        if log:
            ax[0].set_yscale('log')
        if show:
            plt.show()
        
    def expected_vs_observed(self,show=False,legend_loc='upper right',xlabel='',xunits='',ylabel='95% CL Upper Limits',log=True):
        
        fig,ax = plt.subplots(1,len(self.plotters),sharey=True,figsize=(7.5*len(self.plotters),20))
        plt.subplots_adjust(wspace=0)
        for i , plotter in enumerate(self.plotters):
            xdim = int(i%3)
            ydim = int(i/3)

            plotter['plotter'].expected_vs_observed(col=plotter['color'],ax=ax[i],show=False,quiet=True)
            
            
        ax[-1].legend(loc=legend_loc)
        if xlabel!="":
            ax[-1].set_xlabel(f"{xlabel} ({xunits})")
        if ylabel!="":
            ax[0].set_ylabel(ylabel)
        mh.cms.label(self.label,data=self.data,rlabel="", ax=ax[0], loc=0)
        if self.data==False:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel=f'{self.com} TeV')
        else:
            mh.cms.label(None,exp='',data=self.data,llabel="", ax=ax[-1], loc=0,lumi=self.lumi,com=self.com)
        if log:
            ax[0].set_yscale('log')
        if show:
            plt.show()
        
        




            
            
        



        
