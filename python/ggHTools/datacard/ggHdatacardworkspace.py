import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json
from ggHparameters import signal_window, bernstein_coeff_card


class DatacardWorkspace:
    def __init__(self,finalstate,category,period,lifetime,mass,luminosity=1.0,physics="ggH"):
        self.physics=str(physics)
        self.lifetime=str(lifetime)
        self.mass=str(mass)
        self.finalstate=str(finalstate)
        self.category=str(category)
        self.period=str(period)
        self.contributions=[]
        self.systematics=[]

        self.tag=self.physics+"_"+self.finalstate+"_m"+self.mass+"_ct"+self.lifetime+"_"+self.category+"_"+self.period
        self.rootFile = ROOT.TFile("datacardInputs_"+self.tag+".root","RECREATE")
        self.rootFile.cd()
        self.w=ROOT.RooWorkspace("w","w")
        self.luminosity=luminosity
        self.w.factory(self.physics+"_"+period+"_lumi["+str(luminosity)+"]")


    def addSystematic(self,name,kind,values):
        self.systematics.append({'name':name,'kind':kind,'values':values })


    def addDCB(self, name, variable, jsonFile, scale={}, resolution={}):
        self.w.factory("MH[125.0]")
        self.w.var("MH").setConstant(1)
        scaleStr = '0'
        resolutionStr = '0'
        scaleSysts = []
        resolutionSysts = []

        for syst, factor in scale.items():
            self.w.factory(f"{syst}[0,-0.1,0.1]")
            scaleStr += f"+{factor}*{syst}"
            scaleSysts.append(syst)
        for syst, factor in resolution.items():
            self.w.factory(f"{syst}[0,-5,5]")
            resolutionStr += f"+{factor}*{syst}"
            resolutionSysts.append(syst)
        self.w.factory(f"{variable}[{signal_window[0]},{signal_window[1]}]")

        with open(jsonFile) as f:
            info = json.load(f)

            SCALEVar="_".join(["MEAN",name,self.tag])
            self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(
            name=SCALEVar,
            param=info['mean']['value'],
            vv_syst=scaleStr,
            vv_systs=','.join(scaleSysts)))

            SIGMAVar="_".join(["SIGMA",name,self.tag])
            self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(
            name=SIGMAVar,
            param=info['sigma']['value'],
            vv_syst=resolutionStr,
            vv_systs=','.join(resolutionSysts)))

            ALPHA1Var="_".join(["ALPHA1",name,self.tag])
            self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
            name=ALPHA1Var,
            param=info['alpha1']['value']))

            ALPHA2Var="_".join(["ALPHA2",name,self.tag])
            self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
            name=ALPHA2Var,
            param=info['alpha2']['value']))

            N1Var="_".join(["N1",name,self.tag])
            self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
            name=N1Var,
            param=info['n1']['value']))

            N2Var="_".join(["N2",name,self.tag])
            self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
            name=N2Var,
            param=info['n2']['value']))

            pdfName="_".join([name,self.tag])


        fourPhotonMass = ROOT.RooDoubleCB(pdfName, pdfName,
                                          self.w.var(variable),
                                          self.w.function(SCALEVar),
                                          self.w.function(SIGMAVar),
                                          self.w.function(ALPHA1Var),
                                          self.w.function(N1Var),
                                          self.w.function(ALPHA2Var),
                                          self.w.function(N2Var))
        self.w.Import(fourPhotonMass)


    def addBernstein(self, name, variable, jsonFile):
        self.w.factory(f"{variable}[{signal_window[0]},{signal_window[1]}]")
        with open(jsonFile) as f:
            params = json.load(f)
        order = 0
        while f"c_{order}" in params:
            order += 1
        if order == 0:
            raise ValueError(f"No Bernstein coefficients (c_0, c_1, ...) found in {jsonFile}")
        coeffNames = [f"c_{i}_bkg" for i in range(order)]
        keys = [f"c_{i}" for i in range(order)]
        clist=ROOT.RooArgList()
        for coeffName, key in zip(coeffNames, keys):
            self.w.factory("{name}[{val}, {lo}, {hi}]".format(name=coeffName, val=params[key]['value'], lo=bernstein_coeff_card[0], hi=bernstein_coeff_card[1]))
            clist.add(self.w.var(coeffName))
        pdfName = "_".join([name, self.tag])
        bernstein_pdf = ROOT.RooBernsteinFast(order)(pdfName,pdfName,self.w.var(variable),clist)
        self.w.Import(bernstein_pdf)


    def addFixedYield(self,name,ID,value):
        pdfName="_".join([name,self.tag])
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':value})


    def addFixedYieldFromFile(self,name,ID,filename,histoName,lumi=True):
        pdfName="_".join([name,self.tag])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        if not lumi:
            events=histogram.Integral()
            self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})
        else:
            events=histogram.Integral()*self.luminosity
            self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})
        f.Close()


    def makeCard(self):

        f = open("datacard_"+self.tag+'.txt','w')
        f.write('imax 1\n')
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        for c in self.contributions:
            f.write('shapes {name} {channel} {file}.root w:{pdf}\n'.format(name=c['name'],channel=self.tag,file="datacardInputs_"+self.tag,pdf=c['pdf']))
        f.write('shapes {name} {channel} {file}.root w:{name}\n'.format(name="data_obs",channel=self.tag,file="datacardInputs_"+self.tag))
        f.write('-------------------------\n')
        f.write('bin '+self.tag+'\n')
        f.write('observation  -1\n')
        f.write('-------------------------\n')
        f.write('bin\t')

        for shape in self.contributions:
            f.write(self.tag+'\t')
        f.write('\n')

        #Sort the shapes by ID

        shapes = sorted(self.contributions,key=lambda x: x['ID'])
        #print names
        f.write('process\t')
        for shape in shapes:
            f.write(shape['name']+'\t')
        f.write('\n')

        #Print ID
        f.write('process\t')
        for shape in shapes:
            f.write(str(shape['ID'])+'\t')
        f.write('\n')

        #print rates
        f.write('rate\t')
        for shape in shapes:
            if shape['yield']==0:
                f.write("0.00000000001" + '\t')
            else:
                f.write(str(shape['yield'])+'\t')
        f.write('\n')


        #Now systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')

            elif syst['kind'] == 'discrete':
                f.write(syst['name']+'\t'+'discrete\n')

            elif syst['kind'] == 'lnN':
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].items():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
            elif syst['kind'] == 'lnU':
                f.write(syst['name']+'\t'+ 'lnU\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].items():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )

            elif syst['kind'] == 'rateParam':
                f.write(syst['name']+'\t'+'rateParam\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\t'+str(syst['values'][2])+'\t'+str(syst['values'][3])+'\n')



        f.close()


        self.rootFile.cd()
        self.w.Write()
        self.rootFile.Close()


    def importBinnedData(self,filename,histoname,poi,name = "data_obs",scale=1):
        f=ROOT.TFile(filename)
        histogram=f.Get(histoname)
        histogram.Scale(scale)
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            self.w.factory("{}[{},{}]".format(p, signal_window[0], signal_window[1]))
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print ('Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead')
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)

        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))
        f.Close()
