#Datacard maker for cut and count
import numpy as np
class cnc_datacard_maker(object):
    def __init__(self,outDir,binname,cuts):
        self.rates={}
        self.nuisances={}
        self.types={}
        self.data=-1        
        self.binname=binname
        self.cuts=cuts
        self.outDir=outDir


    def add(self,name,typeP,plotter,uncertainties,error_mode='w2'):
        plotter.define("dummyDCVar",'1.0')
        edges,data,w2=plotter.array1d('dummyDCVar',self.cuts,(name,name,2,0,2),error_mode=error_mode)
        rate=float(np.sum(data))
        if typeP=='signal' and rate==0.0:
            print (f"No yield for {name} in this card, skipping this contribution")
            return
        if typeP=='data':
            self.data=int(rate)
        elif typeP=='signal':
            self.types[name]='signal'
            self.rates[name]=rate
        elif typeP=='background':
            self.types[name]='background'
            self.rates[name]=rate
        error=np.sqrt(w2)
        for unc_name,unc in uncertainties.items():
            nData=None                
            block=False
            if unc['type'] =='statSym':
                if rate==0.0:
                    continue
                lower=np.sum(error[0,:])
                upper=np.sum(error[1,:])
                err=0.5*(lower+upper)
                err=float(np.divide(err,data))
                if err==0.0:
                    block=True
                nData = ('lnN',str(1+err))
                
            elif unc['type'] =='statAsym':
                if rate==0.0:
                    continue

                lower=np.sum(error[0,:])
                upper=np.sum(error[1,:])
                if upper==rate and lower==rate:
                    print("I detected the case of a single event where the Poisson bootstrap is applied\n This is fine but I am changing the uncertainty to gmN for this card\n")
                    nData = ('gmN 1',str(rate))
                else:
                    errUp=float(np.divide(upper,rate))
                    errDown=float(np.divide(lower,rate))
                    if errUp==0.0 or errDown==0.0:
                        block=True
                    nData = ('lnN',str(1-errDown)+'/'+str(1+errUp))
            elif unc['type'] =='weightSym':
                if rate==0.0:
                    continue                
                weight=unc['weight']
                weightOrig = unc['weightOrig']
                edges2,data2,w22=plotter.array1d('dummyDCVar',f"({self.cuts})*({weight}/({weightOrig}))",(name,name,2,0,2),error_mode='w2')
                err  = float(np.sum(data2)/rate)
                if err==0.0:
                    block=True
                nData = ('lnN',str(err))
            elif unc['type'] =='replication':
                if rate==0.0:
                    continue                
                cutsUp=self.cuts
                cutsDown=self.cuts
                for o,u,d in zip(unc['originals'],unc['replacementsUp'],unc['replacementsDown']):
                    cutsUp=cutsUp.replace(o,u)
                    cutsDown=cutsDown.replace(o,d)
#                print(self.cuts,cutsUp,cutsDown)
                edgesUp,up,wUp=plotter.array1d('dummyDCVar',f"({cutsUp})",(name,name,2,0,2),error_mode='w2')
                edgesDown,down,wDown=plotter.array1d('dummyDCVar',f"({cutsDown})",(name,name,2,0,2),error_mode='w2')
                errorUp = float(np.sum(up)/rate)
                errorDown = float(np.sum(down)/rate)
                if errorUp==0.0 or errorDown==0.0:
                    block=True
                nData = ('lnN',str(errorDown)+'/'+str(errorUp))
            elif unc['type'] =='weightAsymm':
                if rate==0.0:
                    continue
                
                weightUp=unc['weightUp']
                weightDown=unc['weightDown']
                #we also need to cancel out the original weight:
                weightOrig = unc['weightOrig']
                edges2,up,w2=plotter.array1d('dummyDCVar',f"({self.cuts})*({weightUp}/({weightOrig}))",(name,name,2,0,2),error_mode='w2')
                edges2,down,w2=plotter.array1d('dummyDCVar',f"({self.cuts})*({weightDown}/({weightOrig}))",(name,name,2,0,2),error_mode='w2')
#                print(f"({self.cuts})*({weightUp}/{weightOrig})")
                errorUp = float(np.sum(up)/rate)
                errorDown = float(np.sum(down)/rate)
                if errorUp==0.0 or errorDown==0.0:
                    block=True
                nData = ('lnN',str(errorDown)+'/'+str(errorUp))
            elif unc['type'] =='adhoc':
                nData = (unc['kind'],unc['value'])
            elif unc['type'] =='zeroRate' and rate==0.0:
                nData = ('gmN 0',unc['value'])    
            if block:
                continue
            if nData==None:
                continue
            if unc_name in self.nuisances.keys():
                self.nuisances[unc_name]['contribs'][name] = nData[1]
            else:
                self.nuisances[unc_name]={'kind':nData[0],
                                     'contribs':{name:nData[1]}
                                     }
    def write(self):

        with open(self.outDir+'/datacard_'+self.binname+'.txt','w') as file:
            file.write(f"# {self.cuts}]\n")
            file.write('imax 1 number of channels\n')
            nbackgrounds = sum([1 for x in self.types.keys() if self.types[x]=='background'])
            nsignals     = sum([1 for x in self.types.keys() if self.types[x]=='signal'])
            file.write(f"jmax *\n")
            file.write(f'kmax {len(self.nuisances.keys())} number of nuisance parameters\n')            
            file.write(f"bin {self.binname}\n")
            file.write(f"observation {self.data}\n")
            file.write('bin\t')
            for s in self.rates.keys():
                file.write(f"{self.binname}\t")
            file.write('\n')
            
            file.write("process\t")
            for s in self.rates.keys():
                file.write(f"{s}\t")
            file.write('\n')
            file.write("process\t")
            SIG_COUNTER=-nsignals+1
            BKG_COUNTER=1
            for name,t in self.types.items():
                if t=='background':
                    file.write(f"{BKG_COUNTER}\t")
                    BKG_COUNTER=BKG_COUNTER+1
                if t=='signal':
                    file.write(f"{SIG_COUNTER}\t")
                    SIG_COUNTER=SIG_COUNTER+1
            file.write('\n')
            file.write('rate\t')
            for s,v in self.rates.items():
                file.write(f"{v}\t")
            file.write('\n')

            for name,data in self.nuisances.items():
                kind=data['kind']
                file.write(f'{name} {kind}\t')
                for contrib in self.types.keys():
                    if contrib in data['contribs'].keys():
                        v = data['contribs'][contrib]
                    else:
                        v='-'
                    file.write(f"{v}\t")
                file.write('\n')
                            
                            
                
            
                




            
