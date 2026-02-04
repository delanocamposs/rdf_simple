import subprocess
import ROOT
import json, os
ROOT.gInterpreter.Declare('#include "common/lumiFilter.h"')

def make_jsonHelper(fjson):
    jsondata = json.load(open(fjson))
    runs = []
    firstlumis = []
    lastlumis = []

    for run in jsondata:
        for pair in jsondata[run]:
            runs.append(int(run))
            firstlumis.append(int(pair[0]))
            lastlumis.append(int(pair[1]))
    
    jsonhelper = ROOT.JsonHelper(runs, firstlumis, lastlumis)
    return jsonhelper

def load_meta_data(data):
    dataframe = {}
    #Declare dataframe
    dataframe['Events'] =ROOT.RDataFrame('Events',data['files'])   
    # Apply golden JSON
    if not data['isMC']:
        jsonhelper = None
        if "era" in data.keys():
            if os.path.exists("data/JSON_{}.txt".format(data['era'][:4])): # Non-elegant way to remove pre/postVFP from era definition
                jsonhelper = make_jsonHelper("data/JSON_{}.txt".format(data['era'][:4]))
        if jsonhelper is not None:
            dataframe['Events'] = dataframe['Events'].Define("isGoodLumi", jsonhelper, ["run", "luminosityBlock"])
        else:
            dataframe['Events'] = dataframe['Events'].Define("isGoodLumi", "1")
    else:
        dataframe['Events'] = dataframe['Events'].Define("isGoodLumi", "1")

    #read the HLT string from the sample
    #dataframe=dataframe.DefinePerSample('HLTstring','rdfsampleinfo_.GetS("trigger")')
    dataframe['Events']=dataframe['Events'].Define('HLT_passed',data['trigger']) #need to fix to suppport trigger/amples!!!!
    dataframe['Events']=dataframe['Events'].Define('sample_isMC',str(data['isMC']))
    #get the list of meta keys
    meta_keys =data.keys()
    for key in meta_keys:
        if key=='files' or key=='isMC' or key=='trigger' or key=='jobs' or key=='era' or key == 'customNanoAOD':
            continue;
        dataframe['Events']=dataframe['Events'].Define('sample_'+key,str(data[key]))
    dataframe['Runs'] = ROOT.RDataFrame("Runs", data['files'])
    dataframe['isMC'] = data['isMC']
    dataframe['customNanoAOD'] = data['customNanoAOD']
    dataframe['isScouting']=data['isScouting']
    return dataframe


def loadSample(info,locator='root://cms-xrd-global.cern.ch//'):
    meta={}
    files=[]
    print(info['dataset'])
    if 'local:' in info['dataset']:
        files=info['dataset'].split('local:')[1].split(',')
    else:
        p = subprocess.Popen('/cvmfs/cms.cern.ch/common/dasgoclient -dasmaps=./ -query="file dataset={}"'.format(info['dataset']), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            rootfile = line.decode('ASCII').split('\n')[0]
            if '.root' not in rootfile:
                continue
            files.append(locator+rootfile)
        retval = p.wait()
    triggerStr=[]
    for t in info['triggers']:
        triggerStr.append('{}==1'.format(t))
    for t in info['veto_triggers']:
        triggerStr.append('{}==0'.format(t))
    triggerDecision = '&&'.join(triggerStr)
    meta['trigger']=triggerDecision
    if not ('SIM' in info['dataset'] or 'local' in info['dataset']):
        meta['isMC']=0
    else:
        meta['isMC']=1
    for tag, data in info.items():
        if tag not in ['dataset','triggers','veto_triggers']:
            meta[tag]=data
    meta['files'] = files        
    if 'era' in info.keys():
        meta['era'] = info['era']
    else:
        meta['era'] = '2018' # Assume 2018 if not included
    if 'customNanoAOD' in info.keys():
        meta['customNanoAOD'] = info['customNanoAOD']
    else:
        meta['customNanoAOD'] = False
    if 'isScouting' in info.keys():
        meta['isScouting']=info['isScouting']
    else:
        meta['isScouting']=0
    print(meta['isScouting'])
    print(type(meta['isScouting']))
    return meta


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


def createDataSet(sample,splitFactor=0,processPart=0,locator='root://cms-xrd-global.cern.ch//'):
    spec = loadSample(sample,locator)
    if splitFactor==0:
        return spec
    else:
        #split the sample
        fileNames = spec['files']
        if len(fileNames)<splitFactor:
            print('Number of files smaller than the number of job splits.you should only split in  {} parts.Job will not generate output'.format(len(fileNames)))
        splitFileNames = list(split(fileNames,splitFactor))[processPart]
        print('split file names')
        for s in splitFileNames:
            print(s)
        spec['files']=splitFileNames
        return spec
        
