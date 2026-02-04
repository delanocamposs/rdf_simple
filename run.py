import ROOT
import importlib
#Enable multi-threading
#ROOT.ROOT.EnableImplicitMT()
#ROOT.ROOT.DisableImplicitMT()
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
from common.pyhelpers import *
from optparse import OptionParser





parser = OptionParser()
parser.add_option("-a", "--analysis", dest="analysis",
                  help="Split the sample in N sub samples",type='str',default='ggH')
parser.add_option("-s", "--splitFactor", dest="splitFactor",
                  help="Split the sample in N sub samples",type='int',default=0)
parser.add_option("-p", "--processPart", dest="processPart",
                  help="process the i-th sample if you decide to split ",type='int',default=0)
(options, args) = parser.parse_args()



an = importlib.import_module('analysis.{}'.format(options.analysis))
samp = importlib.import_module('analysis.{}samples'.format(options.analysis))

Data_redirectors={
"Europe": "root://xrootd-cms.infn.it//",
"US":     "root://cmsxrootd.fnal.gov//",
"Global": "root://cms-xrd-global.cern.ch//",
"Asia":   "root://xrootd-cms-kr.kisti.re.kr//"
}
actions=[]
for sample in args:
    data=createDataSet(samp.samples[sample],options.splitFactor,options.processPart,Data_redirectors["US"])
    actions.extend(an.analysis(data,sample))
ROOT.RDF.RunGraphs(actions)

