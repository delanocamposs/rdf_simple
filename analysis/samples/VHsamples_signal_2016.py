#PREVFP scale:0.54
#postVFP scale:0.46

samples_signal_2016 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [15, 20, 30, 40, 50, 55]:
        for ct in [0, 10, 20, 50, 100, 1000]:
            for br in ['2G2Q']:
                samples_signal_2016['{}H{}_M{}_ctau{}_2016{}'.format(v,br,m,ct,'preVFP')] = {
                    'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/{}H_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}cm.root'.format('preVFP',v,br, m, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2016preVFP',
                    'sigma': 0.54,
                    'customNanoAOD': True}
                samples_signal_2016['{}H{}_M{}_ctau{}_2016{}'.format(v,br,m,ct,'postVFP')] = {
                    'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/{}H_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}cm.root'.format('postVFP',v,br, m, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2016postVFP',
                    'sigma': 0.46,
                    'customNanoAOD': True}
                
#ttH samples
for m in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100]:
        for br in ['2G2Q']:
            samples_signal_2016['ttH{}_M{}_ctau{}_2016{}'.format(br,m,ct,'preVFP')] = {
                'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm_v2.root'.format('preVFP', br, m, ct),
                'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2016preVFP',
                'sigma': 0.54,
                'customNanoAOD': True}
            samples_signal_2016['ttH{}_M{}_ctau{}_2016{}'.format(br,m,ct,'postVFP')] = {
                'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm_v2.root'.format('postVFP', br, m, ct),
                'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2016postVFP',
                'sigma': 0.46,
                'customNanoAOD': True}
            
for m in [15, 20, 30, 40, 50, 55]:
    for ct in [1000]:
        for br in ['2G2Q']:
            samples_signal_2016['ttH{}_M{}_ctau{}_2016{}'.format(br,m,ct,'preVFP')] = {
                'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm.root'.format('preVFP', br, m, ct),
                'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2016preVFP',
                'sigma': 0.54,
                'customNanoAOD': True}
            samples_signal_2016['ttH{}_M{}_ctau{}_2016{}'.format(br,m,ct,'postVFP')] = {
                'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm.root'.format('postVFP', br, m, ct),
                'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2016postVFP',
                'sigma': 0.46,
                'customNanoAOD': True}


