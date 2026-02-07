samples_signal_2018 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [15, 20, 30, 40, 50, 55]:
        for ct in [0, 10, 20, 50, 100, 1000]:
            for br in ['2G2Q']:
                samples_signal_2018["{}H{}_M{}_ctau{}_2018".format(v, br, m, ct)] = {
                    'dataset': "local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2018/{}H_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}cm.root".format(v, br, m, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2018',
                    'sigma': 1.0,
                    'customNanoAOD': True}

# ttH samples:
for m in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100,1000]:
        for br in ['2G2Q']:
            samples_signal_2018["ttH{}_M{}_ctau{}_2018".format(br, m, ct)] = {
                'dataset': "local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2018/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm_v2.root".format(br, m, ct),
                'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2018',
                'sigma': 1.0,
                'customNanoAOD': True}

## ttH samples:
#or m in [15, 20, 30, 40, 50, 55]:
 #   for ct in [1000]:
 #       for br in ['2G2Q']:
 #           samples_signal_2018["ttH{}_M{}_ctau{}_UL18".format(br, m, ct)] = {
 #               'dataset': "local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2018/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm.root".format(br, m, ct),
 #               'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
 #               'veto_triggers': [],
 #               'era': '2018',
 #               'sigma': 1.0,
 #               'customNanoAOD': True}
