samples_signal_2017 = {}
for mass in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100, 1000]:
        samples_signal_2017[f"ggH4g_M{mass}_ctau{ct}_2017"] = {
                    'dataset': f"/GluGluH_HTo2LongLivedTo4G_MH-125_MFF-{mass}_ctau-{ct}mm_TuneCP5_13TeV_powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
                    'triggers': ['HLT_TriplePhoton_20_20_20_CaloIdLV2'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2017',
                    'sigma': 1.0,
                    'customNanoAOD': True}

