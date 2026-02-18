samples_signal_2016 = {}
for mass in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100, 1000]:
        samples_signal_2016[f"ggH4g_M{mass}_ctau{ct}_2016preVFP"] = {
                    'dataset': f"/GluGluH_HTo2LongLivedTo4G_MH-125_MFF-{mass}_ctau-{ct}mm_TuneCP5_13TeV_powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton60'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2016preVFP',
                    'sigma': 1.0,
                    'customNanoAOD': True}

        samples_signal_2016[f"ggH4g_M{mass}_ctau{ct}_2016postVFP"] = {
                    'dataset': f"/GluGluH_HTo2LongLivedTo4G_MH-125_MFF-{mass}_ctau-{ct}mm_TuneCP5_13TeV_powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton60'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2016postVFP',
                    'sigma': 1.0,
                    'customNanoAOD': True}

