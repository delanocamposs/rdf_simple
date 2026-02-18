samples_signal_2024 = {}
for mass in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100, 1000]:
        samples_signal_2024[f"ggH4g_M{mass}_ctau{ct}_2024"] = {
                    'dataset': f"/GluGluHTo2LongLivedTo4G_Par-ctau-{ct}mm-MFF-{mass}-MH-125_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM",
                    'jobs':1,
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'veto_triggers': [],
                    'era': '2024',
                    'sigma': 1.0,
                    'customNanoAOD': True}

