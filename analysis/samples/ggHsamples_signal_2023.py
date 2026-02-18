samples_signal_2023 = {}
#slight version changes from nominal dataset naming for select mass,lifetime points. not sure why MC&I did it
datasets={
    (20, 50, "2023preBPix"):"/GluGluHTo2LongLivedTo4G-ctau-50mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM",
    (20, 50, "2023postBPix"): "/GluGluHTo2LongLivedTo4G-ctau-50mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v4/NANOAODSIM",

    (20,1000,"2023preBPix"):"/GluGluHTo2LongLivedTo4G-ctau-1000mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM",
    (20,1000,"2023postBPix"):"/GluGluHTo2LongLivedTo4G-ctau-1000mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v4/NANOAODSIM",

    (40, 20, "2023preBPix"):"/GluGluHTo2LongLivedTo4G-ctau-20mm-MFF-40-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v4/NANOAODSIM",
    (40, 20, "2023postBPix"):"/GluGluHTo2LongLivedTo4G-ctau-20mm-MFF-40-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v3/NANOAODSIM",

    (55, 0, "2023preBPix"):"/GluGluHTo2LongLivedTo4G-ctau-0mm-MFF-55-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v4/NANOAODSIM",
    (55, 0, "2023postBPix"):"/GluGluHTo2LongLivedTo4G-ctau-0mm-MFF-55-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v3/NANOAODSIM",

    (55, 20, "2023preBPix"):"/GluGluHTo2LongLivedTo4G-ctau-20mm-MFF-55-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v4/NANOAODSIM",
    (55, 20, "2023postBPix"):"/GluGluHTo2LongLivedTo4G-ctau-20mm-MFF-55-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v3/NANOAODSIM", 

    (55, 1000, "2023preBPix"):"/GluGluHTo2LongLivedTo4G-MH125_MFF55-ctau1000mm_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM",
    (55, 1000, "2023postBPix"):"/GluGluHTo2LongLivedTo4G-MH125_MFF55-ctau1000mm_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v3/NANOAODSIM"}



for mass in [15,20,30,40,50,55]:
    for ct in [0,10,20,50,100,1000]:
        if (mass==20 and ct==50) or (mass==20 and ct==1000) or (mass==40 and ct==20) or (mass==55 and ct==0) or (mass==55 and ct==20) or (mass==55 and ct==1000):
            samples_signal_2023[f"ggH4g_M{mass}_ctau{ct}_2023preBPix"] = {
                    'dataset': datasets[(mass,ct,"2023preBPix")],
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2023preBPix',
                    'sigma': 1.0,
                    'customNanoAOD': True}

            samples_signal_2023[f"ggH4g_M{mass}_ctau{ct}_2023postBPix"] = {
                    'dataset': datasets[(mass,ct,"2023postBPix")],
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2023postBPix',
                    'sigma': 1.0,
                    'customNanoAOD': True}
        else:
            samples_signal_2023[f"ggH4g_M{mass}_ctau{ct}_2023preBPix"] = {
                    'dataset': f"/GluGluHTo2LongLivedTo4G-ctau-{ct}mm-MFF-{mass}-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2023preBPix',
                    'sigma': 1.0,
                    'customNanoAOD': True}

            samples_signal_2023[f"ggH4g_M{mass}_ctau{ct}_2023postBPix"] = {
                    'dataset': f"/GluGluHTo2LongLivedTo4G-ctau-{ct}mm-MFF-{mass}-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v3/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2023postBPix',
                    'sigma': 1.0,
                    'customNanoAOD': True}

