samples_signal_2022 = {}
#slight version changes from nominal dataset naming for select mass,lifetime points. not sure why MC&I did it
datasets={
    (55, 1000, "2022preEE"):"/GluGluHTo2LongLivedTo4G-MH125_MFF55-ctau1000mm_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM", 
    (55, 1000, "2022postEE"): "/GluGluHTo2LongLivedTo4G-MH125_MFF55-ctau1000mm_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM",

    (20,50,"2022preEE"):"/GluGluHTo2LongLivedTo4G-ctau-50mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM",
    (20,50,"2022postEE"):"/GluGluHTo2LongLivedTo4G-ctau-50mm-MFF-20-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v4/NANOAODSIM",

    (30, 0, "2022preEE"):"/GluGluHTo2LongLivedTo4G-ctau-0mm-MFF-30-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v4/NANOAODSIM",
    (30, 0, "2022postEE"):"/GluGluHTo2LongLivedTo4G-ctau-0mm-MFF-30-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM",

    (30, 100, "2022preEE"):"/GluGluHTo2LongLivedTo4G-ctau-100mm-MFF-30-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v4/NANOAODSIM",
    (30, 100, "2022postEE"):"/GluGluHTo2LongLivedTo4G-ctau-100mm-MFF-30-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM"}

for mass in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100, 1000]:
        if (mass==55 and ct==1000) or (mass==20 and ct==50) or (mass==30 and ct==0) or (mass==30 and ct==100):
            samples_signal_2022[f"ggH4g_M{mass}_ctau{ct}_2022preEE"] = {
                    'dataset': datasets[(mass,ct,"2022preEE")],
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2022preEE',
                    'sigma': 1.0,
                    'customNanoAOD': True}

            samples_signal_2022[f"ggH4g_M{mass}_ctau{ct}_2022postEE"] = {
                    'dataset': datasets[(mass,ct,"2022postEE")],
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2022postEE',
                    'sigma': 1.0,
                    'customNanoAOD': True}
        
        else:
            samples_signal_2022[f"ggH4g_M{mass}_ctau{ct}_2022preEE"] = {
                    'dataset': f"/GluGluHTo2LongLivedTo4G-ctau-{ct}mm-MFF-{mass}-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2022preEE',
                    'sigma': 1.0,
                    'customNanoAOD': True}

            samples_signal_2022[f"ggH4g_M{mass}_ctau{ct}_2022postEE"] = {
                    'dataset': f"/GluGluHTo2LongLivedTo4G-ctau-{ct}mm-MFF-{mass}-MH-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM",
                    'triggers': ['HLT_DoublePhoton33_CaloIdL'],
                    'isMC':1,
                    'jobs':1,
                    'veto_triggers': [],
                    'era': '2022postEE',
                    'sigma': 1.0,
                    'customNanoAOD': True}

