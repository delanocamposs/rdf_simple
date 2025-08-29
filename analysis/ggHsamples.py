samples = {}

# Import samples for all eras
for era in ['2016', '2017', '2018', '2022','2023','2024']:
    for s in ["data"]:
        cmd = """
from analysis.samples.ggHsamples_{s}_{era} import *
for s in samples_{s}_{era}:
   if s not in samples:
        samples[s] = samples_{s}_{era}[s]
""".format(s=s, era=era)
        exec(cmd)
