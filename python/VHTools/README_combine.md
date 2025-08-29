Setup Instructions
==================
1. install combine

cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.2
scramv1 b clean; scramv1 b

2. Install combine tools

cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b

3. Tyler's original scripts

git clone https://github.com/Tyler-Lam/DDPcombine.git python/

4. Replace the CombineHarvester/CombineTools/python/plotting.py script with the one here (only a few minor changes, mainly to axis ranges and legend labeling)