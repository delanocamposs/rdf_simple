from optparse import OptionParser
import ROOT
import os
import time
import importlib
parser = OptionParser()

parser.add_option("-a", "--analysis", dest="analysis",
                  help="analysis to run",type='str',default='ggH')

parser.add_option("-o", "--eosdir", dest="eos",default='root://cmseos.fnal.gov//store/user/dacampos/analysis/MC_processed/',
                  help="EOS output Directory")
parser.add_option("-d", "--donotsubmit", dest="nosubmit",default=0,type=int,
                  help="Do not submit on condor, just show jobs")
parser.add_option("-v", "--verbose", dest="verbose",default=0,type=int,
                  help="Verbose")

(options, args) = parser.parse_args()

samp = importlib.import_module('analysis.{}samples'.format(options.analysis))


datasets= args

# Only keep selected datasets if args provided
if datasets:
    samp.samples = {k: v for k, v in samp.samples.items() if k in datasets}


def check_file(filename):
    f=0
    try:
        f=ROOT.TFile.Open(filename)
    except:
        return False
    if f.IsZombie()>0:
        return False
    else:
        if f.GetNkeys()>0:
            return True
        else:
            return False

    
    
os.system('tar --overwrite --exclude=".git" --exclude="common/__pycache__" --exclude="toys/*" --exclude="local_samples/*" --exclude="pandas_dataframes/*" --exclude="*.ipynb" --exclude="datacard*" --exclude="*.pdf" --exclude="condor_output_files/*" --exclude="DDP/*" --exclude="*.png" --exclude="plots/*" --exclude="higgsCombine*" --exclude="analysis/__pycache__"  -czvf /tmp/sandbox.tar.gz .')
os.system('mv /tmp/sandbox.tar.gz .')


for d,info in samp.samples.items():    
    if (not ('jobs' in info)):
        filename="{redirect}/{sample}_{sample}.root".format(redirect=options.eos,sample=d)
        if check_file(filename):
            if options.verbose:
                print ('Filename exists = {}.Will not submit'.format(filename)) 
            continue
        shell="""#!/bin/sh
        echo starting script
        tar -xzvf sandbox.tar.gz
        mkdir .dasmaps
        source ./setup.sh
        python3 run.py -a {analysis} {dataset}
        echo python done
        
        OUTDIR={eos}
        echo "xrdcp output for condor to "
        echo $OUTDIR
        for FILE in *.root
        do
        xrdcp -f ${{FILE}} ${{OUTDIR}}/{dataset}_${{FILE}} 2>&1
        XRDEXIT=$?
        if [[ $XRDEXIT -ne 0 ]]; then
        rm *.root
        echo "exit code $XRDEXIT, failure in xrdcp"
        exit $XRDEXIT
        fi
        rm ${{FILE}}
        done
        """.format(dataset=d,eos=options.eos,analysis=options.analysis)
        #print(shell)
        f=open("{dataset}_condor.sh".format(dataset=d),"w")
        f.write(shell)
        f.close()
        os.system('chmod +x {dataset}_condor.sh'.format(dataset=d))
        
        condor="""
        universe = vanilla
        Executable = {dataset}_condor.sh
        should_transfer_files = YES
        Transfer_Input_Files = sandbox.tar.gz
        when_to_transfer_output = ON_EXIT
        Output = condor_{dataset}_$(Cluster)_$(Process).stdout
        Error = condor_{dataset}_$(Cluster)_$(Process).stderr
        Log = condor_{dataset}_$(Cluster)_$(Process).log
        Queue 1
        """.format(dataset=d)
        f=open("{dataset}_condor.jdl".format(dataset=d),"w")
        f.write(condor)
        f.close()
        if options.nosubmit:
            print('condor_submit {dataset}_condor.jdl'.format(dataset=d))
        else:
            os.system('condor_submit {dataset}_condor.jdl'.format(dataset=d))
        time.sleep(0.1)
    else:
        for i in range(0,info['jobs']):
            filename="{redirect}/{sample}_{i}_{sample}.root".format(redirect=options.eos,i=i,sample=d)
            if check_file(filename):
                if options.verbose:
                    print ('Filename exists = {}.Will not submit'.format(filename)) 
                continue

            shell="""#!/bin/sh
            echo starting script
            tar -xzvf sandbox.tar.gz
            mkdir .dasmaps
            source ./setup.sh
            python3 run.py {dataset} -s {splitFactor} -p {part} -a {analysis}
            echo python done
        
            OUTDIR={eos}
            echo "xrdcp output for condor to "
            echo $OUTDIR
            for FILE in *.root
            do
            xrdcp -f ${{FILE}} ${{OUTDIR}}/{dataset}_{part}_${{FILE}} 2>&1
            XRDEXIT=$?
            if [[ $XRDEXIT -ne 0 ]]; then
            rm *.root
            echo "exit code $XRDEXIT, failure in xrdcp"
            exit $XRDEXIT
            fi
            rm ${{FILE}}
            done
            """.format(dataset=d,eos=options.eos,splitFactor=info['jobs'],part=i,analysis=options.analysis)
         #   print(shell)
            f=open("{dataset}_{part}_condor.sh".format(dataset=d,part=i),"w")
            f.write(shell)
            f.close()
            os.system('chmod +x {dataset}_{part}_condor.sh'.format(dataset=d,part=i))
            
            condor="""
            universe = vanilla
            Executable = {dataset}_{part}_condor.sh
            should_transfer_files = YES
            Transfer_Input_Files = sandbox.tar.gz
            when_to_transfer_output = ON_EXIT
            Output = condor_{dataset}_{part}_$(Cluster)_$(Process).stdout
            Error = condor_{dataset}_{part}_$(Cluster)_$(Process).stderr
            Log = condor_{dataset}_{part}_$(Cluster)_$(Process).log
            Queue 1
            """.format(dataset=d,part=i)
            f=open("{dataset}_{part}_condor.jdl".format(dataset=d,part=i),"w")
            f.write(condor)
            f.close()
            if options.nosubmit:
                print('condor_submit {dataset}_{part}_condor.jdl'.format(dataset=d,part=i))
            else:
                os.system('condor_submit {dataset}_{part}_condor.jdl'.format(dataset=d,part=i))
            time.sleep(0.1)
            
