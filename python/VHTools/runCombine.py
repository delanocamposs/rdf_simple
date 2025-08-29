import os, math, glob
import ROOT

import optparse
parser = optparse.OptionParser()
parser.add_option("-d", "--date", dest="date", default="08_12_24",help="date for plot directory")
parser.add_option("-y","--year",dest="year",default="2018",help="2016 or 2017 or 2018 (or Run2)")
parser.add_option("-b", "--blind", dest = "blind", action="store_true", default = False, help="Use option --run blind for combine")
parser.add_option("-m", "--mass", dest ="mass")
(options,args) = parser.parse_args()

dir_out = "/uscms/home/gfavila/nobackup/rdf_simple_/" # Change this to wherever you keep your datacards
#dir_out = "/home/tyler/DDP/rdf_simple/"
if not os.path.isdir(dir_out+"datacards_{}".format(options.date)):
    os.mkdir(dir_out+"datacards_{}".format(options.date))
os.chdir(dir_out+"datacards_{}".format(options.date))


#m = options.mass
masses = [options.mass]
#ctaus = [0, 3, 10, 14, 20, 32, 50, 70, 100, 316, 1000]

ctaus = [0,10,20,50,100,1000]

V = ['Z','W']
L = ['ELE','MU']
#V = ['Z']

modelIndependent = False


if modelIndependent: 
    sfx = "_modelIndependent"
else:
    sfx = ""
print(sfx)

years = []
if options.year == "Run2":
    years = ['2016', '2017', '2018']
else:
    years.append(options.year)

#signal strength 
#rMax of 1 = what is input into the datacard
rMax = 5

# First combine bins for the individual categories
for year in years:
    for m in masses:
        for ct in ctaus:
            for v in V:
                for l in L:
                    name = "{}H_{}_m{}_ctau{}_{}{}".format(v,l,m,ct,year,sfx) #NOTE: model independent
                    cmd = "combineCards.py "
                    files = glob.glob("datacard_{}_bin*.txt".format(name))
                    for f in files:
                        cmd += "{} ".format(f)
                    cmd += " > datacard_{}.txt".format(name)
                    os.system(cmd)

# If full run2, then combine all years
if options.year == "Run2": 
    rMax = 2.0 # Need lower rMax due to higher sensitivity
    for m in masses:
        for ct in ctaus:
            for v in V:
                for l in L:
                    name = "{}H_{}_m{}_ctau{}".format(v,l,m,ct)
                    cmd = "combineCards.py "
                    for y in years:
                        cmd += "datacard_{}_{}{}.txt ".format(name,y,sfx)  #NOTE: model independent
                    cmd += " > datacard_{}_Run2{}.txt".format(name,sfx)    #NOTE: model independent
                    os.system(cmd)

for m in masses:
    for ct in ctaus:
        cmd_combine = "combineCards.py " # Command to combine all categories into 1 card
        for v in V:
            cmd_v = "combineCards.py " # Command to combine e+mu categories for a given V type
            for l in L:
                # Run combine to get limits for each individual category
                print ("Running combine for {}->{} m{} ctau{}".format(v,l,m,ct))
                name = "{}H_{}_m{}_ctau{}_{}{}".format(v,l,m,ct,options.year,sfx) #NOTE: model independent
                cmd = "combine -M AsymptoticLimits datacard_{}.txt -m 125 -n .{} --rMax {}".format(name, name, rMax)
                                                    #datacard_WH_ELE_m15_ctau0_2016_modelIndependent.txt
                if options.blind:
                    cmd += " --run blind"
                os.system(cmd)

                cmd_v += "datacard_{}.txt ".format(name)

            # Run combine on combined e+mu datacard
            name = "{}H_m{}_ctau{}_{}{}".format(v,m,ct,options.year,sfx) #NOTE: model independent
            cmd_v += "> datacard_{}.txt".format(name)
            os.system(cmd_v)

            print ("Running combine for all {}H m{} ctau{}".format(v, m, ct))
            cmd = "combine -M AsymptoticLimits datacard_{}.txt -m 125 -n .{} --rMax {}".format(name,name,rMax)
            if options.blind:
                cmd += " --run blind"
            os.system(cmd)

            cmd_combine += "datacard_{}.txt ".format(name)
        if len(V) == 2:
            # combineCards all V types and run combine
            name = "VH_m{}_ctau{}_{}{}".format(m,ct,options.year,sfx) #NOTE: model independent
            cmd_combine += "> datacard_{}.txt".format(name)
            os.system(cmd_combine)

            print ("Running combine for all VH m{} ctau{}".format(m,ct))
            cmd = "combine -M AsymptoticLimits datacard_{}.txt -m 125 -n .{} --rMax {}".format(name, name, rMax)
            if options.blind:
                cmd += " --run blind"
            os.system(cmd)
