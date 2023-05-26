import sys
import subprocess
from mkShapesRDF.processor.framework.Steps_cfg import Steps
from mkShapesRDF.processor.framework.Productions_cfg import Productions
import os
from pathlib import Path


def getFiles(das, folder, process, redirector="root://cms-xrd-global.cern.ch/"):
    files = []
    if das:
        procString = f'dasgoclient --query="file dataset={process}"'
        proc = subprocess.Popen(
            procString, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = proc.communicate()
        out = out.decode("utf-8")
        out = out.split("\n")
        out = list(filter(lambda k: k.strip() != "", out))
        print(out, len(out))
        err = err.decode("utf-8")
        if len(err) != 0:
            print("There were some errors in retriexing file:")
            print(err)
            sys.exit()
        files = list(map(lambda k: redirector + k, out))
        return files[1:]

    elif redirector != "":
        # not das and use redirector provided 
        import fnmatch

        files = []
        listOfFiles = []
        print("Need to query for files for folder", folder)

        proc = subprocess.Popen(
            f"xrdfs {redirector} ls {folder}", shell=True, stdout=subprocess.PIPE
        )
        listOfFiles = proc.communicate()[0].decode("utf-8").split("\n")
        files = list(
            filter(lambda k: fnmatch.fnmatch(k, folder + process), listOfFiles)
        )
        files = list(map(lambda k: redirector + k, files))
        return files[:]

    else:
        import glob
        # not das and no redirector provided -> use direct access
        files = []
        listOfFiles = []
        print("Need to query for files for folder", folder)
        files = glob.glob(folder + '/' + process)
        return files[:]

def getFiles_cfg(sampleName):
    if inputFolder == "":
        return  {
                "das": True,
                "folder": "",
                'process': Samples[sampleName]["nanoAOD"],
            }
    else:
        return {
            "das": False,
            "redirector": redirector,
            "folder": inputFolder,
            "process": "*" + sampleName + "*.root",
        }


condorDir = "/".join(os.path.abspath(".").split("/")[:-1]) + "/condor"
eosDir = '/eos/user/l/lviliani/mkShapesRDF_nanoAOD'

procName = "Run2018_UL2018_nAODv9_Full2018v9"
sampleName = "DoubleMuon_Run2018A-UL2018-v1"
step = "DATAl1loose2018v9"

procName = "Summer20UL18_106x_nAODv9_Full2018v9"
#sampleName = "DYJetsToLL_M-50"
sampleName = "TTTo2L2Nu"
step = "MCFull2018v9_JES_l2Kin"
#sampleName = "EWKZ2Jets_ZToLL_M-50_MJJ-120"

submit = False

default = False
if default:
    inputFolder = ""
    redirector = ""
    useProxy = True
    maxFilesPerJob = 1
    limitFiles = 10

else:
    inputFolder = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9//MCl1loose2018v9__MCCorr2018v9NoJERInHorn__l2tightOR2018v9/"
    redirector = ""
    useProxy = False
    maxFilesPerJob = 1
    limitFiles = 1



fPy = """
import ROOT
ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
"""
fPy += "from mkShapesRDF.processor.framework.mRDF import mRDF\n"

if Productions[procName]["isData"]:
    fPy += 'lumiFile = "' + os.path.abspath(Productions[procName]["jsonFile"]) + '"\n'
else:
    xsFile = Productions[procName]["xsFile"]
    with open(xsFile) as file:
        exec(file.read())
    fPy += "xs_db = " + str(xs_db) + "\n"

with open(Productions[procName]["samples"]) as file:
    exec(file.read())


files = getFiles(**getFiles_cfg(sampleName))[:limitFiles]

if len(files) == 0:
    print("No files found for", sampleName)
    sys.exit()

print(files[0])


fSh = ''
with open('../../../start.sh') as file:
    fSh += file.read()

if useProxy:
    cmd = "voms-proxy-info"
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    out, err = proc.communicate()

    if "Proxy not found" in err.decode("utf-8"):
        print("WARNING: No GRID proxy -> Get one first with:")
        print("voms-proxy-init -voms cms -rfc --valid 168:0")
        sys.exit()

    proxypath = " xxx "
    for line in out.decode("utf-8").split("\n"):
        if "path" in line:
            proxypath = line.split(":")[1].strip()
    print(proxypath)

    os.system("cp " + proxypath + " " + os.environ["HOME"] + "/.proxy")

    fSh += f"export X509_USER_PROXY={os.environ['HOME']}/.proxy\n"



fSh += "time python script.py\n"

folderPathEos = eosDir + '/' + procName + "/" + step 
Path(folderPathEos).mkdir(parents=True, exist_ok=True)

fSh += 'cp output.root ' + folderPathEos + "/" + "${1}.root\n"
fSh += 'rm output.root\n'

print(fSh)

jobDir = condorDir + "/" + procName + "/" + step + "/" 
Path(jobDir).mkdir(parents=True, exist_ok=True)

with open(jobDir + "run.sh", "w") as file:
    file.write(fSh)

os.system("chmod +x " + jobDir + "run.sh")

from math import ceil
nParts = ceil(len(files)/maxFilesPerJob)

fJdl = '''\
universe = vanilla
executable = run.sh
arguments = $(Folder)

should_transfer_files = YES
transfer_input_files = $(Folder)/script.py

output = $(Folder)/out.txt
error  = $(Folder)/err.txt
log    = $(Folder)/log.txt

request_cpus   = 1
+JobFlavour = "workday"

queue 1 Folder in ''' + ', '.join(list(map(lambda k: sampleName + '__part' + str(k), list(range(nParts)))))

with open(jobDir + '/submit.jdl', 'w') as f:
    f.write(fJdl)



frameworkPath = "/".join(os.path.abspath(".").split("/")[:-2])

fPy += "sampleName = '" + sampleName + "'\n"

fPy += "files = RPLME_FILES\n"
fPy += f"ROOT.gInterpreter.Declare('#include \"{frameworkPath}/include/headers.hh\"')\n"

fPy += "df = mRDF()\n"
fPy += 'df = df.readRDF("Events", files)\n'

fPy += "values = []\n"


def addDeclareLines(step):
    global fPy
    if step not in Steps.keys():
        print(f"Error, step {step} not found in Steps")
        sys.exit()
    if Steps[step]["isChain"]:
        for subtarget in Steps[step]["subTargets"]:
            addDeclareLines(subtarget)
    elif "declare" in Steps[step].keys():
        fPy += "from " + Steps[step]["import"] + " import *\n"
        fPy += Steps[step]["declare"] + "\n"
        fPy += "module = " + Steps[step]["module"] + "\n"
        fPy += "df = module.run(df, values) \n"


addDeclareLines(step)

fPy += """
snapshots = []
for val in values:
    if "snapshot" == val[0]:
        snapshots.append(val[1])

if len(snapshots) != 0:
    ROOT.RDF.RunGraphs(snapshots)

def sciNot(value):
    # scientific notation
    return "{:.3e}".format(value)

data = []
for val in values:
    if "snapshot" == val[0]:
        continue
    if "list" in str(type(val)):
        if str(type(val[0])) == "<class 'function'>":
            data.append(val[0](*val[1:]))
        else:
            data.append([val[1], sciNot(val[0].GetValue())])
    else:
        data.append("", sciNot(val.GetValue()))

from tabulate import tabulate

print(tabulate(data, headers=["desc.", "value"]))
"""



fPy = fPy.replace('RPLME_FW', frameworkPath)

if 'RPLME_genEventSumw' in fPy:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.EnableImplicitMT()
    df2 = ROOT.RDataFrame("Runs", files)
    genEventSumw = df2.Sum("genEventSumw").GetValue()
    fPy = fPy.replace('RPLME_genEventSumw', str(genEventSumw))

for part in range(nParts):
    _files = files[ part * maxFilesPerJob : (part +1) * maxFilesPerJob]
    _fPy = fPy.replace('RPLME_FILES', str(_files))


    jobDirPart = jobDir + sampleName + "__part" + str(part) + "/"
    Path(jobDirPart).mkdir(parents=True, exist_ok=True)

    with open(jobDirPart + "/script.py", "w") as f:
        f.write(_fPy)


if submit:
    proc = subprocess.Popen(f'cd {jobDir}; condor_submit submit.jdl;', shell=True)
    proc.wait()
