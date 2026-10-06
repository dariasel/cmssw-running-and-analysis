#!/bin/bash
ene=${1} # mandatory for running the analysis, the script needs to know which energy to analyse
sig=${2} # could be useful for the specification of data location; not mandatory for the analysis script, since this info is parse from the filename
boa=${3} # could be useful for the specification of data location; not mandatory for the analysis script, since this info is parse from the filename
#part=211 # no need to specify, since the provided analysis script parses this out of the filename

path2src="" # need to specify path to cmssw release src, since the dataformats are required for analysis, something like /afs/cern.ch/user/u/username/Simulation/CMSSW_17_0_0_pre1/src

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${path2src}
cmsenv
echo "Working in /tmp/job$1${sig}/"
echo "Working on Ene $1 Sigma $2"
mkdir /tmp/job$1${sig}
cd /tmp/job$1${sig}

# need to specify
dirpy="/path/to/analysis/script/"sim_analysis_digis-sums-layersums.py
dirdat="/path/to/location/of/cmssw/root/file/created/after/cmsRun/"
path2save="/path/to/save/analysis/root/file/"
python3 ${dirpy} ${dirdat} ${ene}
mv *.root ${path2save}
rm -r /tmp/job$1${sig}
