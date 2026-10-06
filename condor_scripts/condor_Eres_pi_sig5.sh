#!/bin/bash
ev=100   # using small number of events with many seeds to combine later
ene=$1   # in GeV
boa=$2   # module latter; avaliable: A, B, D, E, G
lay=$3   # later number in HGCAL terms, so layer 34 is the first with scintillators (corresponds to layer # 8 for CMSSW); available 34-47
see=$4   # seed number
sig=5    # sigma in per cent;  5% is nominal for production LY variation
part=211 # particle id; Available options 13 (muon), 211 (pi+), 111 (pi0), 11 (e-)
savedir="" # need to specify where to move generated root file (as they can be quite large, need to make sure to have enough space); e.g your eos space
path2src="" # need to specify where your cmssw release is, path to CMSSW_blabla/src/, somehting like /afs/cern.ch/user/u/username/Simulation/CMSSW_17_0_pre1/src

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${path2src}/Geometry/HGCalTBCommonData/test/python/
cmsenv
echo "Working in /tmp/job$1$2$3$4${sig}/"
mkdir /tmp/job$1$2$3$4${sig}
cd /tmp/job$1$2$3$4${sig}

cmsRun ${path2src}/Geometry/HGCalTBCommonData/test/python/newgun_testHGCalTB24DESYV2_heback_cfg.py seed=${see} energy=${ene} board=${boa} layer=${lay} ev=${ev} sigma=${sig} particle=${part}

echo "Writing output to ${savedir}"
mv *.root ${savedir}
rm -r /tmp/job$1$2$3$4${sig}
