# CMSSW setup
The updated scintillator digitiser for the HGCAL was setup in release CMSSW_17_0_0_pre1, which requires the new version of geometry (v19). 
A patch was added to a new release (CMSSW_20_0_0_patch1), without which there will be issues for the digitiser, however the new digitizer have not been tested with the new release yet. 
However, the patch is implemented in the packages from the updated usercode for the digitiser, so it can be used with release CMSSW_17_0_0_pre1.

## Setting up on lxplus

`source /cvmfs/cms.cern.ch/cmsset_default.sh`  
`cmsrel CMSSW_17_0_0_pre1`  
`cd CMSSW_17_0_0_pre1/src/`  
`cmsenv`  
`scram b -j10`  

Get condor and analysis scripts and the list of required packages:  
`git clone https://github.com/dariasel/cmssw-running-and-analysis.git`  

To add packages from the file:  
`git cms-addpkg -f list_of_packages`    

Get the package updates with realistic sci digitiser:  
`git remote add dariasel git@github.com:dariasel/cmssw`  
`git fetch dariasel sci-digitiser`  
`git checkout sci-digitiser`  
`scram b -j10`  

## Running the simulation

### Manually/locally
Assuming the following already done:  
`cd CMSSW_17_0_0_pre1/src/`  
`source /cvmfs/cms.cern.ch/cmsset_default.sh`  
`cmsenv`  
`scram b -j10`  

Following instructions are for running the newgun_testHGCalTB24DESYV2_heback_cfg.py script, which allows to set number of events, particle energy, particle id, sigma of LY, module type and layer (for the particle gun) and the seed number. The simulation can be run manually:

`cmsRun Geometry/HGCalTBCommonData/test/python/newgun_testHGCalTB24DESYV2_heback_cfg.py seed=1 energy=10 board=E layer=34 ev=1 sigma=5 particle=13`  

With arguments specified. The example above gives the default values for the arguments if not given.
This will run the simulation and produce a root file starting with "gensimdigireco_" and additional info about selected particle, energy, etc. (all arguments).

Condor scripts are also available. If, for whatever reason, you want to run with varying params, but don't want to use condor, the .sh script can be used instead.
To do that:  
`cd condor_scripts`  
`bash condor_Eres_pi_sig5.sh 10 E 34 1`  

Within the script the number of events needs to be changed, also the particle id. The script also has the paths for where to save the data, in order not to fill up the working directory. They need to be specified. 
Same is required if you want to run this on condor.

### Running on condor
`cd condor_scripts`  
`condor_submit run_condor_Eres_pi_sig5`

This will submit jobs as specified in the run_condor script, it should also be set up as needed. The example file is used for the energy resolution studies, so it loops over a list of energies and specified module, layer and seed.

## Analysis

***Important note:***
The current version of the digitiser is set up such that the subsequent Uncalib Rec Hits and Rec Hits are probably not calculated properly. Currently it is expected that the analysis is done with the digis, in the same way as invisioned for the (Uncalib) Rec Hits. The analysis script provided goes through necesary steps to get the digis and the values that are equivalent to (uncalib)RecHits, but not within CMSSW yet.

### Running analysis
Analysis can either be run manually or on condor (or semi-manually with the use of the condor .sh script instead of submitting).
All steps from running the simulation are valid here, except the filenames. Any relevant info in the comments within these scripts.
The analysis relevant scripts: condor_Eres_pi_sig5_analysis.sh, run_condor_Eres_pi_sig5 and sim_analysis_digis-sums-layersums.py

The result of the analysis script is a root file with three trees: digis, sums and layersums.
