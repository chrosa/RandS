#!/bin/bash 
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/TKinFitter
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
cd /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8
lsetup rcsetup
cd /nfs/dust/atlas/user/csander/RandS/Output/jobdir/527
ln -s /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/TKinFitter/libkinfitter.so .
/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/RunRandS.x
