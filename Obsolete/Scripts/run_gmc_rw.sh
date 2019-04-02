cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root
root -l -b -q 'GetPhotonReweighting.C+("ee",0,"2016")' 
root -l -b -q 'GetPhotonReweighting.C+("mm",0,"2016")' 
