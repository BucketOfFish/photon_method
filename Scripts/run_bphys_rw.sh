cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root
root -l -b -q 'GetJpsiReweighting.C+("ee",0,"2016")' 
root -l -b -q 'GetJpsiReweighting.C+("mm",0,"2016")' 
root -l -b -q 'GetJpsiReweighting.C+("ee",1,"2016")' 
root -l -b -q 'GetJpsiReweighting.C+("mm",1,"2016")' 
root -l -b -q 'GetJpsiReweighting.C+("ee",2,"2016")' 
root -l -b -q 'GetJpsiReweighting.C+("mm",2,"2016")' 
