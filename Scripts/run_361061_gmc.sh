cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
root -l -b -q 'GetPhotonEvents.C+("361061","gmc","/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gmc/",0)' 
