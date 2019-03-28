cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
root -l -b -q 'GetBaseLineEvents.C+("data","data","/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Data/","mm",true)' 
