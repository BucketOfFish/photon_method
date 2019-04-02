cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
root -l -b -q 'GetBaseLineEvents.C+("361406","zjets","/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/zmc/","mm",false)' 
