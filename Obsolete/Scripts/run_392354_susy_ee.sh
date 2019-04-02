cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
root -l -b -q 'GetBaseLineEvents.C+("392354","susy","root://eosatlas//eos/atlas/user/r/rshang/ew_susy_highpt/","ee",false)' 
