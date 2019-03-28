cd ../Root 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup ROOT
root -l -b -q 'GetBPhysicsEvents.C+("bdata","bphys","/eos/atlas/user/r/rshang/bphy_test/bdata/",1)' 
