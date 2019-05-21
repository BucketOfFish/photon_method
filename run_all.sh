# Setup environment

lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"

# Make ntuples

cd NtupleMaker/

#root -l -b -q 'GetPhotonEvents.C("SinglePhoton222_merged_processed","gmc","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_mc16a/","MC","SinglePhoton222_NoSys")'

#root -l -b -q 'GetBaseLineEvents.C("Zjets_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","Zjets_NoSys")'
#root -l -b -q 'GetBaseLineEvents.C("ttbar_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","ttbar_NoSys")'
#root -l -b -q 'GetBaseLineEvents.C("diboson_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","diboson_NoSys")'

#root -l -b -q 'GetPhotonEvents.C("data15-16_merged_processed","gdata","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_Data/","Data","data15-16")'
#root -l -b -q 'GetBaseLineEvents.C("data15-16_merged_processed","zdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Data/","Data","data")'

root -l -b -q 'MakeNtuple.C("SinglePhoton222_merged_processed","gmc","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_mc16a/","MC","SinglePhoton222_NoSys","photon")'

root -l -b -q 'MakeNtuple.C("Zjets_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","Zjets_NoSys","non-photon")'
root -l -b -q 'MakeNtuple.C("ttbar_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","ttbar_NoSys","non-photon")'
root -l -b -q 'MakeNtuple.C("diboson_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/","MC","diboson_NoSys","non-photon")'

root -l -b -q 'MakeNtuple.C("data15-16_merged_processed","gdata","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_Data/","Data","data15-16","photon")'
root -l -b -q 'MakeNtuple.C("data15-16_merged_processed","zdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Data/","Data","data","non-photon")'

cd ..

# Smear photons

cd PhotonSmearing/

#root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222_merged_processed","mm","MC","data15-16",4)'
root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222_merged_processed","ee","MC","data15-16",0)'
root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222_merged_processed","mm","MC","data15-16",0)'
#root -l -b -q 'GetPhotonSmearing.C("data15-16_merged_processed","mm","Data","data15-16",5)'
root -l -b -q 'GetPhotonSmearing.C("data15-16_merged_processed","ee","Data","data15-16",0)'
root -l -b -q 'GetPhotonSmearing.C("data15-16_merged_processed","mm","Data","data15-16",0)'

cd ..

# Copy smearing folder to reweighting folder

# Do reweighting

cd Reweighting/

root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","Data","NoSmear",1)'
root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","Data","NoSmear",2)'

root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","Data","NoSmear",1)'
root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","Data","NoSmear",2)'

#root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","Data","DataSmear",1)' 
#root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","Data","DataSmear",2)' 

root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","MC","NoSmear",1)'
root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","MC","NoSmear",2)'

#root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","MC","McSmear",1)'
#root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","MC","McSmear",2)'

root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","MC","NoSmear",1)'
root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","MC","NoSmear",2)'

cd ..

# Make analysis plots

cd Plotting/

root -l -b -q 'quickDraw_MC.C("data15-16","ee","MET","NoSmear")'
root -l -b -q 'quickDraw_MC.C("data15-16","mm","MET","NoSmear")'
root -l -b -q 'quickDraw_Data.C("data15-16","ee","MET","NoSmear")'
root -l -b -q 'quickDraw_Data.C("data15-16","mm","MET","NoSmear")'
root -l -b -q 'quickDraw_photonPredictionsSR.C("data15-16","ee","MET","NoSmear","gdata")'
root -l -b -q 'quickDraw_photonPredictionsSR.C("data15-16","mm","MET","NoSmear","gdata")'

cd ..
