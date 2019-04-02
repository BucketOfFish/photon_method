root -l -b -q 'GetPhotonEvents.C("SinglePhoton222_merged_processed","gmc","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_mc16a/",0,"SinglePhoton222_NoSys")'

root -l -b -q 'GetBaseLineEvents.C("Zjets_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/",0,"Zjets_NoSys")'
root -l -b -q 'GetBaseLineEvents.C("ttbar_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/",0,"ttbar_NoSys")'
root -l -b -q 'GetBaseLineEvents.C("diboson_merged_processed","ZMC16a","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_mc16a/",0,"diboson_NoSys")'

root -l -b -q 'GetPhotonEvents.C("data15-16_merged_processed","gdata","/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_Data/",1,"data15-16")'
root -l -b -q 'GetBaseLineEvents.C("data15-16_merged_processed","zdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Data/",1,"data")'
