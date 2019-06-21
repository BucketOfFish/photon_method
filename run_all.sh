# Setup environment

lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"

# Make ntuples

cd NtupleMaker/

export PHOTON_MC_PATH='/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_'
export PHOTON_DATA_PATH='/eos/atlas/user/y/ycao/SUSY_dataset/JETM4_v1.6/JETM4_Data/'
export BKG_MC_PATH='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Bkgs_'
export BKG_DATA_PATH='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.6/SUSY2/SUSY2_Data/'

for PERIOD in "mc16a" "mc16cd" "mc16e"
do
    root -l -b -q 'MakeNtuple.C("g_mc","'$PERIOD'","'$PHOTON_MC_PATH$PERIOD'/","SinglePhoton222","MC","photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","Zjets","MC","non-photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","ttbar","MC","non-photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","diboson","MC","non-photon")'
done

for PERIOD in "data15-16" "data17" "data18"
do
    root -l -b -q 'MakeNtuple.C("g_data","'$PERIOD'","'$PHOTON_DATA_PATH'","photon","Data","photon","'$PERIOD'")'
    root -l -b -q 'MakeNtuple.C("bkg_data","'$PERIOD'","'$BKG_DATA_PATH'","photon","Data","non-photon","data")'
done

cd ..

# Smear photons

cd PhotonSmearing/

#root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222_merged_processed","mm","MC","data15-16",4)'
root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222","mc16a","ee",0)'
root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222","mc16a","mm",0)'
#root -l -b -q 'GetPhotonSmearing.C("data15-16_merged_processed","mm","Data","data15-16",5)'
root -l -b -q 'GetPhotonSmearing.C("photon","data15-16","ee",0)'
root -l -b -q 'GetPhotonSmearing.C("photon","data15-16","mm",0)'

cd ..

# Copy smearing folder to reweighting folder

# Do reweighting

cd Reweighting/

root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","Data","NoSmear","Ptll")'
root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","Data","NoSmear","Ptll")'
root -l -b -q 'GetPhotonReweighting.C("data15-16","ee","MC","NoSmear","Ptll")'
root -l -b -q 'GetPhotonReweighting.C("data15-16","mm","MC","NoSmear","Ptll")'

cd ..

# Make analysis plots

cd Plotting/

for CUT in "1" "mll>81&&mll<101" "Ptll>200" "Ptll>400" "HT>200" "HT>400"
do
    for CHANNEL in "ee" "mm"
    do
        for FEATURE in "met_Et" "METt" "METl" "nJet30" "bjet_n" "HT"
        do
            for PHOTON in "Data" "MC"
            do
                sem -j 6 root -l -b -q \''quickDraw.C("data15-16","'$CHANNEL'","'$FEATURE'","NoSmear","'$PHOTON'","'$CUT'")'\'
            done
        done
    done
done

cd ..
