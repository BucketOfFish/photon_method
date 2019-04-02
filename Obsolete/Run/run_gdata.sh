#lsetup ROOT

echo "Doing data15-16..."
rm -f logfiles/run_gdata_data15-16.log
root -l -b -q '../Root/GetPhotonEvents.C+("data15-16","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.2/JETM4/JETM4_Data/JETM4_v1.2_v5/merged_manual/",1)' > logfiles/run_gdata_data15-16.log 2>&1 &
sleep 20

echo "Doing data17..."
rm -f logfiles/run_gdata_data17.log
root -l -b -q '../Root/GetPhotonEvents.C+("data17","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.2/JETM4/JETM4_Data/JETM4_v1.2_v5/merged_manual/",1)' > logfiles/run_gdata_data17.log 2>&1 &
sleep 20


#echo "Doing data18..."
rm -f logfiles/run_gdata_data18.log
root -l -b -q '../Root/GetPhotonEvents.C+("data18","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.2/JETM4/JETM4_Data/JETM4_v1.2_v5/merged_manual/",1)' > logfiles/run_gdata_data18.log 2>&1 &
#sleep 20
