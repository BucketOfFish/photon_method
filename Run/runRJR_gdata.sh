cd ../Root 
#lsetup ROOT

echo "Doing data15-16..."
rm -f logfiles/runRJR_gdata_Data15-16.log
root -l -b -q 'GetPhotonEvents.C+("JETM4_Data15-16","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/JETM4/MAY16_manualDownload/merged/",1)' > logfiles/runRJR_gdata_Data15-16.log 2>&1 &
sleep 20

echo "Doing data17..."
rm -f logfiles/runRJR_gdata_Data17.log
root -l -b -q 'GetPhotonEvents.C+("JETM4_Data17","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/JETM4/MAY16_manualDownload/merged/",1)' > logfiles/runRJR_gdata_Data17.log 2>&1 & 
sleep 20

