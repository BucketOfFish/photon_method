for PERIOD in "mc16a" "mc16cd" "mc16e"
do
    root -l -b -q 'MakeNtuple.C("g_mc","'$PERIOD'","'$PHOTON_MC_PATH$PERIOD'/","SinglePhoton222","photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","Zjets","non-photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","ttbar","non-photon")'
    root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","diboson","non-photon")'
done

for PERIOD in "data15-16" "data17" "data18"
do
    root -l -b -q 'MakeNtuple.C("g_data","'$PERIOD'","'$PHOTON_DATA_PATH'","data","photon")'
    root -l -b -q 'MakeNtuple.C("bkg_data","'$PERIOD'","'$BKG_DATA_PATH'","data","non-photon")'
done
