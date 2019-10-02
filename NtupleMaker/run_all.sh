# Example usages:
# source run_all.sh all
# source run_all.sh Data

TYPES=("Data" "MC")

if [ $# == 1 ]; then
    if [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
else
    echo "Unrecognized arguments"
    return
fi

for PHOTON in "${TYPES[@]}"
do
    if [ $PHOTON == "MC" ]; then
        for PERIOD in "mc16a" "mc16cd" "mc16e"
        do
            root -l -b -q 'MakeNtuple.C("g_mc","'$PERIOD'","'$PHOTON_MC_PATH$PERIOD'/","SinglePhoton222","photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","Zjets","non-photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","ttbar","non-photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH$PERIOD'/","diboson","non-photon")'
        done
    else
        for PERIOD in "data15-16" "data17" "data18"
        do
            root -l -b -q 'MakeNtuple.C("g_data","'$PERIOD'","'$PHOTON_DATA_PATH'","data","photon")'
            root -l -b -q 'MakeNtuple.C("bkg_data","'$PERIOD'","'$BKG_DATA_PATH'","data","non-photon")'
        done
    fi
done
