# Example usages:
# source run_all.sh all
# source run_all.sh Data
# source run_all.sh data17
# source run_all.sh data18 MC

PERIODS=("data15-16" "data17" "data18")
TYPES=("Data" "MC")

if [ $# == 1 ]; then
    if [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 2 ]; then
    if [[ ${PERIODS[*]} =~ $(echo $1) ]] && [[ ${TYPES[*]} =~ $(echo $2) ]]; then
        PERIODS=( $1 )
        TYPES=( $2 )
    else
        echo "Unrecognized arguments"
        return
    fi
else
    echo "Unrecognized arguments"
    return
fi

for TYPE in "${TYPES[@]}"
do
    for PERIOD in "${PERIODS[@]}"
    do
        if [ $TYPE == "MC" ]; then
            root -l -b -q 'MakeNtuple.C("g_mc","'$PERIOD'","'$PHOTON_MC_PATH'","SinglePhoton222","photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH'","Zjets","non-photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH'","ttbar","non-photon")'
            root -l -b -q 'MakeNtuple.C("bkg_mc","'$PERIOD'","'$BKG_MC_PATH'","diboson","non-photon")'
        else
            root -l -b -q 'MakeNtuple.C("g_data","'$PERIOD'","'$PHOTON_DATA_PATH'","data","photon")'
            root -l -b -q 'MakeNtuple.C("bkg_data","'$PERIOD'","'$BKG_DATA_PATH'","data","non-photon")'
        fi
    done
done
