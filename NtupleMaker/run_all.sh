# Example usages:
# source run_all.sh all
# source run_all.sh Data
# source run_all.sh data17
# source run_all.sh data18 MC

PERIODS=("data15-16" "data17" "data18")
TYPES=("Data" "MC")
EVERYNENTRIES=1

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
            root -l -b -q 'MakeNtuple.C("'$PHOTON_MC_PATH'","g_mc","'$PERIOD'","SinglePhoton222","photon",'$EVERYNENTRIES')'
            root -l -b -q 'MakeNtuple.C("'$BKG_MC_PATH'","bkg_mc","'$PERIOD'","Zjets","non-photon",'$EVERYNENTRIES')'
            root -l -b -q 'MakeNtuple.C("'$BKG_MC_PATH'","bkg_mc","'$PERIOD'","ttbar","non-photon",'$EVERYNENTRIES')'
            root -l -b -q 'MakeNtuple.C("'$BKG_MC_PATH'","bkg_mc","'$PERIOD'","diboson","non-photon",'$EVERYNENTRIES')'
        else
            root -l -b -q 'MakeNtuple.C("'$PHOTON_DATA_PATH'","g_data","'$PERIOD'","data","photon",'$EVERYNENTRIES')'
            root -l -b -q 'MakeNtuple.C("'$BKG_DATA_PATH'","bkg_data","'$PERIOD'","data","non-photon",'$EVERYNENTRIES')'
        fi
    done
done
