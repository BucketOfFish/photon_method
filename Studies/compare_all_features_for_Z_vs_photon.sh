# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh mc16a
# source run_all.sh MT2
# source run_all.sh mm mc16cd
# source run_all.sh mm mc16cd MT2

CHANNELS=("ee" "mm")
PERIODS=("mc16a" "mc16cd")
FEATURES=("lepPt" "lepEta" "lepPhi" "MT2" "mll")

# non-uniform, uniform, Drell-Yan, sin3
SAMPLING="sin3"

if [ $# == 1 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]]; then
        CHANNELS=( $1 )
    elif [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${FEATURES[*]} =~ $(echo $1) ]]; then
        FEATURES=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 2 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
    else
        echo "Unrecognized arguments"
        return
    fi
elif [ $# == 3 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${FEATURES[*]} =~ $(echo $3) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        FEATURES=( $3 )
    else
        echo "Unrecognized arguments"
        return
    fi
else
    echo "Unrecognized arguments"
    return
fi

for CHANNEL in "${CHANNELS[@]}"
do
    for PERIOD in "${PERIODS[@]}"
    do
        for FEATURE in "${FEATURES[@]}"
        do
            root -l -b -q 'compare_feature_for_Z_vs_photon.C("'$PERIOD'","'$CHANNEL'","baseline","'$SAMPLING'","'$FEATURE'")'
        done
    done
done
