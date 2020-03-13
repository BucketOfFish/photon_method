# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh mc16a
# source run_all.sh MT2
# source run_all.sh VRcom
# source run_all.sh mm mc16cd VRcom
# source run_all.sh mm mc16cd baseline MT2

CHANNELS=("ee" "mm")
PERIODS=("mc16a" "mc16cd" "mc16e")
FEATURES=("met_Et" "METl" "METl" "nJet30" "bjet_n" "HT" "lepPt" "lepEta" "lepPhi" "dPhiMetJet1" "dPhiMetJet2" "MT2" "mll" "Ptll")
REGIONS=("inclusive" "baseline" "VRcom")

# non-uniform, uniform, Drell-Yan, sin3, histogram
SAMPLING="histogram"

if [ $# == 1 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]]; then
        CHANNELS=( $1 )
    elif [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${FEATURES[*]} =~ $(echo $1) ]]; then
        FEATURES=( $1 )
    elif [[ ${REGIONS[*]} =~ $(echo $1) ]]; then
        REGIONS=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 3 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${REGIONS[*]} =~ $(echo $3) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        REGIONS=( $2 )
    else
        echo "Unrecognized arguments"
        return
    fi
elif [ $# == 4 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${REGIONS[*]} =~ $(echo $3) ]] && [[ ${FEATURES[*]} =~ $(echo $4) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        REGIONS=( $3 )
        FEATURES=( $4 )
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
            for REGION in "${REGIONS[@]}"
            do
                root -l -b -q 'compare_feature_for_Z_vs_photon.C("'$PERIOD'","'$CHANNEL'","'$REGION'","'$SAMPLING'","'$FEATURE'")'
            done
        done
    done
done
