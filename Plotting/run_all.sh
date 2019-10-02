# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh data17
# source run_all.sh MC
# source run_all.sh "Ptll>200"
# source run_all.sh nJet30
# source run_all.sh mm data18 MC "HT>200" METt

CHANNELS=("ee" "mm")
PERIODS=("data15-16" "data17" "data18")
TYPES=("Data" "MC")
CUTS=("1" "mll>81&&mll<101" "Ptll>200" "Ptll>400" "HT>200" "HT>400")
FEATURES=("met_Et" "METt" "METl" "nJet30" "bjet_n" "HT")

if [ $# == 1 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]]; then
        CHANNELS=( $1 )
    elif [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    elif [[ ${FEATURES[*]} =~ $(echo $1) ]]; then
        FEATURES=( $1 )
    elif [[ ${CUTS[*]} =~ $(echo $1) ]]; then
        CUTS=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 5 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${TYPES[*]} =~ $(echo $3) ]] && [[ ${CUTS[*]} =~ $(echo $4) ]] && [[ ${FEATURES[*]} =~ $(echo $5) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        TYPES=( $3 )
        CUTS=( $4 )
        FEATURES=( $5 )
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
        for PHOTON in "${TYPES[@]}"
        do
            for CUT in "${CUTS[@]}"
            do
                for FEATURE in "${FEATURES[@]}"
                do
                    sem -j 6 root -l -b -q \''quickDraw.C("'$PERIOD'","'$CHANNEL'","'$FEATURE'","NoSmear","'$PHOTON'","'$CUT'")'\'
                done
            done
        done
    done
done
