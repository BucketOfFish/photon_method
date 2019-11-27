# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh MC
# source run_all.sh data17
# source run_all.sh ee data17
# source run_all.sh data17 MC
# source run_all.sh mm data18 Data

CHANNELS=("ee" "mm")
TYPES=("Data" "MC")
PERIODS=("data15-16" "data17" "data18")

if [ $# == 1 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]]; then
        CHANNELS=( $1 )
    elif [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 2 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
    elif [[ ${PERIODS[*]} =~ $(echo $1) ]] && [[ ${TYPES[*]} =~ $(echo $2) ]]; then
        PERIODS=( $1 )
        TYPES=( $2 )
    else
        echo "Unrecognized arguments"
        return
    fi
elif [ $# == 3 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${TYPES[*]} =~ $(echo $3) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        TYPES=( $3 )
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
    for TYPE in "${TYPES[@]}"
    do
        for PERIOD in "${PERIODS[@]}"
        do
            root -l -b -q 'GetPhotonSmearing.C("'$PERIOD'","'$CHANNEL'","'$TYPE'")'
        done
    done
done
