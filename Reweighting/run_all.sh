# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh data17
# source run_all.sh MC
# source run_all.sh data18 MC
# source run_all.sh mm data18 MC

CHANNELS=("ee" "mm")
PERIODS=("data15-16" "data17" "data18")
TYPES=("Data" "MC")
SMEARING="McSmear"

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
    if [[ ${PERIODS[*]} =~ $(echo $1) ]] && [[ ${TYPES[*]} =~ $(echo $2) ]]; then
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
    for PERIOD in "${PERIODS[@]}"
    do
        for TYPE in "${TYPES[@]}"
        do
            root -l -b -q 'GetPhotonReweighting.C("'$PERIOD'","'$CHANNEL'","'$TYPE'","'$SMEARING'","Ptll")'
        done
    done
done
