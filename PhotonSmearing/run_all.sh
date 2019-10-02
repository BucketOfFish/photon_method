# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh MC
# source run_all.sh mm Data

CHANNELS=("ee" "mm")
TYPES=("Data" "MC")

if [ $# == 1 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]]; then
        CHANNELS=( $1 )
    elif [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
        return
    fi
elif [ $# == 2 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${TYPES[*]} =~ $(echo $2) ]]; then
        CHANNELS=( $1 )
        TYPES=( $2 )
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
    for PHOTON in "${TYPES[@]}"
    do
        if [ $PHOTON == "MC" ]; then
            for PERIOD in "mc16a" "mc16cd" "mc16e"
            do
                root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222","'$PERIOD'","'$CHANNEL'",0)'
            done
        else
            for PERIOD in "data15-16" "data17" "data18"
            do
                root -l -b -q 'GetPhotonSmearing.C("photon","'$PERIOD'","'$CHANNEL'",0)'
            done
        fi
    done
done
