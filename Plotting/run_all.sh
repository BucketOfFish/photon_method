# Example usages:
# source run_all.sh
# source run_all.sh data17
# source run_all.sh MC
# source run_all.sh data15-16 Data

PERIODS=("data15-16" "data17" "data18")
TYPES=("MC" "Data")
REGIONS="SR SRC SRCZ SRLow4 SRLowZ SRMed4 SRMedZ SRHigh4 SRHighZ VRC VRCZ VRLow4 VRLowZ VRMed4 VRMedZ VRHigh4 VRHighZ"
FEATURES="METl METt met_Et"
GETPHOTONYIELDONLY="false"
BLIND="true"

if [ $# == 1 ]; then
    if [[ ${PERIODS[*]} =~ $(echo $1) ]]; then
        PERIODS=( $1 )
    elif [[ ${TYPES[*]} =~ $(echo $1) ]]; then
        TYPES=( $1 )
    else
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
fi

for PERIOD in "${PERIODS[@]}"
do
    for TYPE in "${TYPES[@]}"
    do
        root -l -b -q 'quickDraw.C("'"$PERIOD"'","'"$TYPE"'","'"$FEATURES"'","'"$REGIONS"'",'"$BLIND"','"$GETPHOTONYIELDONLY"')' | tee 'Output/'$PERIOD'_'$TYPE'_output.txt'
    done
done
