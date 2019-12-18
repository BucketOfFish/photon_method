# Example usages:
# source run_all.sh all
# source run_all.sh ee
# source run_all.sh data17
# source run_all.sh MC
# source run_all.sh SR
# source run_all.sh "Ptll>200"
# source run_all.sh nJet30
# source run_all.sh mm data18 MC
# source run_all.sh mm data18 MC SR "HT>200"
# source run_all.sh mm data18 MC VR "1" METt

CHANNELS=("ee" "mm")
#PERIODS=("data15-16" "data17" "data18")
#TYPES=("Data" "MC")
#REGIONS=("VR" "SR" "VRcom" "SRZ2016" "SRlow2016" "SRmed2016" "SRhigh2016")
#CUTS=("1" "mll>81&&mll<101" "Ptll>200" "Ptll>400" "HT>200" "HT>400")
#FEATURES=("met_Et" "METt" "METl" "nJet30" "HT" "lepPt[0]" "lepPt[1]" "lep_eta[0]" "lep_eta[1]" "dPhiMetJet1" "dPhiMetJet2" "dPhiMetJet12Min" "mll" "mt2leplsp_0")
PERIODS=("data18")
TYPES=("MC")
REGIONS=("SR")
#CUTS=("1" "Ptll>50&&Ptll<75" "Ptll>100&&Ptll<125" "Ptll>300&&Ptll<400" "Ptll>1000&&Ptll<1200")
#FEATURES=("met_Et" "METt" "METl" "mll")
#CUTS=("1" "mll<50" "METt>100||METt<-100" "METl>100" "met_Et>150")
CUTS=("met_Et>150")
FEATURES=("Ptll")

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
    elif [[ ${REGIONS[*]} =~ $(echo $1) ]]; then
        REGIONS=( $1 )
    elif [ $1 != "all" ]; then
        echo "Unrecognized argument"
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
elif [ $# == 5 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${TYPES[*]} =~ $(echo $3) ]] && [[ ${REGIONS[*]} =~ $(echo $4) ]] && [[ ${CUTS[*]} =~ $(echo $5) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        TYPES=( $3 )
        REGIONS=( $4 )
        CUTS=( $5 )
    else
        echo "Unrecognized arguments"
        return
    fi
elif [ $# == 6 ]; then
    if [[ ${CHANNELS[*]} =~ $(echo $1) ]] && [[ ${PERIODS[*]} =~ $(echo $2) ]] && [[ ${TYPES[*]} =~ $(echo $3) ]] && [[ ${REGIONS[*]} =~ $(echo $4) ]] && [[ ${CUTS[*]} =~ $(echo $5) ]] && [[ ${FEATURES[*]} =~ $(echo $6) ]]; then
        CHANNELS=( $1 )
        PERIODS=( $2 )
        TYPES=( $3 )
        REGIONS=( $4 )
        CUTS=( $5 )
        FEATURES=( $6 )
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
            for FEATURE in "${FEATURES[@]}"
            do
                for CUT in "${CUTS[@]}"
                do
                    for REGION in "${REGIONS[@]}"
                    do
                        #sem -j 6 root -l -b -q \''quickDraw.C("'$PERIOD'","'$CHANNEL'","'$FEATURE'","'$PHOTON'","'$REGION'","'$CUT'")'\'
                        root -l -b -q 'quickDraw.C("'$PERIOD'","'$CHANNEL'","'$FEATURE'","'$PHOTON'","'$REGION'","'$CUT'")'
                    done
                done
            done
        done
    done
done
