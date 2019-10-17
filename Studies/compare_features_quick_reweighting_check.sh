PERIOD="mc16e"
#CHANNELS=("ee" "mm")
CHANNELS=("ee")
REGIONS=("baseline" "VRcom")
#FEATURES=("Ptll" "met_Et" "HT")
FEATURES=("Ptll" "met_Et")

for CHANNEL in "${CHANNELS[@]}"
do
    for REGION in "${REGIONS[@]}"
    do
        for FEATURE in "${FEATURES[@]}"
        do
            root -l -b -q 'compare_feature_for_Z_vs_photon.C("'$PERIOD'","'$CHANNEL'","'$REGION'","histogram","'$FEATURE'")'
        done
    done
done
