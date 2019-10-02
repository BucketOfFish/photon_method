for CHANNEL in "ee" "mm"
do
    for PERIOD in "mc16a" "mc16cd" "mc16e"
    do
        root -l -b -q 'GetPhotonSmearing.C("SinglePhoton222","'$PERIOD'","'$CHANNEL'",0)'
    done

    for PERIOD in "data15-16" "data17" "data18"
    do
        root -l -b -q 'GetPhotonSmearing.C("photon","'$PERIOD'","'$CHANNEL'",0)'
    done
done
