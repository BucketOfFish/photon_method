for CHANNEL in "ee" "mm"
do
    for PERIOD in "data15-16" "data17" "data18"
    do
        for PHOTON in "Data" "MC"
        do
            root -l -b -q 'GetPhotonReweighting.C("'$PERIOD'","'$CHANNEL'","'$PHOTON'","NoSmear","Ptll")'
        done
    done
done

