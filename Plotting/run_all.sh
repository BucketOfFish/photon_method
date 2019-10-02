for CUT in "1" "mll>81&&mll<101" "Ptll>200" "Ptll>400" "HT>200" "HT>400"
do
    for CHANNEL in "ee" "mm"
    do
        for FEATURE in "met_Et" "METt" "METl" "nJet30" "bjet_n" "HT"
        do
            for PHOTON in "Data" "MC"
            do
                for PERIOD in "data15-16" "data17" "data18"
                do
                    sem -j 6 root -l -b -q \''quickDraw.C("'$PERIOD'","'$CHANNEL'","'$FEATURE'","NoSmear","'$PHOTON'","'$CUT'")'\'
                done
            done
        done
    done
done

