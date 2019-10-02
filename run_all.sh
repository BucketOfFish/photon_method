# Setup environment
source setup_env.sh

# Make ntuples

cd NtupleMaker/
source run_all.sh
cd ..

# Smear photons

cd PhotonSmearing/
source run_all.sh
cd ..

# Copy smearing folder to reweighting folder

# Do reweighting

cd Reweighting/
source run_all.sh
cd ..

# Make analysis plots

cd Plotting/
source run_all.sh
cd ..
