# Setup environment
source setup_env.sh

# Make ntuples

cd NtupleMaker/
source run_all.sh all
cd ..

# Smear photons

cd PhotonSmearing/
source run_all.sh all
cd ..

# Copy smearing folder to reweighting folder

# Do reweighting

cd Reweighting/
source run_all.sh all
cd ..

# Make analysis plots

cd Plotting/
source run_all.sh all
cd ..
