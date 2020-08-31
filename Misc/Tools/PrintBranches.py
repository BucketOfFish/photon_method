import ROOT
import sys


my_file_name = sys.argv[1]

# my_file_name = "/public/data/Photon/NewSamples/ReducedNtuples/data17_data_bkg.root"
# my_file_name = "/public/data/Photon/NewSamples/ReducedNtuples/data17_data_photon.root"
# my_file_name = "/public/data/Photon/NewSamples/ReducedNtuples/mc16a_diboson.root"
# my_file_name = "/public/data/Photon/NewSamples/ReducedNtuples/mc16a_SinglePhoton222.root"

# my_file_name = "/public/data/Photon/OldSamples/Ntuples/bkg_data/data17_bkg.root"
# my_file_name = "/public/data/Photon/OldSamples/Ntuples/g_data/data17_photon.root"
# my_file_name = "/public/data/Photon/OldSamples/Ntuples/bkg_mc/mc16a_diboson.root"
# my_file_name = "/public/data/Photon/OldSamples/Ntuples/g_mc/mc16a_SinglePhoton222.root"

my_file = ROOT.TFile(my_file_name)
my_tree = my_file.BaselineTree

for branch in my_tree.GetListOfBranches():
    print(branch.GetName())
