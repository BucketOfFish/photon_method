import ROOT
import pickle as pkl
import numpy as np


sample_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/Samples/"
ntuple_path = sample_folder + "Ntuples/"
mc_period = "mc16e"
pt_bins = [40, 45, 50, 55, 60, 70, 80, 100, 120, 140, 160, 180, 200, 220, 260, 280, 300, 350, 400, 600, 1000]
# max_n_events = 10000000
max_n_events = -1

# ----------
# Get Z METl
# ----------

Z_mc_file = ROOT.TFile(ntuple_path + "/bkg_mc/" + mc_period + "_Zjets.root")
Z_mc_tree = Z_mc_file.BaselineTree

Z_METl = [[] for _ in range(len(pt_bins))]
for i, event in enumerate(Z_mc_tree):
    if max_n_events > -1 and i > max_n_events:
        break
    if i % 100000 == 0:
        print(str(i) + " events processed")
    if event.Ptll > 40 and event.nJet30 >= 2 and event.jet_pT[0] > 30 and event.jet_pT[1] > 30 and event.nLep_signal >= 2 and event.lepPt[0] > 25 and event.lepPt[1] > 25:
        bin_n = np.digitize([event.Ptll], pt_bins)[0] - 1
        Z_METl[bin_n].append(event.METl)

Z_mc_file.Close()

# ---------------
# Get Photon METl
# ---------------

photon_mc_file = ROOT.TFile(ntuple_path + "/g_mc/" + mc_period + "_SinglePhoton222.root")
photon_mc_tree = photon_mc_file.BaselineTree

photon_METl = [[] for _ in range(len(pt_bins)-1)]
for i, event in enumerate(photon_mc_tree):
    if max_n_events > -1 and i > max_n_events:
        break
    if event.Ptll > 40 and event.nJet30 >= 2 and event.jet_pT[0] > 30 and event.jet_pT[1] > 30:
        bin_n = np.digitize([event.Ptll], pt_bins)[0]-1
        photon_METl[bin_n].append(event.METl)

photon_mc_file.Close()

# ---------
# Save File
# ---------

with open("binned_METl.pkl", "w") as output_file:
    data = {}
    data["Z_METl"] = Z_METl
    data["photon_METl"] = photon_METl
    pkl.dump(data, output_file)
