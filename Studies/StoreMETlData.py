from __future__ import print_function
import ROOT
import pickle as pkl
import numpy as np
import sys


sample_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/Samples/"
ntuple_path = sample_folder + "Ntuples/"
mc_period = "mc16e"
pt_bins = [40, 45, 50, 55, 60, 70, 80, 100, 120, 140, 160, 180, 200, 220, 260, 280, 300, 350, 400, 600, 1000]
max_n_events = -1

# ----------
# Get Z METl
# ----------

Z_mc_file = ROOT.TFile(ntuple_path + "/bkg_mc/" + mc_period + "_Zjets.root")
Z_mc_tree = Z_mc_file.BaselineTree

Z_mc_tree.SetBranchStatus("*", 0)
used_branches = ['Ptll', 'nJet30', 'jet_pT', 'nLep_signal', 'lepPt']
for branch in used_branches:
    Z_mc_tree.SetBranchStatus(branch, 1)

def process_Z_event(i, event):
    if i % 100000 == 0:
        print(str(i) + " events processed", end="\r")
        sys.stdout.flush()
    if event.Ptll > 40 and event.nJet30 >= 2 and event.jet_pT[0] > 30 and event.jet_pT[1] > 30 and event.nLep_signal >= 2 and event.lepPt[0] > 25 and event.lepPt[1] > 25:
        bin_n = np.digitize([event.Ptll], pt_bins)[0] - 1
        return(bin_n, event.METl)
    else:
        return(-1, -1)

print("Processing " + str(Z_mc_tree.GetEntries()) + " Z events")
Z_METl = [(-1, -1) for _ in range(Z_mc_tree.GetEntries())]
for i, event in enumerate(Z_mc_tree):
    if max_n_events > -1 and i > max_n_events:
        break
    Z_METl[i] = process_Z_event(i, event)
print()
Z_mc_file.Close()

print("Binning Z events by pT")
Z_METl_binned = [[] for _ in range(len(pt_bins)-1)]
for i in range(len(Z_METl_binned)):
    Z_METl_binned[i] = [j[1] for j in Z_METl if j[0]==i]
del Z_METl

# ---------------
# Get Photon METl
# ---------------

photon_mc_file = ROOT.TFile(ntuple_path + "/g_mc/" + mc_period + "_SinglePhoton222.root")
photon_mc_tree = photon_mc_file.BaselineTree

photon_mc_tree.SetBranchStatus("*", 0)
used_branches = ['Ptll', 'nJet30', 'jet_pT']
for branch in used_branches:
    photon_mc_tree.SetBranchStatus(branch, 1)

def process_photon_event(i, event):
    if max_n_events > -1 and i > max_n_events:
        return(-1, -1)
    if i % 100000 == 0:
        print(str(i) + " events processed", end="\r")
        sys.stdout.flush()
    if event.Ptll > 40 and event.nJet30 >= 2 and event.jet_pT[0] > 30 and event.jet_pT[1] > 30:
        bin_n = np.digitize([event.Ptll], pt_bins)[0]-1
        return(bin_n, event.METl)
    else:
        return(-1, -1)

print("Processing " + str(photon_mc_tree.GetEntries()) + " photon events")
photon_METl = [(-1, -1) for _ in range(photon_mc_tree.GetEntries())]
for i, event in enumerate(photon_mc_tree):
    photon_METl[i] = process_photon_event(i, event)
print()
photon_mc_file.Close()

print("Binning photon events by pT")
photon_METl_binned = [[] for _ in range(len(pt_bins)-1)]
for i in range(len(photon_METl_binned)):
    photon_METl_binned[i] = [j[1] for j in photon_METl if j[0]==i]
del photon_METl

# ---------
# Save File
# ---------

with open("binned_METl.pkl", "w") as output_file:
    data = {}
    data["Z_METl"] = Z_METl_binned
    data["photon_METl"] = photon_METl_binned
    pkl.dump(data, output_file)
