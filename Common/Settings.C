#ifndef COMMON_SETTINGS
#define COMMON_SETTINGS

#include "CommonLibraries.C"
#include "CommonCuts.C"
#include "CommonBins.C"
#include "CommonFunctions.C"

std::string sample_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/TinySamples/";
std::string sampling_method = "HistogramSampling";

std::string ntuple_path =  sample_folder + "Ntuples/";
std::string smearing_path = sample_folder + "/" + sampling_method + "/SmearedNtuples/";
std::string reweighting_path = sample_folder + "/" + sampling_method + "/ReweightedNtuples/";
std::string plots_path = sample_folder + "/" + sampling_method + "/Plots/";

#endif
