#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../../Main/Settings.cpp"

string oldpath = "/public/data/SUSY_Systematics/Unskimmed/Diboson/";
string newpath = "/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/";
string treename = "diboson_NoSys";

//string selection = cuts::selections["bkg_baseline"].GetTitle();
//string selection = cuts::selections["strong_preselection"].GetTitle();
string selection = cuts::strong_preselection_noDPhi.GetTitle();

map<string, vector<string>> filename_sets {
    {"reweighted",
        {"data15-16_data_photon_ee.root", "data18_data_photon_ee.root", "mc16cd_SinglePhoton222_ee.root",
        "data15-16_data_photon_mm.root", "data18_data_photon_mm.root", "mc16cd_SinglePhoton222_mm.root",
        "data17_data_photon_ee.root", "mc16a_SinglePhoton222_ee.root", "mc16e_SinglePhoton222_ee.root",
        "data17_data_photon_mm.root", "mc16a_SinglePhoton222_mm.root", "mc16e_SinglePhoton222_mm.root",}},
    {"reweighted_noMC",
        {"data15-16_data_photon_ee.root", "data18_data_photon_ee.root",
        "data15-16_data_photon_mm.root", "data18_data_photon_mm.root",
        "data17_data_photon_ee.root", "data17_data_photon_mm.root",}},
    {"reduced",
        {
        "data15-16_data_bkg.root", "mc16a_singleTop.root", "mc16e_lowMassDY.root",
        "mc16a_topOther.root", "mc16cd_singleTop.root",
        "data17_data_bkg.root", "mc16a_triboson.root", "mc16cd_topOther.root", "mc16e_singleTop.root",
        "mc16a_ttbar.root", "mc16cd_triboson.root", "mc16e_topOther.root",
        "data18_data_bkg.root", "mc16cd_ttbar.root", "mc16e_triboson.root",
        "mc16a_Wjets.root", "mc16e_ttbar.root",
        "mc16a_diboson.root", "mc16a_Zjets.root", "mc16cd_Wjets.root",
        "mc16a_higgs.root", "mc16cd_diboson.root", "mc16cd_Zjets.root", "mc16e_Wjets.root",
        "mc16a_lowMassDY.root", "mc16cd_higgs.root", "mc16e_diboson.root", "mc16e_Zjets.root",
        "mc16cd_lowMassDY.root", "mc16e_higgs.root",
        }},
    {"reduced_diboson",
        {
        "mc16a_diboson.root", "mc16cd_diboson.root", "mc16e_diboson.root",
        }},
    {"Zjets",
        {
        "mc16a_Zjets_merged_processed.root", "mc16cd_Zjets_merged_processed.root", "mc16e_Zjets_merged_processed.root",
        }},
    {"diboson",
        {
        "diboson_merged_processed.root",
        }},
    {"ttbar",
        {
        "ttbarDilep_410472_merged_processed.root", "ttbar_merged_processed.root",
        }},
};
vector<string> filenames = filename_sets["Zjets"];

vector<string> branches_to_copy = vector<string> {
    "Ht30", "LHE3Weight_MUR0.5_MUF0.5_PDF261000", "LHE3Weight_MUR0.5_MUF1_PDF261000",
    "LHE3Weight_MUR1_MUF0.5_PDF261000", "LHE3Weight_MUR1_MUF1_PDF13000", "LHE3Weight_MUR1_MUF1_PDF25300",
    "LHE3Weight_MUR1_MUF1_PDF261000", "LHE3Weight_MUR1_MUF2_PDF261000", "LHE3Weight_MUR2_MUF1_PDF261000",
    "LHE3Weight_MUR2_MUF2_PDF261000", "Ptll", "RandomRunNumber", "bTagWeight", "eventWeight", "genWeight",
    "globalDiLepTrigSF", "jvtWeight", "lepCharge", "lepFlavor", "lepPt", "leptonWeight", "met_Et", "met_Sign",
    "minDPhi2JetsMet", "mll", "mt2leplsp_0", "nJet30", "nLep_base", "nLep_signal", "pileupWeight",
    "trigMatch_2LTrigOR", "nBJet20_MV2c10_FixedCutBEff_77", "mjj", "jetPt", "Rll", "dPhiMetJet1", "dPhiPllMet",
};

void skim_ntuples_with_selection() {
    TreeCreator *reducer = new TreeCreator();

    for (auto filename : filenames) {
        reducer->read(oldpath + filename, treename);

        reducer->setBranchesToCopy(branches_to_copy);
        reducer->setCut(selection);

        reducer->write(newpath + filename, treename);
    }
}
