#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../../Main/Settings.cpp"

string oldpath = "/public/data/SUSY_Systematics/Unskimmed/Zjets/";
string newpath = "/public/data/SUSY_Systematics/Skimmed/EWKPreselection/";

//string selection = cuts::selections["bkg_baseline"].GetTitle();
//string selection = cuts::selections["strong_preselection"].GetTitle();
//string selection = cuts::strong_preselection_noDPhi.GetTitle();
TCut selection_cut = cuts::lep2 + cuts::is_SF + cuts::is_OS + cuts::jet2 + cuts::mll_12 + "met_Sign > 6" +
    "nBJet20_MV2c10_FixedCutBEff_77 == 0" + "mll < 150" + "met_Et > 100";
string selection = selection_cut.GetTitle();

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
        "SUSY2_Bkgs_mc16a/diboson_merged_processed.root", "SUSY2_Bkgs_mc16cd/diboson_merged_processed.root",
        "SUSY2_Bkgs_mc16e/diboson_merged_processed.root",
        }},
    {"ttbar",
        {
        "SUSY2_Bkgs_mc16a/ttbar_merged_processed.root", "SUSY2_Bkgs_mc16cd/ttbar_merged_processed.root",
        "SUSY2_Bkgs_mc16e/ttbar_merged_processed.root",
        }},
};
string treename = "Zjets_NoSys";
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
//vector<string> branches_to_copy = vector<string> { // ttbar
    //"Ht30", "LHE3Weight_muR0p5,muF0p5", "LHE3Weight_muR0p5,muF1", "LHE3Weight_muR1,muF0p5",
    //"LHE3Weight_muR0p5,muF2", "LHE3Weight_muR2,muF0p5", "LHE3Weight_muR1,muF2", "LHE3Weight_muR2,muF1",
    //"LHE3Weight_muR2,muF2", "LHE3Weight_nominal", "LHE3Weight_PDFset265000", "LHE3Weight_PDFset266000",
    //"Ptll", "RandomRunNumber", "bTagWeight", "eventWeight", "genWeight",
    //"globalDiLepTrigSF", "jvtWeight", "lepCharge", "lepFlavor", "lepPt", "leptonWeight", "met_Et", "met_Sign",
    //"minDPhi2JetsMet", "mll", "mt2leplsp_0", "nJet30", "nLep_base", "nLep_signal", "pileupWeight",
    //"trigMatch_2LTrigOR", "nBJet20_MV2c10_FixedCutBEff_77", "mjj", "jetPt", "Rll", "dPhiMetJet1", "dPhiPllMet",
//};

void skim_ntuples_with_selection() {
    TreeCreator *reducer = new TreeCreator();

    for (auto filename : filenames) {
        reducer->read(oldpath + filename, treename);

        reducer->setBranchesToCopy(branches_to_copy);
        reducer->setCut(selection);

        reducer->write(newpath + filename, treename);
    }
}
