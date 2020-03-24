#ifndef COMMON_SETTINGS
#define COMMON_SETTINGS

#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
//#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iomanip> 
#include <algorithm>
#include <math.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TSpectrum.h"
#include "TVirtualFFT.h"
#include "TGraphErrors.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TRandom.h"
#include <ROOT/RDataFrame.hxx>
//#include <ROOT/RDF/RInterface.hxx>
#include "TInterpreter.h"

//--- printing colors
#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define PRED(x) KRED x RST
#define PGRN(x) KGRN x RST
#define PYEL(x) KYEL x RST
#define PBLU(x) KBLU x RST
#define PMAG(x) KMAG x RST
#define PCYN(x) KCYN x RST
#define PWHT(x) KWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

namespace cuts {

    double leading_lep_pt_cut = 25.; // also used for smearing
    double second_lep_pt_cut = 25.; // also used for smearing

    //TCut bkg_baseline("nJet30>=2 && is_OS && lepPt[0]>25.0 && lepPt[1]>25.0 && lepIsoFCTight[0] && lepIsoFCTight[1]");
    //TCut photon_baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");
    TCut bkg_baseline("nJet30>=1 && nLep_signal==2 && nLep_base==2 && (lepCharge[0]!=lepCharge[1]) && lepPt[0]>25.0 && lepPt[1]>25.0 && lepIsoFCTight[0] && lepIsoFCTight[1] && trigMatch_2LTrigOR");
    TCut photon_baseline_ntuples("nJet30>=1 && PhotonPt>15 && nLep_base==0");
    TCut photon_baseline("nJet30>=1 && gamma_pt>15 && nLep_base==0");
    TCut baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");

    TCut mm("channel==0");
    TCut ee("channel==1");
    TCut em("channel==2");
    TCut me("channel==3");
    TCut SF("channel==0 || channel==1");
    TCut DF("channel==2 || channel==3");

    TCut bkg_weight("totalWeight");
    TCut photon_weight("totalWeight");
    TCut photon_weight_rw("totalWeight*reweight_Ptll");

    TCut reweight_region("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2");

    std::unordered_map<std::string, TCut> plot_regions = {
        {"Inclusive", "1"},
        {"SRTest", "nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2"},
        {"VRTest", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && nLep_signal==2 && mll<100"},
        {"VR", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100"},
        {"VRcom", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30"},

        // from https://arxiv.org/pdf/1611.05791.pdf
        {"SRZ2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && (mll>81 && mll<101) && dPhiMetJet12Min>0.4 && Ht30>600"},
        {"SRlow2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4"},
        {"SRmed2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && Ht30>400"},
        {"SRhigh2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && Ht30>700"},

        // from https://indico.cern.ch/event/883484/contributions/3722767/attachments/1984038/3305237/20-02-10-2l.pdf
        {"SRC", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=2 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && mll<50"},
        {"SRCZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && (mll>81 && mll<101)"},
        {"SRLow4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll<150 && !(mll>81 && mll<101))"},
        {"SRLowZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll>81 && mll<101)"},
        {"SRMed4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>500 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>100 && mll<550)"},
        {"SRMedZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>500 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>81 && mll<101)"},
        {"SRHigh4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>800 && mt2leplsp_0>75 && Ptll>40 && (mll>150 && mll<950)"},
        {"SRHighZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>800 && mt2leplsp_0>75 && Ptll>40 && (mll>81 && mll<101)"},

        {"VRC", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=2 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && mll<50"},
        {"VRCZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && (mll>81 && mll<101)"},
        {"VRLow4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll<150 && !(mll>81 && mll<101))"},
        {"VRLowZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll>81 && mll<101)"},
        {"VRMed4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>500 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>100 && mll<550)"},
        {"VRMedZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>500 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>81 && mll<101)"},
        {"VRHigh4", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=4 && Ht30>800 && mt2leplsp_0>75 && Ptll>40 && (mll>150 && mll<950)"},
        {"VRHighZ", "nLep_signal==2 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=6 && Ht30>800 && mt2leplsp_0>75 && Ptll>40 && (mll>81 && mll<101)"},
    };

    TCut CR("met_Et<100.0");

    std::unordered_map<std::string, TCut> plot_region_met_portions = {
        {"Inclusive", "1"},
        {"SRTest", "1"},
        {"VRTest", "1"},
        {"VRcom", "(met_Et>100 && met_Et<200)"},

        // from https://arxiv.org/pdf/1611.05791.pdf
        {"SRZ2016", "met_Et>225"},
        {"SRlow2016", "met_Et>200"},
        {"SRmed2016", "met_Et>200"},
        {"SRhigh2016", "met_Et>200"},

        // from https://indico.cern.ch/event/883484/contributions/3722767/attachments/1984038/3305237/20-02-10-2l.pdf
        {"SRC", "met_Et>250"},
        {"SRCZ", "met_Et>250"},
        {"SRLow4", "met_Et>250"},
        {"SRLowZ", "met_Et>250"},
        {"SRMed4", "met_Et>300"},
        {"SRMedZ", "met_Et>300"},
        {"SRHigh4", "met_Et>300"},
        {"SRHighZ", "met_Et>300"},

        {"VRC", "(met_Et>150 && met_Et<250)"},
        {"VRCZ", "(met_Et>150 && met_Et<250)"},
        {"VRLow4", "(met_Et>150 && met_Et<250)"},
        {"VRLowZ", "(met_Et>150 && met_Et<250)"},
        {"VRMed4", "(met_Et>200 && met_Et<300)"},
        {"VRMedZ", "(met_Et>200 && met_Et<300)"},
        {"VRHigh4", "(met_Et>200 && met_Et<300)"},
        {"VRHighZ", "(met_Et>200 && met_Et<300)"},
    };
}

namespace bins {

    const int n_smearing_bins = 1000;
    double smearing_low = -2000;
    double smearing_high = 2000;

    const int n_pt_bins = 23;
    //double pt_bins[] = {50,75,100,125,150,175,200,250,300,400,500,700,1000,1200,1400,1600,1e10};
    double pt_bins[] = {0,30,35,40,45,50,55,60,70,80,100,120,140,160,180,200,220,260,280,300,350,400,600,1000,1e10};
    double MET_bins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
    double dphi_bin[] = {0,0.5,1.0,1.5,2.0,2.5,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

    const int n_METl_bins = 25;
    double METl_bins[] = {-1e10,-1000,-700,-500,-400,-300,-250,-200,-150,-100,-60,-40,-20,20,40,60,100,150,200,250,300,400,500,700,1000,1e10};

    const int n_mll_bins = 43;
    double mll_bin[] = {12,20,30,40,50,60,70,80,82,84,86,88,90,92,94,96,98,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,440,480,520,560,600,800};

    TH1D *hist_METl_bins, *hist_pt_bins, *hist_MET_bins;
    void init_binning_histograms() {
        hist_METl_bins = new TH1D("hist_METl_bins","",n_METl_bins,METl_bins);
        hist_pt_bins = new TH1D("hist_pt_bins","",n_pt_bins,pt_bins);
        hist_MET_bins = new TH1D("hist_MET_bins","",n_pt_bins,MET_bins); //hist_MET_bins->SetStats(0);
    }
}

//-----------
// FUNCTIONS
//-----------

template<class variableType>
void SetInputBranch(TTree* inputTree, string branchName, variableType variablePointer) {
    // set variablePointer to reference value in branch if it exists
    for (auto branch : *(inputTree->GetListOfBranches())) {
        if (branch->GetName() == branchName) {
            inputTree->SetBranchStatus(branchName.c_str(), 1);
            inputTree->SetBranchAddress(branchName.c_str(), variablePointer);
        }
    }
}

template<class variableType>
void CopyBranch(TTree* inputTree, TTree* outputTree, string inputBranchName, string outputBranchName, variableType variablePointer, string varType) {
    // copy branch and set variablePointer to reference value in branch if it exists
    for (auto branch : *(inputTree->GetListOfBranches())) {
        if (branch->GetName() == inputBranchName) {
            inputTree->SetBranchStatus(inputBranchName.c_str(), 1);
            inputTree->SetBranchAddress(inputBranchName.c_str(), variablePointer);
            if (varType.compare(0, 6, "vector") == 0)
                outputTree->Branch(outputBranchName.c_str(), varType.c_str(), variablePointer);
            else
                outputTree->Branch(outputBranchName.c_str(), variablePointer, (outputBranchName+"/"+varType).c_str());
        }
    }
}

vector<int> int_copy_vars;
vector<double> double_copy_vars;
vector<float> float_copy_vars;

void CopyAllBranches(TTree* inputTree, TTree* outputTree, vector<string> branches) {
    int_copy_vars.clear(); int_copy_vars.reserve(branches.size());
    double_copy_vars.clear(); double_copy_vars.reserve(branches.size());
    float_copy_vars.clear(); float_copy_vars.reserve(branches.size());
    for (string branch : branches) {
        string branch_name = branch.substr(0,branch.length()-2);
        char branch_type = branch.substr(branch.length()-1)[0];
        switch (branch_type) {
            case 'I':
                int int_copy_var; int_copy_vars.push_back(int_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &int_copy_vars.back(), "I");
                break;
            case 'D':
                double double_copy_var; double_copy_vars.push_back(double_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &double_copy_vars.back(), "D");
                break;
            case 'F':
                float float_copy_var; float_copy_vars.push_back(float_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &float_copy_vars.back(), "F");
                break;
            default:
                cout << "Unknown branch type" << endl;
        }
    }
}

//  period: data15-16 (input) -> ZMC16a (source file), data17 -> ZMC16cd, data18 -> ZMC16e
float GetLumi(TString period) {
    float lumi = 1.0;
    if (period.Contains("mc16e") || period.Contains("data18")) lumi = 59900;
    else if (period.Contains("mc16cd") || period.Contains("data17")) lumi = 44000;
    else if (period.Contains("mc16a") || period.Contains("data15-16")) lumi = 36100;
    return lumi;
}

string getMCPeriod(TString period) {
    string mc_period;
    if (period == "data15-16") mc_period = "mc16a";
    else if (period == "data17") mc_period = "mc16cd";
    else if (period == "data18") mc_period = "mc16e";
    else mc_period = period;
    return mc_period;
}

string getMCPeriod(string period) {
    return string(getMCPeriod(TString(period)));
}

string DataPeriod(TString period) {
    string data_period;
    if (period == "mc16a") data_period = "data15-16";
    else if (period == "mc16cd") data_period = "data17";
    else if (period == "mc16e") data_period = "data18";
    else data_period = period;
    return data_period;
}

template <typename T> constexpr string_view type_name() {
    // from https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/56766138#56766138
    // e.g. cout << type_name<decltype(weighted_dataframe)>() << endl;
    std::string_view name, prefix, suffix;
#ifdef __clang__
    name = __PRETTY_FUNCTION__;
    prefix = "std::string_view type_name() [T = ";
    suffix = "]";
#elif defined(__GNUC__)
    name = __PRETTY_FUNCTION__;
    prefix = "constexpr std::string_view type_name() [with T = ";
    suffix = "; std::string_view = std::basic_string_view<char>]";
#elif defined(_MSC_VER)
    name = __FUNCSIG__;
    prefix = "class std::basic_string_view<char,struct std::char_traits<char> > __cdecl type_name<";
    suffix = ">(void)";
#endif
    name.remove_prefix(prefix.size());
    name.remove_suffix(suffix.size());
    return name;
}

void passTest(string msg) {
    cout << BOLD(PGRN("PASSED TEST: ")) << msg << endl;
}

void failTest(string msg) {
    cout << BOLD(PRED("ERROR: ")) << msg << endl;
    exit(0);
}

//-------------------------
// TREE MANIPULATION CLASS
//-------------------------

using BranchRenameOptions = vector<tuple<string, string>>;
using BranchAddOptions = vector<tuple<string, string>>;

class TreeCreator {
public:
    ROOT::RDataFrame *dataframe;
    string out_file_name;
    string out_tree_name;
    string cut;
    string final_cut;
    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;
    BranchAddOptions branches_to_add;

    TreeCreator() {
    }

    void read(string file_name, string tree_name) {
        cout << "Opening read file      : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;

        this->dataframe = new ROOT::RDataFrame(tree_name, file_name);
    }

    void setBranchesToCopy(vector<string> branches_to_copy) {
        this->branches_to_copy = branches_to_copy;
    }

    void setBranchesToRename(BranchRenameOptions branches_to_rename) {
        this->branches_to_rename = branches_to_rename;
    }

    void setBranchesToAdd(BranchAddOptions branches_to_add) {
        this->branches_to_add = branches_to_add;
    }

    void setCut(string cut) {
        this->cut = cut;
    }

    void setFinalCut(string cut) {
        this->final_cut = cut;
    }

    void write(string file_name, string tree_name) {
        this->out_file_name = file_name;
        this->out_tree_name = tree_name;

        cout << "Opening write file     : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;
        cout << endl;

        cout << "Processing" << endl;

        //--- apply cut
        auto reduced_dataframe = this->dataframe->Filter(this->cut.c_str());

        //--- get all branches to save
        vector<string> all_out_branches = this->branches_to_copy;

        //--- rename branches
        for (auto branch : this->branches_to_rename) {
            string old_name = get<0>(branch);
            string new_name = get<1>(branch);
            reduced_dataframe = reduced_dataframe.Define(new_name.c_str(), old_name.c_str());
            all_out_branches.push_back(new_name);
        }

        //--- add branches
        for (auto branch : this->branches_to_add) {
            string branch_name = get<0>(branch);
            //string expression = get<1>(branch);
            string call = get<1>(branch);
            //gInterpreter->Declare(expression.c_str());
            reduced_dataframe = reduced_dataframe.Define(branch_name.c_str(), call.c_str());
            all_out_branches.push_back(branch_name);
        }

        //--- apply cut
        if (!final_cut.empty())
            reduced_dataframe = reduced_dataframe.Filter(this->final_cut.c_str());

        reduced_dataframe.Snapshot(this->out_tree_name.c_str(), out_file_name.c_str(), all_out_branches);
        cout << endl;
    }
};

//---------
// OPTIONS
//---------

struct GlobalOptions {
    bool is_photon; // vs. bkg
    bool is_data; // vs. MC
    string period;
    string sampleID;

    string photon_mc_path;
    string photon_data_path;
    string bkg_mc_path;
    string bkg_data_path;

    string my_samples_folder;
    string sampling_method;
    string reduction_folder;
    string smearing_folder;
    string reweighting_folder;
    string plots_folder;

    string unit_test_folder;

    string channel;
    string type;

    string save_tree_name;
};

struct ReductionOptions {
    string in_file_name;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    string unit_test_folder;

    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;
    BranchAddOptions branches_to_add;

    string cut;
    string final_cut;

    bool unit_testing;
};

struct SmearingOptions {
    string in_file_name;
    string in_file_path;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    string unit_test_folder;

    string period;
    string data_period;
    string mc_period;
    bool is_data;
    string channel;

    vector<string> branches_to_copy;
    BranchAddOptions branches_to_add;

    bool unit_testing;
    bool turn_off_shifting_and_smearing;
    bool diagnostic_plots;
};

struct ReweightingOptions {
    string in_file_name;
    string in_file_path;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    string period;
    string data_period;
    string mc_period;
    bool is_data;
    string channel;
    string reweight_var;

    string reduction_folder;

    vector<string> branches_to_copy;
    BranchAddOptions branches_to_add;

    bool unit_testing;
    bool turn_off_shifting_and_smearing;
    bool diagnostic_plots;
};

struct PlottingOptions {
    string in_file_name;
    string in_file_path;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    string period;
    string data_period;
    string mc_period;
    bool is_data;
    string reweight_var;

    vector<string> regions;
    vector<string> plot_features;
    vector<string> channels;

    bool blinded;
    bool print_photon_yield_only;

    string reduction_folder;
    string reweighting_folder;
    string plots_folder;

    vector<string> branches_to_copy;
    BranchAddOptions branches_to_add;

    bool unit_testing;
    bool turn_off_shifting_and_smearing;
    bool diagnostic_plots;
};

#endif
