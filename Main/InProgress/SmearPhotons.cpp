#include "../Settings.cpp"
#include <fstream>

using namespace std;
using SmearingDists = map<int, pair<float, float>>;
using METlHists = unordered_map<int, ROOT::RDF::RResultPtr<TH1D>>;
using MllHists = unordered_map<int, unordered_map<int, ROOT::RDF::RResultPtr<TH1D>>>;

//-----------
// COPY TREE
//-----------

void copyTree(SmearingOptions options) {
    TreeCreator *tree_copier = new TreeCreator();

    tree_copier->read(options.in_file_name, options.in_tree_name);

    tree_copier->setBranchesToCopy(options.branches_to_copy);
    //tree_copier->setBranchesToRename(options.branches_to_rename);
    tree_copier->setBranchesToAdd(options.branches_to_add);

    //tree_copier->setCut(options.cut);
    tree_copier->setCut("1");

    tree_copier->write(options.out_file_name, options.out_tree_name);

    TFile *in_file = new TFile(options.out_file_name.c_str());
    TTree* BaselineTree = (TTree*)in_file->Get(options.out_tree_name.c_str());
}

//----------
// SMEARING
//----------

class Smearer {
public:
    SmearingOptions options;
    METlHists metl_bin_pt;
    MllHists mll_bin_pt_metl;

    Smearer(SmearingOptions options) {
        this->options = options;
    };

    void setupSmearingHistograms() {
        string data_file_name = this->options.in_file_path + options.data_period + "_data_bkg.root";
        string ttbar_file_name = this->options.in_file_path + options.mc_period + "_ttbar.root";
        string diboson_file_name = this->options.in_file_path + options.mc_period + "_diboson.root";
        string Z_file_name = this->options.in_file_path + options.mc_period + "_Zjets.root";
        string photon_file_name;
        if (options.is_data)
            photon_file_name = this->options.in_file_path + options.mc_period + "_SinglePhoton222.root";
        else
            photon_file_name = this->options.in_file_path + options.data_period + "_data_photon.root";

        TChain *tch_all = new TChain("BaselineTree");
        tch_all->Add(data_file_name.c_str());
        tch_all->Add(ttbar_file_name.c_str());
        tch_all->Add(diboson_file_name.c_str());
        tch_all->Add(Z_file_name.c_str());
        tch_all->Add(photon_file_name.c_str());
        ROOT::RDataFrame *dataframe = new ROOT::RDataFrame(*tch_all);

        ROOT::RDF::TH1DModel hist_metl = ROOT::RDF::TH1DModel("", "E_{T,||}^{miss} [GeV]", bins::n_smearing_bins, bins::smearing_low, bins::smearing_high);
        ROOT::RDF::TH1DModel hist_mll = ROOT::RDF::TH1DModel("", "m_{ll} [GeV]", bins::n_mll_bins, bins::mll_bin);

        string plot_region = "nJet30>=2";

        for (int i=0; i<bins::n_pt_bins; i++) {
            string plot_region = "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && Ptll>" + to_string(bins::pt_bins[i]) + " && Ptll<" + to_string(bins::pt_bins[i+1]);
            this->metl_bin_pt[i] = dataframe->Filter(plot_region).Histo1D(hist_metl, "met_Et", "totalWeight");
            for (int j=0; j<bins::n_METl_bins; j++) {
                plot_region += "METl>" + to_string(bins::METl_bins[i]) + " && METl<" + to_string(bins::METl_bins[i+1]);
                this->mll_bin_pt_metl[i][j] = dataframe->Filter(plot_region).Histo1D(hist_mll, "mll", "totalWeight");
            }
        }
    }
};

//------
// MAIN
//------

void runSmearing(SmearingOptions options) {
    copyTree(options);
    Smearer *smearer = new Smearer(options);
}

//---------------
// MAIN FUNCTION
//---------------

void SmearPhotons(SmearingOptions options) {
    runSmearing(options);
}
