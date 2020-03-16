#include "../Settings.cpp"
#include <fstream>

using namespace std;
using SmearingDists = map<int, pair<float, float>>;
using METlHists = unordered_map<string, unordered_map<int, ROOT::RDF::RResultPtr<TH1D>>>;
using MllHists = unordered_map<string, unordered_map<int, unordered_map<int, ROOT::RDF::RResultPtr<TH1D>>>>;

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
    tree_copier->setCut("nJet30>=0");

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
    METlHists metl_bin_pt; // [process][pt_bin_n]
    MllHists mll_bin_pt_metl; // [process][pt_bin_n][metl_bin_n]
    unordered_map<string, ROOT::RDataFrame*> dataframes;

    Smearer(SmearingOptions options) {
        this->options = options;
    };

    void setupSmearingHistograms() {
        //--- fill metl_bin_pt and mll_bin_pt_metl
        string data_file_name = this->options.in_file_path + options.data_period + "_data_bkg.root";
        string ttbar_file_name = this->options.in_file_path + options.mc_period + "_ttbar.root";
        string diboson_file_name = this->options.in_file_path + options.mc_period + "_diboson.root";
        string Z_file_name = this->options.in_file_path + options.mc_period + "_Zjets.root";
        string photon_file_name;
        if (options.is_data)
            photon_file_name = this->options.in_file_path + options.mc_period + "_SinglePhoton222.root";
        else
            photon_file_name = this->options.in_file_path + options.data_period + "_data_photon.root";

        TChain *tch_data = new TChain("BaselineTree"); tch_data->Add(data_file_name.c_str());
        this->dataframes["data"] = new ROOT::RDataFrame(*tch_data);
        TChain *tch_ttbar = new TChain("BaselineTree"); tch_ttbar->Add(ttbar_file_name.c_str());
        this->dataframes["ttbar"] = new ROOT::RDataFrame(*tch_ttbar);
        TChain *tch_diboson = new TChain("BaselineTree"); tch_diboson->Add(diboson_file_name.c_str());
        this->dataframes["diboson"] = new ROOT::RDataFrame(*tch_diboson);
        TChain *tch_Z = new TChain("BaselineTree"); tch_Z->Add(Z_file_name.c_str());
        this->dataframes["Z"] = new ROOT::RDataFrame(*tch_Z);
        TChain *tch_photon = new TChain("BaselineTree"); tch_photon->Add(photon_file_name.c_str());
        this->dataframes["photon"] = new ROOT::RDataFrame(*tch_photon);

        ROOT::RDF::TH1DModel hist_metl = ROOT::RDF::TH1DModel("", "E_{T,||}^{miss} [GeV]", bins::n_smearing_bins, bins::smearing_low, bins::smearing_high);
        ROOT::RDF::TH1DModel hist_mll = ROOT::RDF::TH1DModel("", "m_{ll} [GeV]", bins::n_mll_bins, bins::mll_bin);

        string bkg_plot_region = "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25";
        string photon_plot_region = "nJet30>=2";

        vector<string> processes{"data", "ttbar", "diboson", "Z", "photon"};
        for (auto process : processes) {
            for (int i=0; i<bins::n_pt_bins; i++) {
                string plot_region;
                if (process == "photon") {
                    plot_region = photon_plot_region +
                        " && gamma_pt>" + to_string(bins::pt_bins[i]) +
                        " && gamma_pt<" + to_string(bins::pt_bins[i+1]);
                    this->metl_bin_pt[process][i] = this->dataframes[process]->Filter(plot_region)
                                                    .Histo1D(hist_metl, "METl_unsmeared", "totalWeight");
                }
                else {
                    plot_region = bkg_plot_region +
                        " && Ptll>" + to_string(bins::pt_bins[i]) +
                        " && Ptll<" + to_string(bins::pt_bins[i+1]);
                    this->metl_bin_pt[process][i] = this->dataframes[process]->Filter(plot_region)
                                                    .Histo1D(hist_metl, "METl", "totalWeight");
                    for (int j=0; j<bins::n_METl_bins; j++) {
                        plot_region += " && METl>" + to_string(bins::METl_bins[i]) +
                            " && METl<" + to_string(bins::METl_bins[i+1]);
                        this->mll_bin_pt_metl[process][i][j] = this->dataframes[process]->Filter(plot_region)
                                                               .Histo1D(hist_mll, "mll", "totalWeight");
                    }
                }
            }
        }
    }
};

//---------------
// MAIN FUNCTION
//---------------

void SmearPhotons(SmearingOptions options) {
    //copyTree(options);
    Smearer *smearer = new Smearer(options);
    smearer->setupSmearingHistograms();
}
