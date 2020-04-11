#include "../../Main/Settings.cpp"

void compare_yields() {
    //--- settings
    //vector<string> regions = vector<string>{"VRDPhiLow6", "VRDPhiMed6", "VRDPhiHigh6", "VR3L"};
    vector<string> regions = vector<string>{"VRDPhiLow6"};
    //string old_file_path = "/public/data/Photon/Samples/HistogramSampling/ReweightedNtuples/";
    //string new_file_path = "/public/data/Photon/NewSamples/HistogramSampling/ReweightedNtuples/";
    //vector<string> filenames{
        //"data15-16_data_photon_ee.root",
        //"data15-16_data_photon_mm.root",
        //"data17_data_photon_ee.root",
        //"data17_data_photon_mm.root",
        //"data18_data_photon_ee.root",
        //"data18_data_photon_mm.root",
    //};
    string old_file_path = "/public/data/Photon/Samples/ReducedNtuples/";
    string new_file_path = "/public/data/Photon/NewSamples/ReducedNtuples/";
    vector<string> filenames{
        "data15-16_data_bkg.root", "mc16a_singleTop.root", "mc16cd_SinglePhoton222.root", "mc16e_lowMassDY.root",
        "data15-16_data_photon.root", "mc16a_topOther.root", "mc16cd_singleTop.root", "mc16e_SinglePhoton222.root",
        "data17_data_bkg.root", "mc16a_triboson.root", "mc16cd_topOther.root", "mc16e_singleTop.root",
        "data17_data_photon.root", "mc16a_ttbar.root", "mc16cd_triboson.root", "mc16e_topOther.root",
        "data18_data_bkg.root", "mc16cd_ttbar.root", "mc16e_triboson.root",
        "data18_data_photon.root", "mc16a_Wjets.root", "mc16e_ttbar.root",
        "mc16a_diboson.root", "mc16a_Zjets.root", "mc16cd_Wjets.root",
        "mc16a_higgs.root", "mc16cd_diboson.root", "mc16cd_Zjets.root", "mc16e_Wjets.root",
        "mc16a_lowMassDY.root", "mc16cd_higgs.root", "mc16e_diboson.root", "mc16e_Zjets.root",
        "mc16a_SinglePhoton222.root", "mc16cd_lowMassDY.root", "mc16e_higgs.root",
    };
    bool is_data = true;
    bool blinded = false;

    //--- get yields
    for (string filename : filenames) {
        cout << filename << endl;
        TFile *old_tfile = TFile::Open((old_file_path + filename).c_str());
        TTree *old_ttree = (TTree*)old_tfile->Get("BaselineTree");
        TFile *new_tfile = TFile::Open((new_file_path + filename).c_str());
        TTree *new_ttree = (TTree*)new_tfile->Get("BaselineTree");
        //for (string region : regions) {
            //cout << "\t" << "old " << region << " entries: ";
            //cout << old_ttree->GetEntries(cuts::plot_regions[region]) << endl;
            //cout << "\t" << "new " << region << " entries: ";
            //cout << new_ttree->GetEntries(cuts::plot_regions[region]) << endl;
        //}
        string region = (filename.find("hoton") != std::string::npos) ? "photon_comparison" : "bkg_baseline";
        cout << "\t" << "old " << region << " entries: ";
        cout << old_ttree->GetEntries(cuts::plot_regions[region]) << endl;
        cout << "\t" << "new " << region << " entries: ";
        cout << new_ttree->GetEntries(cuts::plot_regions[region]) << endl;
    }
}
