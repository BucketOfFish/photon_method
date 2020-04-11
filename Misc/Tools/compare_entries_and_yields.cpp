#include "../../Main/Settings.cpp"

void compare_entries_and_yields() {
    //--- settings
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
    vector<string> filenames {
        "data15-16_data_bkg.root", "data17_data_bkg.root", "data18_data_bkg.root",
        "data15-16_data_photon.root", "data17_data_photon.root", "data18_data_photon.root",
        "mc16a_SinglePhoton222.root", "mc16cd_SinglePhoton222.root", "mc16e_SinglePhoton222.root",
        "mc16a_ttbar.root", "mc16cd_ttbar.root", "mc16e_ttbar.root",
        "mc16a_triboson.root", "mc16cd_triboson.root", "mc16e_triboson.root",
        "mc16a_singleTop.root", "mc16cd_singleTop.root", "mc16e_singleTop.root",
        "mc16a_topOther.root", "mc16cd_topOther.root", "mc16e_topOther.root",
        "mc16a_lowMassDY.root", "mc16cd_lowMassDY.root", "mc16e_lowMassDY.root",
        "mc16a_Wjets.root", "mc16cd_Wjets.root", "mc16e_Wjets.root",
        "mc16a_Zjets.root", "mc16cd_Zjets.root", "mc16e_Zjets.root",
        "mc16a_diboson.root", "mc16cd_diboson.root", "mc16e_diboson.root",
        "mc16a_higgs.root", "mc16cd_higgs.root", "mc16e_higgs.root",
    };

    //--- compare entries
    cout << "Comparing n entries" << endl;
    cout << endl;

    vector<string> matches;
    vector<string> non_matches;
    for (string filename : filenames) {
        TFile *old_tfile = TFile::Open((old_file_path + filename).c_str());
        TTree *old_ttree = (TTree*)old_tfile->Get("BaselineTree");
        TFile *new_tfile = TFile::Open((new_file_path + filename).c_str());
        TTree *new_ttree = (TTree*)new_tfile->Get("BaselineTree");
        if (old_ttree->GetEntries() == new_ttree->GetEntries())
            matches.push_back(filename);
        else
            non_matches.push_back(filename);
    }

    cout << "Matching n entries:" << endl;
    for (string filename : matches) {
        cout << "\t" << filename << endl;
    }
    cout << "Not matching n entries:" << endl;
    for (string filename : non_matches) {
        cout << "\t" << filename << endl;
    }
    cout << endl;

    //--- compare yields
    cout << "Comparing yields for ntuples with mismatched n entries" << endl;
    cout << endl;

    for (string filename : non_matches) {
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
        TH1F* hz = new TH1F("hz", "", 100, 0, 3);
        old_ttree->Draw("mll>>hz", cuts::plot_regions[region]*"totalWeight", "goff");
        cout << hz->Integral() << endl;
        //cout << old_ttree->GetEntries(cuts::plot_regions[region]*"totalWeight") << endl;
        cout << "\t" << "new " << region << " entries: ";
        new_ttree->Draw("mll>>hz", cuts::plot_regions[region]*"totalWeight", "goff");
        cout << hz->Integral() << endl;
        //cout << new_ttree->GetEntries(cuts::plot_regions[region]) << endl;
    }
}
