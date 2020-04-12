#include "../../Main/Settings.cpp"

//--- helper function
tuple<string, string, string> getPlotRegionInfo(string channel, string region) {
    if (cuts::plot_regions.count(region) == 0) {
        cout << "Unrecognized region! Exiting." << endl;
        exit(0);
    }
    TCut plot_region = cuts::plot_regions[region];

    if (TString(channel).EqualTo("ee")) plot_region += cuts::ee;
    else if (TString(channel).EqualTo("mm")) plot_region += cuts::mm;
    else if (TString(channel).EqualTo("em")) plot_region += cuts::em;
    else if (TString(channel).EqualTo("me")) plot_region += cuts::me;
    else if (TString(channel).EqualTo("SF")) plot_region += cuts::SF;
    else if (TString(channel).EqualTo("DF")) plot_region += cuts::DF;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    TCut plot_CR = plot_region + cuts::CR;
    plot_region += cuts::plot_region_met_portions[region];

    string region_name = region + " " + channel;

    return make_tuple(region_name, string(plot_region), string(plot_CR));
}

void get_yields() {
    //--- settings
    //vector<string> regions = vector<string>{"VRC", "VRLow23", "VRLow4", "VRMed23",
                    //"VRMed4", "VRHigh23", "VRHigh4", "VRDPhiLow6", "VRDPhiMed6", "VRDPhiHigh6", "VR3L"};
    vector<string> regions = vector<string>{"VRDPhiLow6", "VRDPhiMed6", "VRDPhiHigh6", "VR3L"};
    vector<string> channels = vector<string>{"SF"};
    bool is_data = true;
    bool blinded = false;

    //--- get photon and data files
    vector<string> periods = vector<string>{"data15-16", "data17", "data18"};
    TChain *tch_photon = new TChain("BaselineTree");
    TChain *tch_data = new TChain("BaselineTree");

    for (auto period : periods) {
        string photon_path = "/public/data/Photon/Samples/HistogramSampling/ReweightedNtuples/";
        string data_period = DataPeriod(period);
        string mc_period = getMCPeriod(period);
        string photon_ee_filename, photon_mm_filename;
        if (is_data) {
            photon_ee_filename = photon_path + data_period + "_data_photon_ee.root";
            photon_mm_filename = photon_path + data_period + "_data_photon_mm.root";
        }
        else {
            photon_ee_filename = photon_path + mc_period + "_SinglePhoton222_ee.root";
            photon_mm_filename = photon_path + mc_period + "_SinglePhoton222_mm.root";
        }
        tch_photon->Add(photon_ee_filename.c_str());
        tch_photon->Add(photon_mm_filename.c_str());

        string data_filename = "/public/data/Photon/Samples/ReducedNtuples/" + data_period + "_data_bkg.root";
        tch_data->Add(data_filename.c_str());
    }

    //--- get yields
    for (string region : regions) {
        for (string channel : channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            cout << region << " " << channel << endl;
            cout << plot_region << endl;
            cout << "\t" << "Data entries: ";
            if (blinded && (region.find("SR") != std::string::npos))
                cout << "BLINDED" << endl;
            else
                cout << tch_data->GetEntries(plot_region.c_str()) << endl;
            //cout << "\t" << "Photon entries: " << tch_photon->GetEntries(plot_region.c_str()) << endl;
            //cout << "\t" << "Data yield: ";
            //if (blinded && (region.find("SR") != std::string::npos))
                //cout << "BLINDED" << endl;
            //else
                //cout << tch_data->GetEntries(plot_region.c_str()) << endl;
            //cout << "\t" << "Photon yield: " << tch_photon->GetEntries(plot_region.c_str()) << endl;
            //cout << "\t" << "Data yield (CR): " << tch_photon->GetEntries(plot_CR.c_str()) << endl;
            //cout << "\t" << "Photon yield (CR): " << tch_photon->GetEntries(plot_CR.c_str()) << endl;
            //cout << endl;
        }
    }
}
