#include "../../Common/Settings.C"

using namespace std;

void higgsino_cut_efficiency_grid() {

    TCut baseline_cut("nJet30>=2 && jetPt[0]>30 && jetPt[1]>30 && nLep_signal>= 2 && lepPt[0]>25 && lepPt[1]>25 && met_Et>200");
    TCut additional_cut("Ptll>50");

    string filePath = "/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.7/SUSY2/SUSY2_Signal_mc16e/";

    // SLN1 is gluino-slepton; do ZN1 (gluino-Z) when available
    string fileNames [] = {"GG_N2_SLN1_1000_100_merged_processed.root", "GG_N2_SLN1_1000_300_merged_processed.root", "GG_N2_SLN1_1000_500_merged_processed.root", "GG_N2_SLN1_1000_700_merged_processed.root", "GG_N2_SLN1_1000_800_merged_processed.root", "GG_N2_SLN1_1000_840_merged_processed.root", "GG_N2_SLN1_1000_880_merged_processed.root", "GG_N2_SLN1_1000_900_merged_processed.root", "GG_N2_SLN1_1000_920_merged_processed.root", "GG_N2_SLN1_1000_940_merged_processed.root", "GG_N2_SLN1_1000_960_merged_processed.root", "GG_N2_SLN1_1000_970_merged_processed.root", "GG_N2_SLN1_1200_1000_merged_processed.root", "GG_N2_SLN1_1200_100_merged_processed.root", "GG_N2_SLN1_1200_1040_merged_processed.root", "GG_N2_SLN1_1200_1080_merged_processed.root", "GG_N2_SLN1_1200_1100_merged_processed.root", "GG_N2_SLN1_1200_1120_merged_processed.root", "GG_N2_SLN1_1200_1140_merged_processed.root", /*"GG_N2_SLN1_1200_1160_merged_processed.root", "GG_N2_SLN1_1200_1170_merged_processed.root",*/ "GG_N2_SLN1_1200_300_merged_processed.root", "GG_N2_SLN1_1200_500_merged_processed.root", "GG_N2_SLN1_1200_700_merged_processed.root", "GG_N2_SLN1_1200_900_merged_processed.root", "GG_N2_SLN1_1400_100_merged_processed.root", "GG_N2_SLN1_1400_1100_merged_processed.root", "GG_N2_SLN1_1400_1200_merged_processed.root", "GG_N2_SLN1_1400_1240_merged_processed.root", "GG_N2_SLN1_1400_1280_merged_processed.root", "GG_N2_SLN1_1400_1300_merged_processed.root", "GG_N2_SLN1_1400_1320_merged_processed.root", "GG_N2_SLN1_1400_1340_merged_processed.root", "GG_N2_SLN1_1400_1360_merged_processed.root", "GG_N2_SLN1_1400_1370_merged_processed.root", "GG_N2_SLN1_1400_300_merged_processed.root", "GG_N2_SLN1_1400_500_merged_processed.root", "GG_N2_SLN1_1400_700_merged_processed.root", "GG_N2_SLN1_1400_900_merged_processed.root", "GG_N2_SLN1_1600_100_merged_processed.root", "GG_N2_SLN1_1600_1100_merged_processed.root", "GG_N2_SLN1_1600_1300_merged_processed.root", "GG_N2_SLN1_1600_1400_merged_processed.root", "GG_N2_SLN1_1600_1440_merged_processed.root", "GG_N2_SLN1_1600_1480_merged_processed.root", "GG_N2_SLN1_1600_1500_merged_processed.root", "GG_N2_SLN1_1600_1520_merged_processed.root", "GG_N2_SLN1_1600_1540_merged_processed.root", "GG_N2_SLN1_1600_1560_merged_processed.root", /*"GG_N2_SLN1_1600_1570_merged_processed.root",*/ "GG_N2_SLN1_1600_300_merged_processed.root", "GG_N2_SLN1_1600_500_merged_processed.root", "GG_N2_SLN1_1600_700_merged_processed.root", "GG_N2_SLN1_1600_900_merged_processed.root", "GG_N2_SLN1_1700_100_merged_processed.root", "GG_N2_SLN1_1700_1100_merged_processed.root", "GG_N2_SLN1_1700_300_merged_processed.root", "GG_N2_SLN1_1700_500_merged_processed.root", "GG_N2_SLN1_1700_700_merged_processed.root", "GG_N2_SLN1_1700_900_merged_processed.root", "GG_N2_SLN1_1800_100_merged_processed.root", "GG_N2_SLN1_1800_1100_merged_processed.root", "GG_N2_SLN1_1800_1300_merged_processed.root", "GG_N2_SLN1_1800_1500_merged_processed.root", "GG_N2_SLN1_1800_1640_merged_processed.root", "GG_N2_SLN1_1800_1680_merged_processed.root", "GG_N2_SLN1_1800_1700_merged_processed.root", "GG_N2_SLN1_1800_1720_merged_processed.root", "GG_N2_SLN1_1800_1740_merged_processed.root", "GG_N2_SLN1_1800_1760_merged_processed.root", "GG_N2_SLN1_1800_1770_merged_processed.root", "GG_N2_SLN1_1800_300_merged_processed.root", "GG_N2_SLN1_1800_500_merged_processed.root", "GG_N2_SLN1_1800_700_merged_processed.root", "GG_N2_SLN1_1800_900_merged_processed.root", "GG_N2_SLN1_1900_100_merged_processed.root", "GG_N2_SLN1_1900_1100_merged_processed.root", "GG_N2_SLN1_1900_300_merged_processed.root", "GG_N2_SLN1_1900_500_merged_processed.root", "GG_N2_SLN1_1900_700_merged_processed.root", "GG_N2_SLN1_1900_900_merged_processed.root", "GG_N2_SLN1_2000_100_merged_processed.root", "GG_N2_SLN1_2000_1100_merged_processed.root", "GG_N2_SLN1_2000_1300_merged_processed.root", "GG_N2_SLN1_2000_1500_merged_processed.root", "GG_N2_SLN1_2000_1700_merged_processed.root", "GG_N2_SLN1_2000_1900_merged_processed.root", "GG_N2_SLN1_2000_300_merged_processed.root", "GG_N2_SLN1_2000_500_merged_processed.root", "GG_N2_SLN1_2000_700_merged_processed.root", "GG_N2_SLN1_2000_900_merged_processed.root", "GG_N2_SLN1_2200_100_merged_processed.root", "GG_N2_SLN1_2200_1100_merged_processed.root", "GG_N2_SLN1_2200_1300_merged_processed.root", "GG_N2_SLN1_2200_1500_merged_processed.root", "GG_N2_SLN1_2200_1700_merged_processed.root", "GG_N2_SLN1_2200_1900_merged_processed.root", "GG_N2_SLN1_2200_300_merged_processed.root", "GG_N2_SLN1_2200_500_merged_processed.root", "GG_N2_SLN1_2200_700_merged_processed.root", "GG_N2_SLN1_2200_900_merged_processed.root", "GG_N2_SLN1_2400_100_merged_processed.root", "GG_N2_SLN1_2400_1100_merged_processed.root", "GG_N2_SLN1_2400_1300_merged_processed.root", "GG_N2_SLN1_2400_1500_merged_processed.root", "GG_N2_SLN1_2400_300_merged_processed.root", "GG_N2_SLN1_2400_500_merged_processed.root", "GG_N2_SLN1_2400_700_merged_processed.root", "GG_N2_SLN1_2400_900_merged_processed.root", "GG_N2_SLN1_600_100_merged_processed.root", "GG_N2_SLN1_600_300_merged_processed.root", "GG_N2_SLN1_600_400_merged_processed.root", "GG_N2_SLN1_600_440_merged_processed.root", "GG_N2_SLN1_600_480_merged_processed.root", "GG_N2_SLN1_600_500_merged_processed.root", "GG_N2_SLN1_600_520_merged_processed.root", "GG_N2_SLN1_600_540_merged_processed.root", "GG_N2_SLN1_600_560_merged_processed.root", "GG_N2_SLN1_600_570_merged_processed.root", "GG_N2_SLN1_800_100_merged_processed.root", "GG_N2_SLN1_800_300_merged_processed.root", "GG_N2_SLN1_800_500_merged_processed.root", "GG_N2_SLN1_800_600_merged_processed.root", "GG_N2_SLN1_800_640_merged_processed.root", "GG_N2_SLN1_800_680_merged_processed.root", "GG_N2_SLN1_800_700_merged_processed.root", "GG_N2_SLN1_800_720_merged_processed.root", "GG_N2_SLN1_800_740_merged_processed.root", "GG_N2_SLN1_800_760_merged_processed.root", "GG_N2_SLN1_800_770_merged_processed.root"};

    vector<tuple<int, int, float>> results;
    for (string fileName : fileNames) {
        TFile* inputFile = TFile::Open((filePath + fileName).c_str());
        stringstream unsplitName(fileName);
        string sampleID;
        int slepton_mass, neutralino_mass;
        for (int i=0; i<5; i++) {
            string substring;
            getline(unsplitName, substring, '_');
            sampleID += (substring);
            if (i==3) slepton_mass = stoi(substring);
            if (i==4) neutralino_mass = stoi(substring);
            if (i<4) sampleID += "_";
        }
        string treeName = sampleID + "_NoSys";
        TTree* inputTree = (TTree*)inputFile->Get(treeName.c_str());

        inputTree->Draw(">> baseline_cut_event_list", baseline_cut);
        inputTree->Draw(">> additional_cut_event_list", baseline_cut + additional_cut);
        TEventList *baseline_cut_event_list = (TEventList*)gDirectory->Get("baseline_cut_event_list");
        TEventList *additional_cut_event_list = (TEventList*)gDirectory->Get("additional_cut_event_list");
        int n_baseline_cut_passed = baseline_cut_event_list->GetN();
        int n_additional_cut_passed = additional_cut_event_list->GetN();
        float efficiency = float(n_additional_cut_passed) / n_baseline_cut_passed;

        tuple<int, int, float> result = make_tuple(slepton_mass, neutralino_mass, efficiency);
        results.push_back(result);
    }

    auto canvas = new TCanvas("canvas", "canvas", 700, 700);
    gStyle->SetOptStat(0);
    auto hist_grid = new TH2F("hist_grid","hist_grid",25,0,2500,200,0,2000);
    auto color_grid = new TH2F("color_grid","Efficiency of Ptll>50 Cut on SLN1 Grid",25,0,2500,200,0,2000);
    for (auto result : results) {
        color_grid->Fill(get<0>(result), get<1>(result), get<2>(result));
        hist_grid->Fill(get<0>(result), get<1>(result), get<2>(result));
    }
    color_grid->GetXaxis()->SetTitle("Slepton Mass");
    color_grid->GetYaxis()->SetTitle("Neutralino Mass");
    canvas->SetLeftMargin(0.15);
    color_grid->GetYaxis()->SetTitleOffset(2.0);
    color_grid->Draw("COL");
    TLine *midline = new TLine(0,0,2000,2000);
    midline->SetLineColor(kRed);
    midline->SetLineWidth(2);
    midline->Draw();
    hist_grid->SetMarkerSize(1.8);
    hist_grid->SetBarOffset(0.2);
    hist_grid->SetMarkerSize(0.7);
    gStyle->SetPaintTextFormat("0.2f");
    hist_grid->Draw("TEXT SAME");
    canvas->Print("higgsino_efficiency_grid.eps");
}
