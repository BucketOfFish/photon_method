#include "../../../Main/Settings.cpp"

void process_systematics_table(vector<string> files, string process) {
    // Options
    TChain chain((process + "_NoSys").c_str());
    for (auto file : files)
        chain.Add(file.c_str());

    //vector<string> regions = {"SRC", "SRLow", "SRMed", "SRHigh", "SRLowZ", "SRMedZ", "SRHighZ",
                              ////"CRC", "CRLow", "CRMed", "CRHigh", "CRLowZ", "CRMedZ", "CRHighZ", 
                              //"VRC", "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ"}; 
    vector<string> regions = {"CRZ", "SRLow_1", "SRLow_2", "SRLow2", "SROffShell_2", "EWK_VRLow", "VRLow2",
                              "VROffShell"};

    TCut lumi = "(RandomRunNumber<320000 ? 36200 : (RandomRunNumber>320000 && RandomRunNumber<348000) ? 44300 : 58500)";
    TCut MC_weight = "genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*globalDiLepTrigSF";

    vector<string> systematics = {"scale", "PDF"};
    if (process == "ttbar")
        systematics = {"ttbar_scale", "ttbar_PDF"};

    // Table Production
    for (auto systematic : systematics) {
        vector<string> LHE3Weights;
        if (systematic == "scale") {
            LHE3Weights.push_back("MUR1_MUF1_PDF261000");
            LHE3Weights.push_back("MUR0p5_MUF0p5_PDF261000");
            LHE3Weights.push_back("MUR0p5_MUF1_PDF261000");
            LHE3Weights.push_back("MUR1_MUF0p5_PDF261000");
            LHE3Weights.push_back("MUR1_MUF2_PDF261000");
            LHE3Weights.push_back("MUR2_MUF1_PDF261000");
            LHE3Weights.push_back("MUR2_MUF2_PDF261000");
        }
        else if (systematic == "PDF") {
            LHE3Weights.push_back("MUR1_MUF1_PDF261000");
            LHE3Weights.push_back("MUR1_MUF1_PDF13000");
            LHE3Weights.push_back("MUR1_MUF1_PDF25300");
            //LHE3Weights.push_back("MUR1_MUF1_PDF269000");
            //LHE3Weights.push_back("MUR1_MUF1_PDF270000");
            //for (int i=261001;i<261101;i++) LHE3Weights.push_back("MUR1_MUF1_PDF" + to_string(i));
        }
        else if (systematic == "ttbar_scale") {
            LHE3Weights.push_back("nominal");
            LHE3Weights.push_back("muR0p5,muF0p5");
            LHE3Weights.push_back("muR0p5,muF1");
            LHE3Weights.push_back("muR1,muF0p5");
            LHE3Weights.push_back("muR0p5,muF2");
            LHE3Weights.push_back("muR2,muF0p5");
            LHE3Weights.push_back("muR1,muF2");
            LHE3Weights.push_back("muR2,muF1");
            LHE3Weights.push_back("muR2,muF2");
        }
        else if (systematic == "ttbar_PDF") {
            LHE3Weights.push_back("nominal");
            LHE3Weights.push_back("PDFset265000");
            LHE3Weights.push_back("PDFset266000");
        }

        ofstream myfile;
        myfile.open((process + ".txt").c_str());

        for (auto region : regions) {
            TCut cut = cuts::selections[region];
            myfile << "Finding Systematics for Region " << region << endl;

            vector<double> yields;
            for (auto LHE3Weight : LHE3Weights) {
                myfile << "Working on Weight " << LHE3Weight << endl;
                TH1F* hist = new TH1F(("hists_" + LHE3Weight).c_str(),"",1,0,1);
                TCut LHE_branch = ("LHE3Weight_" + LHE3Weight).c_str();
                chain.Draw(("mll>>hists_" + LHE3Weight).c_str(), cut*lumi*MC_weight*LHE_branch, "goff");
                float yield = fabs(hist->Integral(0,2));
                myfile << "region yield: " << yield << endl;
                yields.push_back(yield);
                delete hist;
            }

            float uncertainty_max = 0;
            for (int i=0; i<yields.size(); i++) {
                if (i==0) continue;
                float uncertainty = abs((yields[i]-yields[0]) / yields[0]);
                if (uncertainty>uncertainty_max) uncertainty_max = uncertainty;
            }
            myfile << "Region " << region << " " << systematic << " uncertainty: " << uncertainty_max << endl;
        }	

        myfile.close();
    }
}

void bin_yield_systematics_table() {
    vector<string> processes = {"Zjets", "diboson", "ttbar"};
    string mc_path = "/public/data/SUSY_Systematics/Skimmed/EWKPreselection/";
    for (auto process : processes) {
        vector<string> files;
        files.push_back(mc_path + "SUSY2_Bkgs_mc16a/" + process + "_merged_processed.root");
        files.push_back(mc_path + "SUSY2_Bkgs_mc16cd/" + process + "_merged_processed.root");
        files.push_back(mc_path + "SUSY2_Bkgs_mc16e/" + process + "_merged_processed.root");
        process_systematics_table(files, process);
    }
}
