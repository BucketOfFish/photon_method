#include "../Common/Settings.C"

using namespace std;

class TreeReducer {
public:
    TFile *in_file, *out_file;    
    TTree *in_tree, *out_tree;    
    string out_tree_name;
    string cut;
    vector<string> branches_to_copy;
    vector<pair<string, string>> branches_to_rename;
    vector<string> new_branches;

    TreeReducer() {
    }

    void openReadFileGetTree(string file_name, string tree_name) {
        cout << "Opening read file      : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;

        this->in_file = TFile::Open(file_name.c_str());
        this->in_tree = (TTree*)this->in_file->Get(tree_name.c_str());

        cout << "Events in tree         : " << in_tree->GetEntries() << endl;
        cout << endl;
    }

    void openWriteFileSetTree(string file_name, string tree_name) {
        this->out_tree_name = tree_name;

        cout << "Opening write file     : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;
        cout << endl;

        this->out_file = TFile::Open(file_name.c_str(), "recreate");
        this->out_tree = this->in_tree->CloneTree(0);
    }

    void setBranchesToCopy(vector<string> branches_to_copy) {
        this->branches_to_copy = branches_to_copy;
    }

    void setBranchesToRename(vector<pair<string, string>> branches_to_rename) {
        this->branches_to_rename = branches_to_rename;
    }

    void setNewBranches(vector<string> new_branches) {
        this->new_branches = new_branches;
    }

    void setAllBranches() {
        this->in_tree->SetBranchStatus("*", 0);
        for (auto branch : this->branches_to_copy) {
            this->in_tree->SetBranchStatus(branch.c_str(), 1);
        }
        //for (auto branch : this->branches_to_rename) {
        //}
        for (auto branch : this->new_branches) {
        }
    }

    void setCut(string cut) {
        this->cut = cut;
    }

    void write() {
        cout << "Processing" << endl;

        this->setAllBranches();

        this->in_tree->Draw(">>event_list", this->cut.c_str(), "goff");
        TEventList *event_list = (TEventList*)gDirectory->Get("event_list");
        cout << "N events selected      : " << event_list->GetN() << endl;
        cout << endl;

        for (Long64_t i=0; i<event_list->GetN(); i++) {
            this->in_tree->GetEntry(event_list->GetEntry(i));
            this->out_tree->Fill();
        }
        this->out_tree->Write();
        this->out_file->Write();
    }

    void close() {
        in_file->Close();
        out_file->Close();
    }
};

//------------------
// HELPER FUNCTIONS
//------------------

struct Options {
    string in_file_name = "/public/data/Photon/Ntuples/bkg_data/data15-16_bkg.root";
    string in_tree_name = "BaselineTree";
    string out_file_name = "/public/data/Photon/SkimmedSamples/data15-16_bkg_compare.root";
    string out_tree_name = "BaselineTree";

    vector<string> branches_to_copy {
        "channel",
        "met_Et",
    };
    vector<pair<string, string>> branches_to_rename {
        make_pair("lepPt", "lep_pT"),
    };
    vector<string> new_branches {
        "empty_flag",
    };

    string cut = "met_Et>300";

    bool unit_testing = true;
};

Options setUnitTestOptions(Options options);
void performUnitTests(TTree* out_tree);

void makeReduceNtuples(Options options) {
    TreeReducer *reducer = new TreeReducer();

    if (options.unit_testing) {
        cout << "Performing unit testing" << endl;
        cout << endl;
        options = setUnitTestOptions(options);
    }

    reducer->openReadFileGetTree(options.in_file_name, options.in_tree_name);
    reducer->openWriteFileSetTree(options.out_file_name, options.out_tree_name);

    reducer->setBranchesToCopy(options.branches_to_copy);
    reducer->setBranchesToRename(options.branches_to_rename);
    reducer->setNewBranches(options.new_branches);

    reducer->setCut(options.cut);

    reducer->write();
    if (options.unit_testing)
        performUnitTests(reducer->out_tree);
    reducer->close();
}

//------------
// UNIT TESTS
//------------

Options setUnitTestOptions(Options options) {
    options.in_file_name = "/public/data/Photon/UnitTest/data15-16_bkg.root";
    options.out_file_name = "unit_test.root";
    options.branches_to_copy = vector<string> {
        "channel",
        "met_Et",
    };
    options.branches_to_rename = vector<pair<string, string>> {
        make_pair("lepPt", "lep_pT"),
    };
    options.new_branches = vector<string> {
        "empty_flag",
    };
    options.cut = "met_Et>300";

    return options;
}

void throwError(string error) {
    cout << "ERROR: " << error << endl;
    exit(0);
}

void performUnitTests(TTree* out_tree) {
    if (out_tree->GetEntries() == 394)
        cout << "Correct number of events after cut" << endl;
    else
        throwError("Wrong number of events after cut");

    if (out_tree->GetMinimum("met_Et") > 300)
        cout << "Cut performed correcly" << endl;
    else
        throwError("Cut performed incorrectly");

    cout << "Passed all unit tests" << endl;
    cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void ReduceNtuples() {
    Options options;

    options.in_file_name = "/public/data/Photon/Ntuples/bkg_data/data15-16_bkg.root";
    options.in_tree_name = "BaselineTree";
    options.out_file_name = "/public/data/Photon/SkimmedSamples/data15-16_bkg_compare.root";
    options.out_tree_name = "BaselineTree";

    options.branches_to_copy = vector<string> {
        "channel",
        "met_Et",
    };
    options.branches_to_rename = vector<pair<string, string>> {
        make_pair("lepPt", "lep_pT"),
    };
    options.new_branches = vector<string> {
        "empty_flag",
    };

    options.cut = "met_Et>300";

    options.unit_testing = true; makeReduceNtuples(options);
    options.unit_testing = false; makeReduceNtuples(options);
}
