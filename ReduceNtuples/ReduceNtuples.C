#include "../Common/Settings.C"
#include <any>

using namespace std;
using BranchRenameOptions = vector<tuple<string, string>>;

class TreeReducer {
public:
    TFile *in_file, *out_file;    
    TTree *in_tree, *out_tree;    
    string out_file_name;
    string out_tree_name;
    string cut;
    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;
    unordered_map<string, std::any> new_branch_vars; 

    TreeReducer() {
    }

    void openReadFileGetTree(string file_name, string tree_name) {
        cout << "Opening read file      : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;

        this->in_file = TFile::Open(file_name.c_str());
        this->in_tree = (TTree*)this->in_file->Get(tree_name.c_str());

        cout << "Events in tree         : " << this->in_tree->GetEntries() << endl;
        cout << endl;
    }

    void openWriteFileSetTree(string file_name, string tree_name) {
        this->out_file_name = file_name;
        this->out_tree_name = tree_name;

        cout << "Opening write file     : " << file_name << endl;
        cout << "Tree name              : " << tree_name << endl;
        cout << endl;

        this->out_file = TFile::Open(file_name.c_str(), "recreate");
    }

    void setBranchesToCopy(vector<string> branches_to_copy) {
        this->branches_to_copy = branches_to_copy;
    }

    void setBranchesToRename(BranchRenameOptions branches_to_rename) {
        this->branches_to_rename = branches_to_rename;
    }

    void copyBranches() {
        for (auto branch : this->branches_to_copy) {
            this->in_tree->SetBranchStatus(branch.c_str(), 1);
        }
    }

    void renameBranches() {
        for (auto branch : this->branches_to_rename) {
            string old_name = get<0>(branch);
            string new_name = get<1>(branch);
            this->in_tree->SetBranchStatus(old_name.c_str(), 1);
            this->in_tree->GetBranch(old_name.c_str())->SetTitle(new_name.c_str());
            this->in_tree->GetBranch(old_name.c_str())->SetName(new_name.c_str());
        }
    }

    void initAllBranches() {
        this->in_tree->SetBranchStatus("*", 0);
        this->copyBranches();
        this->renameBranches();
        this->out_tree = this->in_tree->CloneTree(0);
    }

    void setCut(string cut) {
        this->cut = cut;
    }

    void write() {
        cout << "Processing" << endl;

        this->initAllBranches();

        this->in_tree->Draw(">>event_list", this->cut.c_str(), "goff");
        TEventList *event_list = (TEventList*)gDirectory->Get("event_list");
        cout << "N events selected      : " << event_list->GetN() << endl;
        cout << endl;

        for (Long64_t i=0; i<event_list->GetN(); i++) {
            this->in_tree->GetEntry(event_list->GetEntry(i));
            for (auto var : this->new_branch_vars)
                var.second = -999;
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
    string in_file_name;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;

    string cut;

    bool unit_testing;
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

    reducer->setCut(options.cut);

    reducer->write();
    if (options.unit_testing) {
        performUnitTests(reducer->out_tree);
        remove(reducer->out_file_name.c_str());
    }
    reducer->close();
}

//------------
// UNIT TESTS
//------------

void throwError(string error) {
    cout << "ERROR: " << error << endl;
    exit(0);
}

Options setUnitTestOptions(Options options) {
    options.in_file_name = "/public/data/Photon/UnitTest/data15-16_bkg.root";
    options.out_file_name = "unit_test.root";
    options.in_tree_name = "BaselineTree";
    options.out_tree_name = "BaselineTree";

    options.branches_to_copy = vector<string> {
        "channel",
        "met_Et",
        "lep_phi",
    };
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("lepPt", "lep_pT"),
    };

    options.cut = "met_Et>300";

    return options;
}

void performUnitTests(TTree* out_tree) {
    if (out_tree->GetEntries() == 394)
        cout << "Correct number of events after cut" << endl;
    else
        throwError("Wrong number of events after cut");

    if (out_tree->GetMinimum("met_Et") > 300)
        cout << "Skimming performed correctly" << endl;
    else
        throwError("Skimming performed incorrectly");

    if (out_tree->GetNbranches() == 3+1)
        cout << "Slimming performed correctly" << endl;
    else
        throwError("Slimming performed incorrectly");

    bool rename_pass = false;
    for (auto branch : *out_tree->GetListOfBranches()) {
        string branch_name = branch->GetName();
        if (branch_name == "lepPt") {
            rename_pass = false;
            break;
        }
        else if (branch_name == "lep_pT") {
            rename_pass = true;
        }
    }
    if (rename_pass)
        cout << "Branch renaming performed correctly" << endl;
    else
        throwError("Branch renaming performed incorrectly");

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
    options.out_file_name = "/public/data/Photon/SkimmedSamples/data15-16_bkg.root";
    options.out_tree_name = "BaselineTree";

    options.branches_to_copy = vector<string> {
        "channel",
        "met_Et",
        "lep_phi",
    };
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("lepPt", "lep_pT"),
    };

    options.cut = "met_Et>300";

    options.unit_testing = true; makeReduceNtuples(options);
    options.unit_testing = false; makeReduceNtuples(options);
}
