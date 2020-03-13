#include "Common/Settings.C"
#include "TInterpreter.h"

using namespace std;

class TreeReducer {
public:
    ROOT::RDataFrame *dataframe;
    string out_file_name;
    string out_tree_name;
    string cut;
    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;
    BranchAddOptions branches_to_add;

    TreeReducer() {
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

        reduced_dataframe.Snapshot(this->out_tree_name.c_str(), out_file_name.c_str(), all_out_branches);
        cout << endl;
    }
};

//------------------
// HELPER FUNCTIONS
//------------------

ReductionOptions setUnitTestOptions(ReductionOptions options);
void performUnitTests(TTree* out_tree);

void runReduction(ReductionOptions options) {
    TreeReducer *reducer = new TreeReducer();

    if (options.unit_testing) {
        cout << "Performing unit testing" << endl;
        cout << endl;
        options = setUnitTestOptions(options);
    }

    reducer->read(options.in_file_name, options.in_tree_name);

    reducer->setBranchesToCopy(options.branches_to_copy);
    reducer->setBranchesToRename(options.branches_to_rename);
    reducer->setBranchesToAdd(options.branches_to_add);

    reducer->setCut(options.cut);

    reducer->write(options.out_file_name, options.out_tree_name);

    if (options.unit_testing) {
        TFile *out_file = TFile::Open(options.out_file_name.c_str());
        TTree *out_tree = (TTree*)out_file->Get(options.out_tree_name.c_str());
        performUnitTests(out_tree);
        remove(options.out_file_name.c_str());
    }
}

//------------
// UNIT TESTS
//------------

void throwError(string error) {
    cout << "ERROR: " << error << endl;
    exit(0);
}

ReductionOptions setUnitTestOptions(ReductionOptions options) {
    options.in_file_name = "/public/data/Photon/UnitTest/data15-16_bkg.root";
    options.out_file_name = "unit_test.root";
    options.in_tree_name = "BaselineTree";
    options.out_tree_name = "BaselineTree";

    options.branches_to_copy = vector<string> {
        "channel",
        "met_Et",
        "lepFlavor",
    };
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("lepPt", "lep_pT"),
    };
    string getChannel =
        "int getChannel(const ROOT::VecOps::RVec<int> &lepFlavor) {"
            "if (lepFlavor[0] == 2 and lepFlavor[1] == 2) return 0;"
            "else if (lepFlavor[0] == 1 and lepFlavor[1] == 1) return  1;"
            "else if (lepFlavor[0] == 1 and lepFlavor[1] == 2) return 2;"
            "else if (lepFlavor[0] == 2 and lepFlavor[1] == 1) return 3;"
            "return -1;"
        "}";
    string countTo3 =
        "vector<int> countTo3() {"
            "return vector<int>{1, 2, 3};"
        "}";
    gInterpreter->Declare(getChannel.c_str());
    gInterpreter->Declare(countTo3.c_str());
    options.branches_to_add = BranchAddOptions {
        make_tuple("test1", "getChannel(lepFlavor)"),
        make_tuple("test2", "countTo3()"),
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

    if (out_tree->GetNbranches() == 3+1+2)
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

    // need to add unit tests for new branches

    cout << "Passed all unit tests" << endl;
    cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void ReduceNtuples(ReductionOptions options) {
    runReduction(options);
}
