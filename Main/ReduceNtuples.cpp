#include "Settings.cpp"

using namespace std;

//------------
// UNIT TESTS
//------------

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
        passTest("Correct number of events after cut");
    else
        failTest("Wrong number of events after cut");

    if (out_tree->GetMinimum("met_Et") > 300)
        passTest("Skimming performed correctly");
    else
        failTest("Skimming performed incorrectly");

    if (out_tree->GetNbranches() == 3+1+2)
        passTest("Slimming performed correctly");
    else
        failTest("Slimming performed incorrectly");

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
        passTest("Branch renaming performed correctly");
    else
        failTest("Branch renaming performed incorrectly");

    // need to add unit tests for new branches

    passTest("Passed all unit tests");
    cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void ReduceNtuples(ReductionOptions options) {
    TreeCreator *reducer = new TreeCreator();

    if (options.unit_testing) {
        cout << BOLD(PBLU("Performing unit testing on reduction step")) << endl;
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

    cout << PBLU("Done with reduction") << endl;
    cout << endl;
}
