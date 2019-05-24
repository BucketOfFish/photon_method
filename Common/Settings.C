#ifndef COMMON_SETTINGS
#define COMMON_SETTINGS

std::string sample_folder = "NonLepIso";
std::string ntuple_path = "/eos/user/m/mazhang/PhotonMethod/v1.6/" + sample_folder + "/Ntuples/";
std::string smearing_path = "/eos/user/m/mazhang/PhotonMethod/v1.6/" + sample_folder + "/SmearedNtuples/";
std::string reweighting_path = "/eos/user/m/mazhang/PhotonMethod/v1.6/" + sample_folder + "/ReweightedNtuples/";
std::string plots_path = "/eos/user/m/mazhang/PhotonMethod/v1.6/" + sample_folder + "/Plots/";

//--- ntuple production
double leading_lep_pt_cut = 25.; // also used for smearing
double second_lep_pt_cut = 25.; // also used for smearing

//--- binning for smearing methods
const int bin_size = 22;
//double sm_pt_bin[bin_size+1] ={0,50,75,100,125,150,175,200,250,300,400,500,700,1000,1200,1400,1600,1e10,1e10,1e10,1e10,1e10,1e10};
double sm_pt_bin[bin_size+1] ={50,75,100,125,150,175,200,250,300,400,500,700,1000,1200,1400,1600,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
double met_bin[bin_size+1] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
double dphi_bin[bin_size+1] = {0,0.5,1.0,1.5,2.0,2.5,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
const int dpt_bin_size = 25;
double dpt_bin[dpt_bin_size+1] = {-1e10,-1000,-700,-500,-400,-300,-250,-200,-150,-100,-60,-40,-20,20,40,60,100,150,200,250,300,400,500,700,1000,1e10};
const int mll_bin_size = 43;
double mll_bin[mll_bin_size+1] = {12,20,30,40,50,60,70,80,82,84,86,88,90,92,94,96,98,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,440,480,520,560,600,800};

//--- binning for reweighting and plotting
const unsigned int nptbins = 16;
double ptbins[nptbins+1] = {40, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 850, 1000};

#endif
