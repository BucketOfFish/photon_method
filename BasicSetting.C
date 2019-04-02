//-----------------------------------------------------------------------------------------------
// this script defines the binning for photon smearing and reweighting.
// bool useDeconvolution allows you to switch between deconvolution method and the standard smearing method
//-----------------------------------------------------------------------------------------------

double lumi = 36100.;

std::string ntuple_path = "/eos/user/m/mazhang/PhotonMethod/v1.6/Ntuples/";
std::string smearingPath = "/eos/user/m/mazhang/PhotonMethod/v1.6/SmearedNtuples/";
std::string oldSmearingPath = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/PhotonMethodOutput/OutputNtuples_EWK2L/";
std::string plotsPath = "/eos/user/m/mazhang/PhotonMethod/v1.6/Plots/";

string photon_tag = "";
int event_interval = 1; // use every nth event
bool useMETtrig = false;

double leading_lep_pt_cut = 25.;
double second_lep_pt_cut = 25.;

const int bin_size = 22;
// binning for smearing methods
double sm_pt_bin[bin_size+1] ={50, 75,100,125,150 ,175 , 200,250 ,300, 400 ,500 ,700 ,1000,1200,1400,1600,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

// binning for reweighting methods
double njet_bin[bin_size+1] = {0, 1 ,2  ,3  ,4  ,5  ,6   ,7   , 8  , 9  , 10,   11,  12,  13,  14,  15,  16,  17,  18,  19,  20,1e10,1e10};
double pt_bin[bin_size+1] =   {50, 75,100,125,150 ,175 , 200,250 ,300,  350, 400, 450, 500, 600, 700, 850,1000,1200,1400,1600,1e10,1e10,1e10};
double et_bin[bin_size+1] =   {0, 20, 40, 60, 80,100,120 ,140 , 160,200 ,300 ,400 ,600 ,800 ,1000,1500,2000,1e10,1e10,1e10,1e10,1e10,1e10};
double lpt_bin[bin_size+1] =  {0, 20, 30, 40, 50, 60, 70 , 80 ,  90,100 ,120 ,140 , 160,200 ,300 ,400 ,600 ,800 ,1000,1500,1e10,1e10,1e10};
double ht_bin[bin_size+1] =   {0, 60, 80,100,120,140,160 ,180 ,200 ,250 ,300 ,350 ,400 ,500 ,600 ,800 ,1000,1200,1400,1800,1e10,1e10,1e10};

// binning for Jpsi/Upsi methods
double bphys_pt_bin[bin_size+1] ={0, 2, 4, 6, 8, 10, 15,  20, 25, 30, 35, 40, 45, 50,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
double bphys_sm_pt_bin[bin_size+1] ={0, 2, 4, 6, 8, 10, 15,  20, 25, 30, 35, 40, 45, 50,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

// not using them
double met_bin[bin_size+1] =  {0, 20,40 ,60 ,80 ,100,120 ,140 ,160 ,180 ,200,250  ,300 ,350 ,400 ,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
double dphi_bin[bin_size+1] = {0,0.5,1.0,1.5,2.0,2.5,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

// binning for mll assignments
const int dpt_bin_size = 25;
double dpt_bin[dpt_bin_size+1] = {-1e10,-1000,-700,-500,-400,-300,-250,-200,-150,-100,-60,-40, -20,20,40,60,100,150,200,250,300,400,500,700,1000,1e10};
const int mll_bin_size = 43;
double mll_bin[mll_bin_size+1] = {12,20,30,40,50,60,70,80,82,84,86,88,90,92,94,96,98,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,440,480,520,560,600,800};

