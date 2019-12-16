#include "CommonLibraries.C"

namespace bins {

    //--- binning for smearing
    const int n_smearing_bins = 1000;
    const double smearing_low = -2000;
    const double smearing_high = 2000;

    const int n_pt_bins = 16;
    const double pt_bins[n_pt_bins+1] = {50,75,100,125,150,175,200,250,300,400,500,700,1000,1200,1400,1600,1e10};
    const double MET_bins[n_pt_bins+1] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,1e10,1e10};
    const double dphi_bin[n_pt_bins+1] = {0,0.5,1.0,1.5,2.0,2.5,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

    const int n_METl_bins = 25;
    const double METl_bins[n_METl_bins+1] = {-1e10,-1000,-700,-500,-400,-300,-250,-200,-150,-100,-60,-40,-20,20,40,60,100,150,200,250,300,400,500,700,1000,1e10};

    const int n_mll_bins = 43;
    const double mll_bin[n_mll_bins+1] = {12,20,30,40,50,60,70,80,82,84,86,88,90,92,94,96,98,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,440,480,520,560,600,800};

    TH1D *hist_METl_bins, *hist_pt_bins, *hist_MET_bins;
    void init_binning_histograms() {
        hist_METl_bins = new TH1D("hist_METl_bins","",n_METl_bins,METl_bins);
        hist_pt_bins = new TH1D("hist_pt_bins","",n_pt_bins,pt_bins);
        hist_MET_bins = new TH1D("hist_MET_bins","",n_pt_bins,MET_bins); //hist_MET_bins->SetStats(0);
    }

    //--- binning for reweighting and plotting (Ptll)
    const unsigned int n_reweighting_bins = 25;
    const double reweighting_bins[n_reweighting_bins+1] = {0,20,25,30,35,40,45,50,55,60,70,80,100,120,140,160,180,200,220,260,280,300,350,400,600,1000};
}
