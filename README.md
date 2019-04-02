Authors: Ruo-Yu Shang, Johny Echevers, Matt Zhang

"source run_all.sh" to run everything.

1. First, we prepare the environment on lxplus.

2. Next, we produce skimmed ntuples for photon/Z data/MC.

3. Then we produce smeared photon data/MC ntuples.

    There are 4 different modes of smearing:
        smearing_method = 0; # no smearing
        smearing_method = 1; # use MC Sherpa Z + 1 jet events and photon + 1 jet events
        smearing_method = 2; # use data Z + 1 jet events and photon + 1 jet events
        smearing_method = 3; # use MC Sherpa Z + jets truth response function (not recommended)
    The smearing functions are sliced in boson pT, as defined in BasicSetting.C (sm_pt_bin).
    The histograms to build smearing functions are sampled here: GetSmearingHistogram.C.
    After we sample the histograms, we use deconvolution to derive smearing function in GetPhotonSmearing.C.
    The deconvolution is done in line 156-245: pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1).
    Once the smearing function is derived, we use it to smear photon resolution in line 472 (gamma_pt_smear = gamma_pt-photon_smear), and MET in line 477 (METl_smear = METl + photon_smear_l).
    Mll modeling is also done in GetPhotonSmearing.C in line 113: GetMllHistogram(ch).
    After Mll, we compute the two lepton kinematics in line 538: GetIndividualLeptonInfo(z_4vec).

4. After this we reweight photon data/MC.

    The reweighting binning is defined in BasicSetting.C (pt_bin).
    The reweighting histograms are sampled in GetReweightingHistogram.C.
    The reweighting histogram with the most inclusive selection is hist_bveto_*, which should be good enough for any possible analysis region.
    However, you can also define dedicated selection for your own analysis regions. For example, look for hist_ht200_*, which is the reweighting factor that we used for strong 2L SR with HT>200 GeV. 
    Keep in mind that you shouldn't have MET cuts or Delta Phi cuts in the selection of reweighting factors.

5. Now we examine results of smearing and reweighting.

    The photon MET distribution (red) is compared with the Z+jets MET distribution (blue).
    If the 1-jet region closure is bad in a pT range, it is often due to the poor statistics in that range.
    We need good statistics to make a shape of MET parallel for deconvolution.
    If a bin does not have enough stats, go to BasicSetting.C and change the binning of sm_pt_bin. 

6. Making plots in the Strong2L 2017 paper
    cd MakeHistogram
    root -l
    .L MakeSelectionHistogram_Strong2LPaper2017.C++
    MakeSelectionHistogram_Strong2LPaper2017()
    This command will produce many histograms in OutputHistogram/
    When the job is done, go to MakePlots/
    cd MakePlots
    python SimplePlot_Strong2LPaper2017.py
    Then you got plots for the 2017 paper!
