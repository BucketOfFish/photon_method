Setup:
	setupATLAS
	lsetup git

Check out the code:
	git clone https://:@gitlab.cern.ch:8443/rshang/PhotonTemplateMethod.git
	git clone https://github.com/ruoyushang/Photon_method

A tutorial of the photon method:


	1. Generate small scripts to produce skimmed ntuples for the photon method:
		cd ScriptMaker
		python makeRunScript.py # this will generate scripts (each contains a single-line command) that run small single jobs in Scripts/
		For example, take a look at Scripts/run_periodA_gdata16.sh:
			root -l -b -q 'GetPhotonEvents.C+("periodA","gdata","/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gdata16/",1)'
		This command takes input from /eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gdata16/, which stores the photon data files that we used for 2017 strong analysis.
		If you wish to change the input file directory, please modify makeRunScript.py

	2. Use the small scripts to produce skimmed ntuples
		Now we need to compile the photon code
		cd ../Scripts
		sh run_361419_zmm.sh # this is a small job that runs over a Sherpa Z+jets sample (DSID361419) and compiles Root/GetBaseLineEvents.C 
		sh run_361062_gmc.sh # this is a small job that runs over a Sherpa photon+jets sample (DSID361062) and compiles Root/GetPhotonEvents.C

		Root/GetBaseLineEvents.C and Root/GetPhotonEvents.C take inputs from QuickAna ntuples
		and check the trigger requirements and analysis selections.
		The input directory of this tutorial is /eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/MC/
		To take inputs from your customized ntuples, you will need to change the variable names in
		Root/InputVariables.C, Root/GetBaseLineEvents.C and Root/GetPhotonEvents.C

		After we compile the code, now we can run
		sh local_z.sh		# this creates Sherpa 2.1 Z+jets ntuples for MC-closure test
		sh local_gmc.sh		# this creates Sherpa 2.1 photon+jets ntuples for MC-closure test
		sh local_data.sh	# this creates dilepton data ntuples for analysis
		sh local_gdata.sh	# this creates photon+jets data ntuples for analysis prediction
		sh local_tt.sh		# this creates Top MC ntuples for analysis prediction
		sh local_vv.sh		# this creates diboson MC ntuples for analysis prediction
		sh local_Vg.sh		# this creates V+gamma MC ntuples for photon sample contamination prediction
		sh local_z221.sh	# this creates Sherpa 2.2.1 Z+jets ntuples for analysis prediction
		sh local_bphys.sh	# this creates B-physics data ntuples for low pT analysis prediction
		Then we will see many ntuples are produced in OutputNtuples/

	3. Merge the ntuples that you just produced
		We need to merge the ntuples of individual DSID into inclusive files:
		cd ../OutputNtuples
		sh hadd.sh
		this will create OutputNtuples/gmc/gmc_raw.root, OutputNtuples/zjets/zjets_ee.root and  OutputNtuples/zjets/zjets_mm.root, etc.
		After merging all ntuples, you can delete the old ntuples with individual DSID to save your storage space. 

	4. Run the photon smearing
		After merging the files, go back to Scripts/
		cd ../Scripts
		sh run_gmc_sm.sh	# this script will derive smearing function from 1-jet events in gmc_raw.root, zjets_ee.root, and zjets_mm.root
					# and smear the photon resolution 
					# and create OutputNtuples/gmc/gmc_ee_McSmear.root and OutputNtuples/gmc/gmc_mm_McSmear.root 
		sh run_gdata_sm.sh	# this script does the smearing for photon data
		There are 4 different modes of smearing in Root/BasicSetting.C
			int smearing_method = 0; # No smearing will be applied.
			int smearing_method = 1; # use the MC smearing funtion from Sherpa Z+ 1 jet events and photon + 1 jet events
			int smearing_method = 2; # use data smearing function from data Z + 1 jet events and photon + 1 jet events
			int smearing_method = 3; # Truth smearing using use the MC smearing funtion from Sherpa Z+jets truth response function. This option is not recommended.
		The smearing functions are sliced in boson pT, the slicing is defined in Root/BasicSetting.C:
			double sm_pt_bin[bin_size+1] ={50, 75,100,125,150 ,175 , 200,250 ,300, 400 ,500 ,700 ,1000,1200,1400,1600,1e10,1e10,1e10,1e10,1e10,1e10,1e10};
			You can change this binning for your analysis. Make sure you have enough statistics for each slice.
		The histograms to build smearing functions are sampled here: Root/GetSmearingHistogram.C
		After we sample the histograms, we use deconvolution to derive smearing function in Root/GetPhotonSmearing.C
		The deconvolution is done in line 156-245: pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1)
		Once the smearing function is derived, we use it to smear photon resolution in line 472 (gamma_pt_smear = gamma_pt-photon_smear), and MET in line 477 (METl_smear = METl + photon_smear_l).
		Mll modeling is also done in Root/GetPhotonSmearing.C in line 113: GetMllHistogram(ch)
		After Mll, we compute the two lepton kinematics in line 538: GetIndividualLeptonInfo(z_4vec)
			

	5. Run the photon reweighting
		sh run_gmc_rw.sh	# this script will derive reweighting factors from gmc_raw.root, zjets_ee.root, and zjets_mm.root
		sh run_gdata_rw.sh	# this script does the reweighting for photon data
		The reweighting binning is defined in Root/BasicSetting.C: pt_bin[bin_size+1]
		The reweighting histograms are sampled in Root/GetReweightingHistogram.C
		The reweighting histogram with the most inclusive selection is hist_bveto_*, which should be good enough for any possible analysis region.
		However, you can also define dedicated selection for your own analysis regions. For example, look for hist_ht200_*, which is the reweighting factor that we used for strong 2L SR with HT>200 GeV. 
		Keep in mind that you shouldn't have MET cuts or Delta Phi cuts in the selection of reweighting factors.

	6. A quick look at your photon method output
		Up to here, the photon ntuple is corrected with smearing and reweighting.
		To have a quick look at the output:
		cd ../OutputNtuples
		root -l quickDraw.C 
		and you will see the photon MET distribution (red) is compared with the Z+jets MET distribution (blue)

		in quickDraw.C, you can see:

		The selection in the plot:
			TCut mycut = "";
			mycut += "bjet_n==0";  # the photon prediction is not restricted to b-veto region, you can remove this cut and see the difference
			mycut += "lep_pT[0]>20 && lep_pT[1]>20";
			mycut += "abs(lep_eta[0])<2.5 && abs(lep_eta[1])<2.5";
			mycut += "jet_n>=2";
			mycut += "Z_pt>50"; 	# the phtoon method does not provide low pT estimation, 
						# low pT (<50GeV) region will be covered by Upsilon template method.

		and the reweighitng factor:
			TCut weight_g = "totalWeight*36.1*1000*ptrw_bveto";	# totalWeight includes xsec, eventWeight, trigger weight, etc.
										# ptrw_bveto is the reweighting factor derived from a b-veto region	

		sanity check the smearing function:
		root -l drawSmearingHistogram.C
		you will see a plot of MET parallel distributions (projection of MET in the Z or photon direction):
			-- blue: Z+jets
			-- green: photon+jets before smearing
			-- red: photon+jets after smearing
			-- pink: smearing function used to correct photon resolution
		these plots are made in 1-jet region (where smearing function is derived),
		so you should see a good agreement between blue and red curves by design.
		In drawSmearingHistogram.C, TString index = "10" gives you the distribution in 500 < Z pT < 700 GeV range.
		Change the number to see other pT ranges.
	

		If the 1-jet region closure is bad in a pT range, it is often due to the poor statistics in that range.
		We need good statistics to make a shape of MET parallel for deconvolution.
		If there is one range does not have enough stat, then you should go to Root/BasicSetting.C and change the binning of sm_pt_bin 
		to make sure you have enough statistics in the relavent pT range. 
	
		In Root/BasicSetting.C, you can also see the variable event_interval = 10, this means this tutorial only runs 10% events in the samples.
		To run full stat job, change event_interval to 1.

		In Root/BasicSetting.C, smearing_method = 1 means the deconvolution is done using MC 1-jets events.
		You can change smearing_method to 2 if you want to derive the smearing function from data 1-jet events, or to 0 if you don't want smearing. 

	7. Making plots in the Strong2L 2017 paper
		cd MakeHistogram
		root -l
		.L MakeSelectionHistogram_Strong2LPaper2017.C++
		MakeSelectionHistogram_Strong2LPaper2017()
		This command will produce many histograms in OutputHistogram/
		When the job is done, go to MakePlots/
		cd MakePlots
		python SimplePlot_Strong2LPaper2017.py
		Then you got plots for the 2017 paper!



If you have further questions, please send an email to ruo-yu.shang@cern.ch

