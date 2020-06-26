## Overview

This code performs Z+jets background estimation for high-energy physics analysis, using alterations on data-driven photon+jets events.

Since Z+jets and photon+jets both involve single well-measured bosons recoiling against one or more jets, their kinematics and sources of mismeasurement are very similar. Thus, we can use photon+jets events to estimate Z+jets backgrounds in analyses, instead of relying on Monte Carlo methods.

## Main

The Main/ folder contains the bulk of the method.

This code should be run in the lxplus environment, for access to the relevant physics libraries. Simply edit options in Main.cpp and run via ROOT with the command "root -l -b -q Main.cpp".

The code can be run in four discrete steps, with intermediate data saved at the end of each step in specified folders. The steps are reduction, smearing, reweighting, and plotting.

The reduction step simply involves taking DAOD or AOD samples, extracting and renaming the relevant branches, and performing an event filter.

Smearing takes each photon event and adds uncertainties to MET values based on differences between photon and Z distributions. This step also artificially creates two leptons for each photon event via a stochastic process by boosting into the photon frame, calculating a two-particle decay in accordance with Z decay distributions, and boosting the daughter particles back into the lab frame.

Reweighting adds a per-event weight to each photon event in order to match kinematic distributions for any number of given features. For example, if Ptll and Ht30 are used as the reweighting parameters, then a branch will be created such that the photon and Z Ptll and Ht30 distributions match under that reweighting.

Plotting creates the final paper plots and tables.

## Misc

This folder contains various scripts for performing tests, manipulating samples, producing skims (for speeding up runtimes), and producing study plots.
