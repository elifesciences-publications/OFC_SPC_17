
### This code is associated with the paper from Sadacca et al., "Orbitofrontal neurons signal sensory associations underlying model-based inference in a sensory preconditioning task". eLife, 2017. http://dx.doi.org/10.7554/eLife.30373


# OFC_SPC_17
Scripts and functions for the analysis of orbitofrontal cortex single-unit data collected during sensory preconditioning

Processing begins from loading sorted single-unit data in bo_run_main script - which calls functions to process behavioral data, bin timestamped spiking data, exclude data on the basis of signal-to-noise ratio or single-unit drift. 

Other scripts respectively plot single-unit and average data (bo_plot_line.m), correlational analysis (bo_rho_resample.m), classification individual trials for preconditioning (bo_preconditioning_classification.m) or in the probe test (bo_probe_classification.m)
