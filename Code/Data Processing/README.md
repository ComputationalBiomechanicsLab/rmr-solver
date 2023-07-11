This folder collects scripts used for data processing. Many of them are general purpose scripts, whose use is documented in other functions or in the scripts themselves.

In particular:
- `SPM_robustness_analysis.m` allows to reproduce the robustness analysis performed in our paper. It requires 4 steps before running it:
		1. installing the SP1D package for Matlab (instructions at https://spm1d.org/install/InstallationMatlab.html)
		2. generating the perturbed models (using the scripts in `Code/Model Perturbation`)
		3. getting the marker data from the open-access study from Seth et al (available at: https://simtk.org/projects/thoracoscapular). The marker data should be copied into the `ExperimentalData/Markers` folder of this repository
		4. running the RMR solver on the models and markers retrieved above. The analysis can be performed through `Code/Compute Muscle Activity/RMR solver/analyse_perturbed_models.m`. Results should be stored in `Results/RMR analysis/robustness analysis`, following the indications provided there.
