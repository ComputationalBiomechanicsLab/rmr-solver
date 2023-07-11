Store here the results of the RMR solver analysis performed on the perturbed models (see `OpenSim Models/perturbed models`)

In order for the SPM analysis script (`Code/Data Processing/SPM_robustness_analysis`) to work as-is, run the `analyse_perturbed_models.m` in `Code/Compute Muscle Activity/RMR solver`. Select all the models that you would like to analyse, and the corresponding `.trc` marker trajectories.

Save in here the results of the analysis performed without the GH constraint on the unloaded flexions (`flx01`, `flx02`, `flx03` - using the correct models with increased hand-mass).