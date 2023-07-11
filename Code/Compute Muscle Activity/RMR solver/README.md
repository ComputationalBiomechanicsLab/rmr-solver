The script to analyse a user-chosen motion is the `main_analyse_dataset.m`.
It can be run from Matlab, and will prompt the user with a selection of which task should be analysed. Then, the RMR solver (implemented in `RMR_analysis.m`) is used to find the optimal muscle activations corresponding to that particular movement.

Other important scripts are:
- `RMR_analysis_withOpenSimJRF.m`, which implements the RMR solver but computes the joint reaction force at the glenoid leveraging a non-closed expression of the force itself (much slower);
- `analyse_perturbed_models`, which is very similar to `main_analyze_dataset` but implemented to perform run the RMR solver on the noise-diversified dataset described in the paper.

The other scripts are documented individually.