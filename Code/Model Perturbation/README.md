In our paper we generated 200 different models to test for the robustness of our conclusions against intrinsic uncertainty of the study. We do not share the models, as this would clutter this repo, but we provide the scripts to generate them.

In order to achieve similar results to ours:
1. run the `perturb_models.m` to produce 100 different models, selecting as input the model in `OpenSim Models/for RMR solver/TSM_subject_noWeight.osim`. Copy-paste the resulting models into  `OpenSim Models/for RMR solver/perturbed models/noWeight`. 
3. Run the script `add_2kgLoad_toModels.m`, selecting the results of the previous step, and copy-paste the resulting models in `OpenSim Models/for RMR solver/perturbed models/2kgWeight`.