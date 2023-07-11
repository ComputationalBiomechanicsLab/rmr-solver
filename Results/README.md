In this folder we store the results of our study. 

- in `CMC analysis` you can find the results of analysing the various shoulder movements with the computed muscle control algorithm
- in `RMR analysis`, we store the results of the same analysis, performed with our solver on the same marker data (with and without including the glenohumeral constraint)
- in `IK solutions`, we share the results of our inverse kinematics analysis (which are the same both for the CMC and RMR tools)


Note that in `RMR analysis/robustness analysis` we provide the structure that you should populate yourself with the results that the RMR solver delivers when considering different perturbed models (see our paper and the content of `OpenSim Models/for RMR solver/perturbed models`). We do not share our results here as they consist of over 3000 files, that would make this repository very heavy to clone and store. However, we hope to provide sufficient instructions so that our results can be reproduced fully.
