# Rapid Muscle Redundancy Solver

## What is it?

**RMR solver** stands for Rapid Muscle Redundancy solver. It is an algorithm which solves the muscle redundancy problem, by selecting the muscles which are recruited by the human body to generate a given motion (leveraging a musculoskeletal model). The solver is presented in detail in our paper ("_Effects of glenohumeral joint stability on estimated shoulder muscle actions: a computational biomechanics study using a rapid muscle redundancy solver_"), and has been tested on the dataset and model available freely at https://simtk.org/projects/thoracoscapular (Seth et al, 2019).

## Features
The RMR solver solves for muscle and joint forces by implementing three new features: 

- Muscle forces include active and passive forces that are length and velocity dependent
- It enforces feasibility of the activation dynamics through linear constraints
- It formulates joint reaction forces as linear functions of the estimated muscle activation.

These features enable the efficient inclusion of active and passive muscle properties in the evaluation of muscle forces and joint stability during the optimization. 

In the paper, we considered specifically shoulder movements and the effect that the inclusion of the glenohumeral (GH) stability condition has on the muscle activations. Using filtered electromyography (EMG) data as ground truth for the muscle activations, the RMR solver proved to achieve the same accuracy as the widely employed computed muscle control (CMC) algorithm, but in a fraction of the time (below, the comparison for the analysis of motion data at 100 Hz).

<p align="center" width="100%">
    <img src='https://user-images.githubusercontent.com/50029203/221613483-3b5178cd-b777-4f54-a2a7-b8ce8f3e16e2.png' width=600>
</p>

## Requirements
In order to use the RMR solver, you will need:

- Software:

    -  **MATLAB 2021b**: to run the code in this repository
    - **OpenSim 4.3**: [Download here](https://simtk.org/frs/?group_id=91). The OpenSim team have published a setup guide [here](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab).
    - **SPM1D 0.4**: Get the MATLAB version for free at https://spm1d.org/ (used for statistical analysis, not within the solver)

    The code in this repository has been written and tested against those versions of the software. Older/different versions might not be supported.

- Data:

    - **EMG and marker data**, available for free (licensed CC BY 4.0) at: https://simtk.org/projects/thoracoscapular# . Following the link, you will be redirected to the SimTK page of the paper where the data is published. Clicking on `Download Latest Releases` and providing a simple motivation for why you are interested in the data, you will be able to download the original study, that contains the dataset in `ThoracoscapularShoulderPaperMaterials\ExperimentalData`. Just copy-paste the content of `EMG` and `Markers` in the same sub-folders of our `ExperimentalData` folder and you are ready to go! 
## Structure
The material is organized as follows:
- `Code` contains the scripts used for general data (pre)processing (`Data Processing`) and for running the RMR solver simulations as well as CMC (`Compute Muscle Activity`)
- `OpenSim Models` stores the biomechanical models that are considered;
- `ExperimentalData` stores the dataset considered in the paper (marker data and filtered EMG values)
- `Results` stores our results, supporting the reproducibility of our findings.

Each of the folders also have a more specific `README`.

## Brief guide to our code

1. In order to reproduce the figures and results commented in the paper, please run:
    - `Code\Image Generation\main_plotMuscleResults.m` (returns figures comparing the RMR solver against the CMC algorithm and filtered EMG data, and mean absolute errors resulting from this comparison);
    - `Code\Data Processing\spm_analysis.m` (performs a statistical analysis to determine the effect of the glenohumeral constraint on the muscle activations predicted by the RMR solver) 
    
      [Note: it does not run by itself, needs to be copy-pasted in the installation folder of the SPM1D package, as detailed in the script itself]

2. In order to use the RMR solver on the dataset, and to estimate realistic muscle activations given user-chosen motions, please run:
    - `Code\Compute Muscle Activation\RMR solver\main_analyze_dataset.m`


3. To perform a similar CMC analysis on the dataset, to reach our same results, please run:
    - `Code\Compute Muscle Activation\CMC analysis\batch_CMC_analysis.m` (by default considers the whole dataset)


If you wish to make use of the RMR solver, you can either download the individual files the you need or (recommended) you can fork this repository into your own, and work freely there! We encourage you to try the method on your own data/model: in case you do, remember to cite our publication!


## License
Our code is licensed under the Apache 2.0 license (see the `LICENSE_code` file), while the data and models are licensed under CC BY 4.0 Use Agreement terms.
```
Technische Universiteit Delft hereby disclaims all copyright interest in the program “RMR solver” developed to solve the muscle redundancy problem in biomechanical models written by the Author(s).
Prof. Dr. Ir. Fred van Keulen, Dean of Faculty of Mechanical, Maritime and Materials Engineering (3mE).
```
## Contributors
Italo Belli, Sagar Joshi, Irene Beck
