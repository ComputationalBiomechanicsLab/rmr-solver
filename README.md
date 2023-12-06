[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8360269.svg)](https://doi.org/10.5281/zenodo.8360269)
# RMR Solver  [<img align="right" src="https://github.com/ComputationalBiomechanicsLab/rmr-solver/assets/50029203/efb31852-625c-4456-91c5-ed7fc3ac80b3" width=150>](https://zenodo.org/record/8359200)

## What is it?

**RMR solver** stands for _Rapid Muscle Redundancy_ solver. It is an algorithm which solves the muscle redundancy problem, by selecting the muscles which are recruited by the human body to generate a given motion (leveraging a musculoskeletal model). The solver is presented in detail in our paper:

```bib
@article{belli2023does,
  title={Does enforcing glenohumeral joint stability matter? A new rapid muscle redundancy solver highlights the importance of non-superficial shoulder muscles},
  author={Belli, Italo and Joshi, Sagar and Prendergast, J Micah and Beck, Irene and Della Santina, Cosimo and Peternel, Luka and Seth, Ajay},
  journal={Plos one},
  volume={18},
  number={11},
  pages={e0295003},
  year={2023},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

To refer specifically to the data and model we use (available freely at https://simtk.org/projects/thoracoscapular), please also cite:

```bib
@article{seth2019muscle,
  title={Muscle contributions to upper-extremity movement and work from a musculoskeletal model of the human shoulder},
  author={Seth, Ajay and Dong, Meilin and Matias, Ricardo and Delp, Scott},
  journal={Frontiers in neurorobotics},
  volume={13},
  pages={90},
  year={2019},
  publisher={Frontiers Media SA}
}
```

<table align="center">
  <tr>
    <td colspan="2" align="center">Funding Institutions</td>
  </tr>
  <tr>
    <td align="center">
      <a>
        <img src="https://user-images.githubusercontent.com/50029203/226883398-97b28065-e144-493b-8a6c-5cbbd9000411.png" alt="TUD logo" height="128">
        <br />
        <a href="https://www.tudelft.nl/3me/over/afdelingen/biomechanical-engineering">Biomechanical Engineering</a> and <br />
        <a href="https://www.tudelft.nl/3me/over/afdelingen/cognitive-robotics-cor">Cognitive Robotics</a> at TU Delft</p>
      </a>
    </td>
    <td align="center">
      <a href="https://chanzuckerberg.com/">
        <img src="https://user-images.githubusercontent.com/50029203/226883506-fbb59348-38a4-43f9-93c9-2c7b8ba63619.png" alt="CZI logo" width="128" height="128">
        <br />
        Chan Zuckerberg Initiative
      </a>
    </td>
  </tr>
</table>

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

    -  **MATLAB 2021b**: to run the code in this repository (Optimization and Signal Processing Toolboxes are needed)
    - **OpenSim 4.3**: [Download here](https://simtk.org/frs/?group_id=91). The OpenSim team have published a setup guide [here](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab) that can be followed to enable MATLAB scripting.
    - **SPM1D 0.4**: Get the MATLAB version for free at https://spm1d.org/ (used for statistical analysis, not within the solver)

    The code in this repository has been written and tested against those versions of the software. Older/different versions might not be supported.

- Data:

    - **EMG and marker data**, available for free (licensed CC BY 4.0) at: https://simtk.org/projects/thoracoscapular# . Following the link, you will be redirected to the SimTK page of the paper where the data is published. Clicking on `Download Latest Releases` and providing a simple motivation for why you are interested in the data, you will be able to download the original study, that contains the dataset in `ThoracoscapularShoulderPaperMaterials\ExperimentalData`. Just copy-paste the content of `EMG` and `Markers` in the same sub-folders of our `ExperimentalData` folder and you are ready to go! 
## Structure
The material is organized as follows:
- `Code` contains the scripts used for general data (pre)processing and statistical analysis (`Data Processing`) and for running the RMR solver simulations as well as CMC (`Compute Muscle Activity`)
- `OpenSim Models` stores the biomechanical models that are considered;
- `ExperimentalData` will store the dataset considered in the paper (marker data and filtered EMG values), freely obtainable as described above;
- `Results` stores our results, supporting the reproducibility of our findings.

Each of the folders also have a more specific `README`.

## Brief guide to our code

1. In order to reproduce the figures and results commented in the paper, please run:
    - `Code\Image Generation\PlotResults.m` (produces figures comparing the RMR solver against the CMC algorithm and filtered EMG data, and mean absolute errors resulting from this comparison);
    - `Code\Data Processing\SPM_robustness_analysis.m` (performs a statistical analysis to determine the effect of the glenohumeral constraint on the muscle activations predicted by the RMR solver) 
    
      [Note: in order to run, it requires installing the SPM1D package]

2. In order to use the RMR solver on the dataset, and to estimate realistic muscle activations given user-chosen motions, please run:
    - `Code\Compute Muscle Activation\RMR solver\main_analyze_dataset.m`
  This scripts produces our numerical results and can be used as a template to apply the RMR solver in other scenarios.


3. To perform a similar CMC analysis on the dataset, please run:
    - `Code\Compute Muscle Activation\CMC analysis\batch_CMC_analysis.m` (by default considers the whole dataset)


If you wish to make use of the RMR solver, you can either download the individual files the you need or (recommended) you can fork this repository into your own, and work freely there! We encourage you to try the method on your own data/model: in case you do, remember to cite our publication!

### Trouble-shooting
If you encounter any troubles or issues with running the code contained in this repository, feel free to open an issue and report it. We will do our best to help you!


## License
Our code is licensed under the Apache 2.0 license (see the `LICENSE_code` file), while the data and models are licensed under CC BY 4.0 Use Agreement terms.
```
Technische Universiteit Delft hereby disclaims all copyright interest in the program “RMR solver”
developed to solve the muscle redundancy problem in biomechanical models written by the Author(s).

Prof. Dr. Ir. Fred van Keulen, Dean of Faculty of Mechanical, Maritime and Materials Engineering (3mE).
```
## Contributors
Italo Belli, Sagar Joshi, Irene Beck

## Acknowledgements
[<img src="https://github.com/ComputationalBiomechanicsLab/rmr-solver/assets/50029203/4227ef6e-cced-4b44-a707-292aec97562c" width=130>](https://codecheck.org.uk) 

We would like to thank the CodeCheck initiative for making sure that our code is reproducible! Their enthusiasm was really contagious, learn more about their great work at https://codecheck.org.uk
