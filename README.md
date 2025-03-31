# Modeling Glenohumeral Stability in Musculoskeletal Simulations  

This repository includes a customized version of the **Rapid Muscle Redundancy (RMR) solver** integrated with different formulations of **glenohumeral stability**.  

## Overview  

The original RMR solver, modeled glenohumeral stability as an **inequality constraint** that keeps the **glenohumeral contact force direction** within the stability border of the glenoid cavity, approximated as a **circle**.  

The original RMR solver is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8360269.svg)](https://doi.org/10.5281/zenodo.8360269)

In this customized version, we integrated different formulations to model glenohumeral stability using:  

- **Inequality constraints**, where the direction of the glenohumeral contact force was constrained within different stability border approximations:  
  - A **point**  
  - An **ellipse**  
  - A **polynomial fit** of empirical values obtained from **concavity compression tests**  
- **Penalty terms** in the objective function,
  - **Conditional penalty** penalty: trigers if the direction of the contact force exceedes a speficied circular stability border
  - **Planar penalty**: increases as the direction of the glenhumeral contact force points further away from the glenoid cavity border. Here, the glenoid cavity is considered as a flat plane.This plane is the plane that contains the glenoid caity border.
  -  **Curve penalty**: increases as the direction of the glenhumeral contact force points further away from the glenoid cavity border. Here, the curvature of the glenoid cavity is approximated as spherical.

 ## Musculoskeletal Model
 All the analysis was performed on the thoracoscapular shoulder model **DOI:** [10.3389/fnbot.2019.0009](https://doi.org/10.3389/fnbot.2019.00090)

## Validation  

Dataset of upper body **kinematics and glenohumeral contact force measurements** from an **instrumented shoulder prosthesis** were used to validate the estimated **magnitude** and **direction** of the glenohumeral contact force.  
Dataset available at **DOI:** [10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38](https://doi.org/10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38)
