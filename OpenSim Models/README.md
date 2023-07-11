The musculoskeletal models considered in this study are modifications of the one available at: https://simtk.org/projects/thoracoscapular#

The models in the two folders are different as the RMR solver does not need to lock any coordinate, and requires coordinate actuators to be present in the model (for each coordinate).
On the other hand, for CMC we lock the x-y-x degrees of freedom (complying to the study which originally analyzed the dataset we consider, "Muscle Contributions to Upper-Extremity Movement and Work From a Musculoskeletal Model of the Human Shoulder", by Seth A et al.).

Since the movements in the dataset are performed with and without a 2-kg load held by the subject, two models are present to represent the two conditions (which differ just for the mass and inertia of the `hand`).

Lastly, to correctly consider the glenohumeral constraint, we add three extra markers on the model considered by the RMR solver, to define:
- the center of the glenoid (`Glenoid_Center`);
- the center of the humeral head (`HumHead_Center`);
- the width of the glenoid (`Glenoid_Edge`).

When considering your own model with the RMR solver, these markers should be added there as well.