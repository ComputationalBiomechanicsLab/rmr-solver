The dataset that we consider here is available in its complete form at https://simtk.org/projects/thoracoscapular

To obtain the `EMG` and `Markers` data, just visit the link provided above, download the original material and simply-copy paste the content of the original `EMG` and `Markers` folders in the two corresponding folders here.
In particular, the data needed is:
- marker data 3D trajectories (in `.trc` format)
- filtered EMG values for superficial muscles (ground-truth for the muscle activations, in `.exp` format)

Inside the `Markers\noise_diversified` folder you can find already the dataset augmented with white Gaussian noise, described in our paper.

The folder `IK setup files` collects the setup files needed to perform the inverse kinematics on the basis of the experimental marker data, and the musculoskeletal model considered.