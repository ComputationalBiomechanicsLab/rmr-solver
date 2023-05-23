## Marker data
In this folder, we collect the marker data corresponding to the various tasks that are analyzed in our paper.
The movements were performed by a single subject, and are:
- unloaded abduction (ABD01, ABD02, ABD03)
- loaded abduction (ABD21, ABD22, ABD23)
- unloaded forward flexion (FLX01, FLX02, FLX03)
- loaded forward flexion (FLX21, FLX22, FLX23)
- unloaded shrugging (SHRUG01, SHRUG02, SHRUG03)
- loaded shrugging (SHRUG21, SHRUG22, SHRUG23)

This data is part of the dataset available at https://simtk.org/projects/thoracoscapular, and we duplicate it here to facilitate reproducing of our results.

In the folder `noise_diversified` you can find the augmented version of this dataset, where Gaussian noise has been injected on the elbow markers (and only the relevant markers have been preserved). These files, grouped in 3 tasks (where the difference between loaded and unloaded conditions ceases to be relevant), are used in our paper to analyze statistically the effect of including the GH constraint in the estimation of the corresponding muscle activations.