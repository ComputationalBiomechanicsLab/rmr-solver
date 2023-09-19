## EMG data
This folder is used to collect the EMG data corresponding to the various tasks that are analyzed in our paper.
The movements were performed by a single subject, and are:
- unloaded abduction (ABD01, ABD02, ABD03)
- loaded abduction (ABD21, ABD22, ABD23)
- unloaded forward flexion (FLX01, FLX02, FLX03)
- loaded forward flexion (FLX21, FLX22, FLX23)
- unloaded shrugging (SHRUG01, SHRUG02, SHRUG03)
- loaded shrugging (SHRUG21, SHRUG22, SHRUG23)

You can get the data and populate your local version of the repository from the dataset available at https://simtk.org/projects/thoracoscapular.

It could be annoying, but we do not provide the dataset directly on purpose.
Why? We want to avoid duplicating it over and over, so that we can save memory and energy in the servers and appropriately credit the researchers who collected the data in the first place!

If you use this data alone, please cite:

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

If you employ our results on the dataset (or refer to them), please also cite:
```bib
@article {Belli2023.07.11.548542,
	author = {Italo Belli and Sagar Joshi and J Micah Prendergast and Irene Beck and Cosimo Della Santina and Luka Peternel and Ajay Seth},
	title = {Does enforcing glenohumeral joint stability matter? A new rapid muscle redundancy solver highlights the importance of non-superficial shoulder muscles},
	elocation-id = {2023.07.11.548542},
	year = {2023},
	doi = {10.1101/2023.07.11.548542},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/07/12/2023.07.11.548542},
	eprint = {https://www.biorxiv.org/content/early/2023/07/12/2023.07.11.548542.full.pdf},
	journal = {bioRxiv}
}
```