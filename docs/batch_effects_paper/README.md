This folder contains the necessary figures and scripts to reproduce the results obtained in Bridgeford et al. (2024). This folder contains the following contents:

- [The drivers](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/data_analysis_scripts), which are the scripts which took connectomes from the CoRR mega-study and produced first-pass processed results (e.g., executed the statistical tests or requisite estimation procedures). 
- [The processed data](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/data), which are the processed derivatives produced by the drivers.
- [Utility functions](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/utilities) which are used by numerous methods in the drivers and figure reproduction.
- [Figure reproduction code](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/Figure_reproduction), which takes the processed data and produces the figures found in the main content of the manuscript. The figures found in the manuscript were further processed via Adobe Illustrator to fine-tune headings, text sizing, titles, and other non-data-related attributes for the purposes of making them more readable and digestable in the manuscript.

To execute the drivers, we recommend use of our [docker container](https://hub.docker.com/r/neurodata/causal_batch), which was used on our local servers for analyzing the pre-processed connectomes. To execute the figure reproduction code locally, we would recommend installation of the following packages. The version in which the figures were originally produced is also provided:

```
grid=4.3.2
gridExtra=2.3
ggpubr=0.6.0
reticulate=1.34.0
MatrixStats=1.2.0
BiocParallel=1.36.0
```
