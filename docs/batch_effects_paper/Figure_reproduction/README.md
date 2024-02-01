The enclosed markdown files will produce first-pass figures for the manuscript:
- [simulations.Rmd](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/Figure_reproduction/simulations.Rmd) will reproduce Figures 1 and 4 up to machine-specific randomness,
- [demographics.Rmd](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/Figure_reproduction/demographics.Rmd) will reproduce Figure 5,
- [pairwise.Rmd](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/Figure_reproduction/pairwise.Rmd) will reproduce Figure 6, and
- [subseq_inference.Rmd](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/Figure_reproduction/subseq.Rmd) will reproduce Figure 7.

The figures found in the manuscript were further post-processed via Adobe Illustrator to fine-tune headings, text sizing, titles, and other non-data-related attributes for the purposes of making them more readable and digestable in the manuscript. Figures 2 and 3 were produced via Adobe Illustrator and are for illustrative use only.

To execute the figure reproduction code locally, we would recommend installation of the following packages. The version in which the figures were originally produced is also provided:

```
grid=4.3.2
gridExtra=2.3
ggpubr=0.6.0
reticulate=1.34.0
MatrixStats=1.2.0
BiocParallel=1.36.0
scales=1.3.0
cowplot=1.1.3
tidyvers=2.0.0
```

