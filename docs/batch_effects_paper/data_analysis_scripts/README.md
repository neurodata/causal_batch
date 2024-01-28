This folder contains two drivers:
- [pairwise_analysis.R](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/data_analysis_scripts/pairwise_analysis.R) which assesses the detectability (or lackthereof) of batch effects for all pairs of datasets in the CoRR study, and
- [subseq_inference.R](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper/data_analysis_scripts/subseq_inference.R) which performs batch effect correction using various techniques on all of the CoRR studies followed by an overlapping subset, and then assesses the presence (or lackthereof) of sex effects at each edge in the corrected connectomes, conditional on age.


These drivers were executed on a high-performance server, and each took approximately 1 week to run on a machine with 100 cores and 1 TB RAM. The necessary pre-processed input data (fMRI connectomes) for these scripts can be downloaded from [the neurodata website](https://neurodata.io/mri/), from the subheading "functional MRI" for each dataset in the CoRR mega-study.
