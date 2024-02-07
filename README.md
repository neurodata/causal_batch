# Causal Effect Detection and Correction

[![arXiv shield](https://img.shields.io/badge/arXiv-2307.13868-red.svg?style=flat)](https://arxiv.org/abs/2307.13868)
[![biorxiv shield](https://img.shields.io/badge/bioRxiv-2021.09.03.458920-blue.svg?style=flat)](https://www.biorxiv.org/content/10.1101/2021.09.03.458920v6)
[![](https://cranlogs.r-pkg.org/badges/causalBatch)](https://cran.rstudio.com/web/packages/causalBatch/index.html)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results and figure reproduction](#results-and-figure-reproduction)
- [License](./LICENSE)
- [Issues](https://github.com/ebridge2/causal_batch/issues)
- [Citation](#citation)


# Overview

Batch effects, undesirable sources of variance across multiple experiments, present significant challenges for scientific and clinical discoveries. Specifically, batch effects can (i) produce spurious signals and/or (ii) obscure genuine signals, contributing to the ongoing reproducibility crisis. Typically, batch effects are modeled as classical, rather than causal, statistical effects. This model choice renders the methods unable to differentiate between biological or experimental sources of variability, leading to unnecessary false positive and negative effect detections and over-confidence. We formalize batch effects as causal effects to address these concerns, and augment existing batch effect detection and correction approaches with causal machinery. Simulations illustrate that our causal approaches mitigate spurious findings and reveal otherwise obscured signals as compared to non-causal approaches. Applying our causal methods to a large neuroimaging mega-study reveals instances where prior art confidently asserts that the data do not support the presence of batch effects when we expect to detect them. On the other hand, our causal methods correctly discern that there exists irreducible confounding in the data, so it is unclear whether differences are due to batches or not. This work therefore provides a framework for understanding the potential capabilities and limitations of analysis of multi-site data using causal machinery.

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): usage of the `causalBatch` package on many real and simulated data examples for scientific articles.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.


# System Requirements

## Hardware Requirements

The `causalBatch` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3 GHz) and internet of speed 100 Mbps.

## Software Requirements

### OS Requirements

The package development version is tested on *Mac* operating systems. The developmental version of the package has been tested on the following systems:

Linux: 
Mac OSX:  Ventura 13.1
Windows:  

Before setting up the `causalBatch` package, users should have `R` version 4.2.0 or higher, and several packages set up from CRAN.

# Installation Guide

## Stable Release

The stable release of the package is available on CRAN, and can be installed from `R` as:

```
install.packages('causalBatch')
```

## Development Version

### Package dependencies

Users should install the following packages prior to installing `lolR`, from an `R` terminal:

```
install.packages(c('cdcsis', 'MatchIt', 'nnet', 'dplyr', 'tidyverse', 'magrittr'))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("sva")
```

which will install in about 1 minute on a machine with the recommended specs.

The `causalBatch` package functions with all packages in their latest versions as they appear on `CRAN` on January 22, 2024. The versions of software are, specifically:

```
sva=3.50.0
cdcsis=2.0.3
tidyverse=2.0.0
dplyr=1.1.4
MatchIt=4.5.5
nnet=7.3.19
magrittr=2.0.3
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/causal_batch/issues). 

### Package Installation

From an `R` session, type:

```
require(devtools)

# install causalBatch with the vignettes
install_github('neurodata/causal_batch', build_vignettes=TRUE, force=TRUE)

require(causalBatch)
# view one of the basic vignettes
vignette("causal_simulations", package="causalBatch") 
```

The package should take approximately 40 seconds to install with vignettes on a recommended computer. 

# Demo

For interactive demos of the functions, please check out the vignettes built into the package. They can be accessed as follows:

```
require(causalBatch)
vignette("causal_simulations", package="causalBatch")
vignette("causal_balancing", package="causalBatch")
vignette("causal_cdcorr", package="causalBatch")
vignette("causal_ccombat", package="causalBatch")
```

# Results and figure reproduction

See [Batch Effects Paper](https://github.com/neurodata/causal_batch/tree/main/docs/batch_effects_paper) for instructions to reproduce figures from Bridgeford et al. (2024). 

# Citation

For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](./citation.bib).
