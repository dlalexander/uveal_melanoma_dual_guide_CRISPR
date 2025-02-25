# Installation

Details of software dependencies, installation instructions and system requirements are specified below.

## Software dependencies

### System requirements

Analysis stages described in [README.md](https://github.com/dlalexander/uveal_melanoma_dual_guide_CRISPR/blob/master/README.md) were run on a high-performance computing (HPC) cluster running Ubuntu 20.04 (Focal Fossa) Linux using Nextflow version 21.04.3. Run time for this stage can be completed in half a day.
[Downstream analyses and plotting](https://github.com/dlalexander/uveal_melanoma_dual_guide_CRISPR/tree/master/SCRIPTS/plotting) can be performed on a local R installation, with library dependencies specified below. Run time for this stage can be completed in 1 hour. 

### C-SAR

Single guide analysis was performed using [C-SAR](https://github.com/cancerit/C-SAR) version 1.3.6, tag release [here](https://github.com/cancerit/C-SAR/releases/tag/1.3.6) including installation instructions.

C-SAR uses the packages:
* [MAGeCK version 0.5.9.3](https://sourceforge.net/projects/mageck/files/0.5/)
* [BAGEL2 version 2.0 build 114 (commit reference f9eedca)](https://github.com/hart-lab/bagel/tree/f9eedca7dc16299943dd1fd499bc1df4350ce8ef)

### R

Analyses were run using [R version 4.2.3](https://cran.r-project.org/) with the following package dependencies:

* data.table
* GGally
* ggplot2
* ggrepel
* gplots
* RColorBrewer
* reshape2
* stringi
* tidyverse
* vroom
* optparse

Details of the R dependencies for this project are detailed within the project `renv.lock` file.
