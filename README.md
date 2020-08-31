# sigDriver
sigDriver is a driver discovery tool base on mutational signature exposures. Currently,  single base substitution signatures are supported.

# Prerequisites
sigDriver works in R (>= 3.4.0) and depends on the following packages: data.table, optparse, dplyr, GenomicRanges, reshape2, SKAT, matrixStats, rtracklayer, trackViewer, qqman, VariantAnnotation

# Installation
Installation in R using devtools:
```R
library(devtools)
install_github("wkljohn/sigDriver",subdir="R_library")
```

Installing the helper Rscript and the example data: 
```console
git clone https://github.com/wkljohn/sigDriver_runner.git
```

# How to use
**Please find the sigDriver helper scripts and example data in this repository:
[sigDriver_runner](https://github.com/wkljohn/sigDriver_runner)**

Run sigDriver using helper Rscripts
```console
Rscript run_sigdriver.R -s SBS9 -t 2 -v ./example/example_SBS9.snv.simple.gz -e ./example/example_signatures.tsv -m ./example/example_SBS9_metadata.tsv -o ./output/
```

Run results annotation using helper Rscripts
```console
Rscript run_sigdriver_annotate.R -g gencode.v19.gtf -t 2 -l ./output/SBS9_results.tsv -s SBS9 -v ./example/example_SBS9.snv.simple.gz -e ./example/example_signatures.tsv -m ./example/example_SBS9_metadata.tsv -o ./output/
```

# Tutorial
...Coming soon

# Contact
For questions and suggestions, please contact:

John Wong: wkljohn@gmail.com

