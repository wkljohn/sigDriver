# sigDriver
sigDriver is an R package for driver discovery using mutational signatures

# Installation
Installing sigDriver from shell
```console
git clone https://github.com/wkljohn/sigDriver.git`
R CMD INSTALL sigDriver-master`
```

Alternative installation in R using devtools:
```R
library(devtools)
install_github("wkljohn/sigDriver")
```

Installing the helper Rscript and the example data
```console
git clone https://github.com/wkljohn/sigDriver_runner.git
```

# How to use

Run sigDriver using helper Rscripts
```console
Rscript run_sigdriver.R -s SBS9 -v ./example/example_SBS9.snv.simple.gz -e ./example/example_signatures.tsv -m ./example/example_SBS9_metadata.tsv -o ./output/
```

Run annotation for sigDriver results using helper Rscripts
```console
Rscript run_sigdriver_annotate.R -g gencode.v19.gtf -l ./output/SBS13_results.tsv -s SBS9 -v ./example/example_SBS9.snv.simple.gz -e ./example/example_signatures.tsv -m ./example/example_SBS9_metadata.tsv -o ./output/
```

# Tutorial
...Coming soon