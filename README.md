# sigDriver


# Installation
Installing sigDriver from shell
```console
$ git clone https://github.com/wkljohn/sigDriver.git`
$ R CMD INSTALL sigDriver-master`
```

Alternative installatio in R using devtools:
```R
library(devtools)
install_github("wkljohn/sigDriver")
```

Installing the helper Rscript and example data
```console
$ git clone https://github.com/wkljohn/sigDriver_runner.git
```

# How to use

Run sigDriver using helper Rscripts
```console
$ Rscript run_sigdriver.R -s SBS13 -v ./example/XXXX.simple.gz -e ./example/XXXX_signatures.tsv -m ./example/XXXX_metadata.tsv -o ./output/
```

Run annotation for sigDriver results using helper Rscripts
```console
$ Rscript run_sigdriver_annotate.R -g gencode.v19.gtf -l ./output/SBS13_results.tsv -s SBS13 -v ./example/XXXX.simple.gz -e ./example/XXXX_signatures.tsv -m ./example/XXXX_metadata.tsv -o ./output/
```