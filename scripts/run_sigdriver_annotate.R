#!/usr/bin/env Rscript

#Libraries

#functions:
#output most affected gene
#output sites used for test
#output importance of each site
#output plots?

#parameters parsing library
library(optparse)

#testing
if (F){
	#Rscript sigdriver_annotate.R /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv SBS2 /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/ /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/SBS2_results.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz
	#PCAWG  cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=6 -l mem=20g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -d /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1_E05/"$0"_var_meta.rds -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1_E05/ -l /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1_E05/"$0"_results.tsv -t 10 -g /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz\" run_sigdriver_annotate.R"}'
	#B:	 cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv | sed '1d' | awk '{print "bsub -m \"abi-cn24u12 abi-cn24u11 abi-cn24u13 abi-cn24u14 abi-cn24u15 abi-cn24u16 abi-cn24u17 abi-cn24u18 abi-cn24u19 abi-cn24u20 abi-cn24u21 abi-cn24u22 odcf-cn34u03s07 odcf-cn34u03s10 odcf-cn34u09s12 odcf-cn34u18s01 odcf-cn34u18s02 odcf-cn34u18s03 odcf-cn34u18s04 odcf-cn34u18s05 odcf-cn34u18s11 odcf-cn34u18s12 odcf-cn31u17\" -q verylong -n 6 -M 20G  -R \"span[hosts=1] rusage[mem=20GB]\" -o \""$0"\".log -e \""$0".err\"  \"module load r/3.6.0;module load gcc/7.2.0;Rscript sigdriver_annotate.R -v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/ -l /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/"$0"_results.tsv -t 10 -g /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz\""}'
	#PCAWGE cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=10 -l mem=20g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/ -l /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/"$0"_results.tsv -g /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz -t 10\" run_sigdriver_annotate.R"}'
	#CATCH	cut -f 1 /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=6 -l mem=10g -F \"/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv "$0" /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/ /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/"$0"_results.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz\" sigdriver_annotate.R"}'
	
	variant_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz"
	signature_file="/b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv"
	covar_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv"
	signature_test="SBS13"
	threads=10
	out_path="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/"
	results_file=paste("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/",signature_test,"_results.tsv",sep="")
	annotation_gtf="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.annotation.gtf.gz"
	source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association_release/sigdriver_annotate.R")

	out_path="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/"
	results_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/SBS9_results.tsv"
	annotation_gtf="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation.gtf"
	source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association_release/sigdriver_annotate.R")

	variant_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz"
	signature_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv"
	covar_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv"
	#signature_test="SBS2,SBS13"
	signature_test="SBS2"
	threads=10
	out_path="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/"
	results_file=paste("/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/",signature_test,"_results.tsv",sep="")
	#annotation_gtf="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation.gtf"
	annotation_gtf="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/gencode.v29lift37.basic.annotation.gtf.gz"
	source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association_release/sigdriver_annotate.R")

}


if (T){
	#Rscript run_sigdriver_annotate.R -l /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/SBS2_results.tsv -g /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation.gtf -s SBS13 -v /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz -t 8 -e /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv -o ./output/

	#Rscript run_sigdriver_annotate.R -l /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/SBS2_results.tsv -g /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation.gtf -s SBS13 -v /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz -t 8 -e /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv -o /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/test/
	help_message <-
"usage: Rscript sigdriver_annotate.R [options] <-l results file to load> <-g gtf file> <-s signature test> <-v simple variant file> <-d rds variant metadata> <-e signature exposures file> <-m metadata file> <-o out path>

Required:
    --gtf       | -g        Path of the gtf file for annotation
    --signature | -s        Name of the signature to test
    --variant   | -v        Location of variants in simple format
    --rds       | -d        Location of variants metadata in RDS format
    --exposures | -e        Location of signature exposures file
                            Column: signature names
                            Rows  : sample ID
    --metadata  | -m        Path to the metadata file
                            <required columns: ID, entity, gender>
    --load      | -l        Path to the results file to load
    --out       | -o        Path for output
    
Options:
    --threads   | -t        Number of threads (default:1)
    --correct   | -C        Correction for sigProfiler output (default:1)
    --entity    | -E        Cut-off for entity filter (default:0.05)
    --help      | -h        Show this help message
"
                            
	option_list <- list(
	  # Base file
	  make_option(c("-l", "--load"), type = "character", dest = "results_file"),
	  make_option(c("-g", "--gtf"), type = "character", dest = "annotation_gtf"),
	  make_option(c("-s", "--signature"), type = "character", dest = "signature_test"),
	  make_option(c("-v", "--variant"), type = "character", dest = "variant_file"),
	  make_option(c("-d", "--rds"), type = "character", dest = "variant_meta"),
	  make_option(c("-e", "--exposures"), type = "character", dest = "signature_file"),
	  make_option(c("--correct-exposures"), type = "character", dest = "correction_signature_file",default=NA),
	  make_option(c("-m", "--metadata"), type = "character", dest = "covar_file"),
	  make_option(c("-o", "--out"), type = "character", dest = "out_path"),
	  make_option(c("-C", "--correct"), type = "character", dest = "correction",default=1),
	  make_option(c("-E", "--entity"), type = "numeric", dest = "minentityposcasespct",default=0.05),
	  make_option(c("-L", "--outliers"), type = "numeric", dest = "outliersThreshold",default=100),
	  make_option(c("-t", "--threads"), type = "numeric", dest = "threads",default=1)
	)
	
	args <- commandArgs(trailingOnly = TRUE)
	help <- (sum(c("--help", "-h") %in% args) >= 1)
	if (help) {
	    cat(help_message)
	    quit()
	}

	argv <- parse_args(OptionParser(option_list = option_list))

	
	showhelp=F
	#print( names(argv))
	if (!("results_file" %in% names(argv))){
		cat("Please provide results from sigDriver\n")
		showhelp=T
	}
	if (!("signature_test" %in% names(argv))){
		cat("Please provide a signature to test\n")
		showhelp=T
	}
	if (!("variant_file" %in% names(argv))){
		cat("Please provide variants metadata file in RDS format\n")
		showhelp=T
	}
	if (!("signature_file" %in% names(argv))){
		cat("Please provide a signature exposures file\n")
		showhelp=T
	}
	if (!("covar_file" %in% names(argv))){
		cat("Please provide a metadata file\n")
		showhelp=T
	}
	if (!("annotation_gtf" %in% names(argv))){
		cat("Please provide a gtf file\n")
		showhelp=T
	}
	if (!("out_path" %in% names(argv))){
		cat("Please provide an output directory\n")
		showhelp=T
	}
	if (showhelp) {
			cat("=======================================\n")
	    cat(help_message)
	    quit()
	}
	results_file = argv$results_file
	signature_test = argv$signature_test
	variant_file = argv$variant_file
	variant_meta = argv$variant_meta
	signature_file = argv$signature_file
	covar_file = argv$covar_file
	annotation_gtf = argv$annotation_gtf
	out_path = argv$out_path
	threads = argv$threads
	sigProfilerInput = as.numeric(argv$correction)
	minentityposcasespct = argv$minentityposcasespct
	outliersThreshold = argv$outliersThreshold

        #check file exists
        if (!file.exists(annotation_gtf)){ stop("GTF file not found") }
        if (!file.exists(results_file)){ stop("Results file not found") }
        if (!file.exists(variant_file)){ stop("Variant file not found") }
        if (!file.exists(covar_file)){ stop("Metadata file not found") }
        if (!file.exists(signature_file)){ stop("Signature file not found") }
	
	#parameters overview
	cat("=======================================\n")
	cat(paste("GTF file         : ",annotation_gtf,"\n",sep=""))
	cat(paste("Results file     : ",results_file,"\n",sep=""))
	cat(paste("Signature to test: ",signature_test,"\n",sep=""))
	cat(paste("Variant file     : ",variant_file,"\n",sep=""))
	cat(paste("Variant metadata : ",variant_meta,"\n",sep=""))
	cat(paste("Signature file   : ",signature_file,"\n",sep=""))
	cat(paste("medata file      : ",covar_file,"\n",sep=""))
	cat(paste("Entity cut-off   : ",minentityposcasespct,"\n",sep=""))
	cat(paste("sigProfiler input: ",sigProfilerInput,"\n",sep=""))
	cat(paste("Outliers cut-off : ",outliersThreshold,"\n",sep=""))
	cat(paste("Output path      : ",out_path,"\n",sep=""))
	cat(paste("Threads          : ",threads,"\n",sep=""))
	cat("=======================================\n")
}else{

	#usage 
	args = commandArgs(trailingOnly=TRUE)

	#called by command
	if (length(commandArgs(trailingOnly=TRUE)) > 0){
		if (length(args) < 5) {
		  stop("Usage: Rscript sigdriver_annotate.R <variant_simple_file> <signature_exposures> <covariates file> <test signature> <output path> <sigdriver results>", call.=FALSE)
		}
		#print variables
		variant_file=args[1]
		signature_file=args[2]
		covar_file=args[3]
		signature_test=args[4]
		out_path=args[5]
		results_file=args[6]
		annotation_gtf=args[7]
	}
}

require(sigDriver)
sigDriver_annotate(signature_test=signature_test,
                   variant_file=variant_file,
                   variant_meta=variant_meta,
                   signature_file=signature_file,
                   covar_file=covar_file,
                   out_path=out_path,
                   results_file=results_file,
                   annotation_gtf=annotation_gtf,
                   threads=threads,
									 entitycuttoff=minentityposcasespct,
									 outliersThreshold=outliersThreshold,
                   sigProfilerInput=sigProfilerInput)
