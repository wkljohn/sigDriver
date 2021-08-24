#!/usr/bin/env Rscript

#parameters parsing library
library(optparse)

#testing
if (F){
	#Rscript sigdriver.R /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv SBS2 /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/
	#APOBECCOMBINE  cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1ORG_combine_NOART_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1ORG_combine_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/ -t 16\"  run_sigdriver.R"}'
	#PCAWG  cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/ -t 16\"  run_sigdriver.R"}'
	#:MELA  cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv "$0" /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/\" sigdriver.R"}' | grep SBS7
	#PCAWGE cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1/ -t 16 \" run_sigdriver.R"}' | grep -v EX7
	#:MELA	cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_E1_raw_exposure.tsv /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv "$0" /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/\" sigdriver.R"}' | grep  EX7
	#PCAWGP cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_PENTA_NOART_raw_exposure.tsv | sed '1d' | awk '{print "bsub -m \"abi-cn24u12 abi-cn24u11 abi-cn24u13 abi-cn24u14 abi-cn24u15 abi-cn24u16 abi-cn24u17 abi-cn24u18 abi-cn24u19 abi-cn24u20 abi-cn24u21 abi-cn24u22 odcf-cn34u03s07 odcf-cn34u03s10 odcf-cn34u09s12 odcf-cn34u18s01 odcf-cn34u18s02 odcf-cn34u18s03 odcf-cn34u18s04 odcf-cn34u18s05 odcf-cn34u18s11 odcf-cn34u18s12 odcf-cn31u17\" -q verylong -n 16 -M 40G  -R \"span[hosts=1] rusage[mem=40GB]\" -o \""$0"\".log -e \""$0".err\"  \"module load r/3.6.0;module load gcc/7.2.0; Rscript run_sigdriver.R -v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_PENTA_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/\""}' | grep -v "SBS7\s"
	#PCAWGPQ cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_PENTA_NOART_raw_exposure.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_PENTA_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/ -t 16\"  run_sigdriver.R"}' | grep -v "SBS7\s"
	#CATCH	Rscript sigdriver.R /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart_consensus.simple.gz /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv SBS13 /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/
	#CATCHP cut -f 1 /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=15g -F \"/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart_consensus.simple.gz /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv "$0" /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/\" sigdriver.R"}'
	#::BSUB cut -f 1 /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv | sed '1d' | awk '{print "bsub  -q verylong -n 16 -m \"abi-cn24u12 abi-cn24u11 abi-cn24u13 abi-cn24u14 abi-cn24u15 abi-cn24u16 abi-cn24u17 abi-cn24u18 abi-cn24u19 abi-cn24u20 abi-cn24u21 abi-cn24u22 odcf-cn34u03s07 odcf-cn34u03s10 odcf-cn34u09s12 odcf-cn34u18s01 odcf-cn34u18s02 odcf-cn34u18s03 odcf-cn34u18s04 odcf-cn34u18s05 odcf-cn34u18s11 odcf-cn34u18s12 odcf-cn31u17\" -M 30G  -R \"span[hosts=1] rusage[mem=30GB]\" -o \"CATCH"$0".log\" -e \"CATCH"$0".err\" \"Rscript sigdriver.R /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart_consensus.simple.gz /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv "$0" /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/\" "}'
	
	#################REVIEW2#####################
	#PCAWG  cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -r /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/SBS1_whitelist_regions -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1/ -t 16\"  run_sigdriver.R"}' | egrep "SBS1 | SBS10| SBS2 |SBS17a|SBS15|SBS16| SBS3 | SBS39 | SBS84 | SBS13"
	#IMMUNEP	cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv | sed '1d' | awk '{print "qsub -l nodes=1:ppn=16 -l mem=40g -F \"-v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.mSBS10.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -r /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/SBS1_whitelist_regions -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1/ -t 16\"  run_sigdriver.R"}' | egrep "SBS84 | SBS9"
	# | egrep " SBS9 | SBS84 | SBS6 | SBS13 | SBS39 | SBS84 | SBS37 "
	# | egrep "SBS1 |SBS10a| SBS2 |SBS17a|SBS15|SBS16| SBS3 | SBS39 | SBS84"
	#bsub_PCAWG: cut -f 1 /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv | sed '1d' | awk '{print "bsub -q verylong -n 16 -M 30G  -R \"span[hosts=1] rusage[mem=30GB]\" -o \""$0"\".log -e \""$0".err\"  \"module load r/3.6.0;module load gcc/7.2.0; Rscript run_sigdriver.R -v /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.noskin.tsv -s "$0" -r /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/SBS1_whitelist_regions -o /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k_cal_binbg1/ -t 16\""}' |  egrep "SBS1 |SBS10a| SBS2 |SBS17a|SBS15|SBS16| SBS3 | SBS39 | SBS84 | SBS9 | SBS13 "
	 
	
	
	variant_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz"
	signature_file="/b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv"
	context_file="/b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv"
	
	covar_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv"
	signature_test="SBS3"
	threads=20
	out_path="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/"
	source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association_release/sigdriver.R")

	variant_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz"
	signature_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv"
	covar_file="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/metadata/CATCH_exposure_SBS_withSPID_forsigdrive.tsv"
	#signature_test="SBS2,SBS13"
	threads=20
	signature_test="SBS13"
	out_path="/b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/test/"
	source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association_release/sigdriver.R")

}

#usage 
args = commandArgs(trailingOnly=TRUE)

#TODO: change script to a parameterized script
#options definitions
if (T){
	#Rscript sigdriver.R -s SBS13 -v /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/somatic_SNVs/simple/CATCH_all_SNV_filtart.simple.gz -e /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/signatures/CATCH/SNVs_consensus/merged_CATCH_exposure_SBS_withSPID_forsigdrive.tsv -m /b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv -o /b06x-isilon/b06x-c/chromothripsis/results/hipo/summary/drivers/CATCH/test/

	help_message <-
"usage: Rscript run_sigdriver.R [options] <-s signature test> <-v simple variant file> <-e signature exposures file> <-m metadata file> <-o out path>

Required:
    --signature | -s        Name of the signature to test
    --variant   | -v        Location of variants file in simple format
    --exposures | -e        Location of signature exposures file
                            Column: signature names
                            Rows  : sample ID
    --metadata  | -m        Path to the metadata file
                            <required columns: ID, entity, gender>
    --out       | -o        Path for output
    
Options:
    --rand      | -a        Random seed (default:1)
    --threads   | -t        Number of threads (default:1)
    --correct   | -C        Correction for sigProfiler output (default:1)
    --regions   | -r        Whitelist file for regions (<chr>:<start>-<end>)
    --help      | -h        Show this help message
"           
	option_list <- list(
	  # Base file
	  make_option(c("-s", "--signature"), type = "character", dest = "signature_test"),
	  make_option(c("-v", "--variant"), type = "character", dest = "variant_file"),
	  make_option(c("-e", "--exposures"), type = "character", dest = "signature_file"),
	  make_option(c("-m", "--metadata"), type = "character", dest = "covar_file"),
	  make_option(c("-o", "--out"), type = "character", dest = "out_path"),
	  make_option(c("-t", "--threads"), type = "numeric", dest = "threads",default=1),
	  make_option(c("-r", "--regions"), type = "character", dest = "testregions",default=NA),
	  make_option(c("-b", "--background"), type = "character", dest = "backgroundsigs",default="SBS1,SBS5,SBS8"),
	  make_option(c("-c", "--context"), type = "character", dest = "context_file",default=NA),
	  make_option(c("-C", "--correct"), type = "character", dest = "correction",default=1),
	  make_option(c("-a", "--rand"), type = "character", dest = "randSeed",default=1),
	  make_option(c("--fp"), type = "character", dest = "fp",default=NA),
	  make_option(c("--fm"), type = "character", dest = "fm",default=NA),
	  make_option(c("--mc"), type = "character", dest = "mc",default=NA)
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
	print(paste(commandArgs()))
	if (!("signature_test" %in% names(argv))){
		cat("Please provide a signature to test\n")
		showhelp=T
	}
	if (!("variant_file" %in% names(argv))){
		cat("Please provide variants file in simple format\n")
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
	if (!("out_path" %in% names(argv))){
		cat("Please provide an output directory\n")
		showhelp=T
	}
	if (showhelp) {
			cat("=======================================\n")
	    cat(help_message)
	    quit()
	}
	
	signature_test = argv$signature_test
	variant_file = argv$variant_file
	signature_file = argv$signature_file
	covar_file = argv$covar_file
	testregions = argv$testregions
	context_file = argv$context_file
	backgroundsigs = argv$backgroundsigs
	out_path = argv$out_path
	threads = argv$threads
	randSeed = argv$randSeed
	sigProfilerInput = as.logical(argv$correction)

  #check file exists
  if (!file.exists(variant_file)){ stop("Variant file not found") }
  if (!file.exists(covar_file)){ stop("Metadata file not found") }
  if (!is.na(testregions) && !file.exists(testregions)){ stop("whitelist regions file not found") }
  if (!is.na(context_file) && !file.exists(context_file)){ stop("context file not found") }
  if (!file.exists(signature_file)){ stop("Signature file not found") }

	#parameters overview
	cat("=======================================\n")
	cat(paste("Signature to test: ",signature_test,"\n",sep=""))
	cat(paste("Variant file     : ",variant_file,"\n",sep=""))
	cat(paste("Signature file   : ",signature_file,"\n",sep=""))
	cat(paste("medata file      : ",covar_file,"\n",sep=""))
	cat(paste("test regions     : ",testregions,"\n",sep=""))
	#cat(paste("Background       : ",backgroundsigs,"\n",sep=""))
	cat(paste("sigProfiler input: ",sigProfilerInput,"\n",sep=""))
	cat(paste("Output path      : ",out_path,"\n",sep=""))
	cat(paste("Threads          : ",threads,"\n",sep=""))
	cat("=======================================\n")
}else{
	#called by command
	if (length(commandArgs(trailingOnly=TRUE)) > 0){
		if (length(args) < 5) {
		  stop("Usage: Rscript sigdriver.R <variant_simple_file> <signature_exposures> <covariates file> <test signature> <output path>", call.=FALSE)
		}
		#print variables
		variant_file=args[1]
		signature_file=args[2]
		covar_file=args[3]
		signature_test=args[4]
		out_path=args[5]
	}
}

if (!file.exists(out_path)){
  print("Creating output directory")
  dir.create(out_path)
}


require(sigDriver)
sigDriver(signature_test=signature_test,
          variant_file,
					signature_file,
					covar_file,
					testregions,
					backgroundsigs,
					out_path,
					threads,
          sigProfilerInput=sigProfilerInput,
					randSeed=randSeed)
