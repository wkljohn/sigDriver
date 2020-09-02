library(sigDriver)
source("sigdriver_annotate.R")
source("./lib/annotate_sites.R")
source("./lib/association_SKAT_hotspot_V4.R")
source("./lib/association_SKAT_hotspot_V4Perm.R")
source("./lib/collapse_hits.R")
source("./lib/collapse_hits.R")
source("sigdriver_functions.R")
source("sigdriver_variables.R")

variant_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Mutation_catalogue/annovar/somatic_annotated_V2/merged.avout.hg19_multianno.nofilt.noartsnv.simple.gz"
signature_file="/b06x-isilon/b06x-c/chromothripsis/software/sigProfiler/python_binary/for_testing/output/ICGC_VCF_ORG_NOART_raw_exposure.tsv"
covar_file="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/datafile/merged_pheno_tab.V3.tsv"
signature_test="SBS13"
threads=1
out_path="/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/testing/release/output/"
results_file=paste("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/bin1k/",signature_test,"_RAD51B_results.tsv",sep="")
annotation_gtf="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation.gtf"


sigDriver_annotate(signature_test=signature_test,
                   variant_file=variant_file,
                   signature_file=signature_file,
                   covar_file=covar_file,
                   out_path=out_path,
                   results_file=results_file,
                   annotation_gtf=annotation_gtf,
                   threads=threads)


	#libraries
	library(GenomicRanges,quietly=T)
	library("VariantAnnotation")
	library(data.table,quietly=T)
	library(dplyr,quietly=T)
	library(rtracklayer,quietly=T)


	#run preparation, small differences with sigdriver main routine
	sigdriver_results = read.table(results_file,sep="\t",header=T,stringsAsFactors=F)
	sigdriver_results = sigdriver_results[which(sigdriver_results$p_adjust_BH < 0.05),]
	somaticvarranges = read_variants_ranges_withGT(variant_file)
	sigexpinfo = read_signature_exposures_matrix(signature_file)
	sampleinfo = read_metadata_matrix(covar_file,covariates)
	sampleinfo = merge_signature_samples(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,signature_test=signature_test,thresholdhypmutation=thresholdhypmutation)
	sampleinfo = filter_sample_info_matrix_by_vrange(sampleinfo=sampleinfo,somaticvarranges=somaticvarranges)
	entities_include = get_signature_positive_entities(sampleinfo=sampleinfo,minentityposcasespct=minentityposcasespct,maxentityposcasespct=maxentityposcasespct)
	sampleinfofiltered = filter_sample_info_matrix(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,entities_include=entities_include)

	#keep only tumors to be tested in variants table
	gns = read_genomic_bins(sigdriver_results)
	#first filt variant range by sample list
	somaticvarranges = samplefilter_somatic_vranges(somaticvarranges,sampleinfofiltered)
	#reduce the bins to test by variant distance
	gns = prefilter_genomic_bins(gns,somaticvarranges,framesize_pruned,frame_pruned_min_nvar)
	#in turn shrink the variant array by bins to test
	somaticvarranges = prefilter_somatic_vranges(gns,somaticvarranges)
	somaticvarranges = split_variants_GR_by_chr(somaticvarranges) #acceleration by splitting chr, only after whole variant file operations finished

	#run main
	outfile = paste(out_path,"/",signature_test,"_annotation_intermediate_results.tsv",sep="")
	write(paste("Region","test_site","subregions","remaining_sites","marker_snvs","n.samples","Q","p.value","importance","entityvar",sep=" "),file=outfile)	 
	fulloutpath=paste(out_path,"/",signature_test,"_perturb_sites_annotated_results.tsv",sep="")


	#init parallelization
	gc()
	if (threads > 1){
		print("Starting association workers...")
		if (grepl("b06",Sys.info()["nodename"])){
		  cl <- parallel::makeCluster(6,useXDR=TRUE)
		  #cl <- parallel::makeCluster(10,useXDR=FALSE,type="PSOCK")
		}else{
		  cl <- parallel::makeCluster(6,useXDR=FALSE,type="FORK")
		}
		resultsSKATanno=parSapply(cl, 1:length(gns$SYMBOL), doassocandwriteSKAThotspotPerm, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)

	}else{
		print("Starting without parallelization...")
		resultsSKATanno=lapply(1:length(gns$SYMBOL), doassocandwriteSKAThotspotPerm, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)

	}
	
	#annotate then write table
	resultsimportancedf = importance_output_to_table(resultsSKATanno)
	#annotate table
	gtfref = readAnnotationGTF(annotation_gtf)
	resultsimportancedf = annotate_importance(resultsimportancedf,gtfref)
	writefields = c("Region","test_site","annotation","annotationtype","remaining_sites","marker_snvs","p.value","importance","entityvar")
	write.table(resultsimportancedf[,writefields],fulloutpath,sep="\t",quote=F,row.names=F)

	#plot importance
	plot_lolli(resultsSKATanno,
		out_path=out_path,
		somaticvarranges=somaticvarranges,
		resultsimportancedf=resultsimportancedf,
		gtfref=gtfref)
	print("Done")