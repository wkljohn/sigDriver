#' sigDriver results refinement and annotation
#' 
#' @param signature_test Name of the signature to test
#' @param variant_file Path to variant simple file
#' @param signature_file Path to tab separated signatures exposure file
#' @param covar_file Path to the metadata file
#' @param out_path Path for output
#' @param results_file Path to the results file to load
#' @param annotation_gtf Path to the gtf file for annotation
#' @param threads Number of threads for running
#' @param entitycuttoff cut-off for entity positivity
#' @param sigProfilerInput(default=TRUE) Apply adjustments for sigProfiler signature exposure outputs
#' @param randSeed(default=1) Random seed
#' @import data.table
#' @export
sigDriver_annotate <- function(signature_test,
											variant_file,
											variant_meta,
											signature_file,
											covar_file,
											out_path,
											results_file,
											annotation_gtf,
											threads,
											entitycuttoff,
											outliersThreshold,
											sigProfilerInput=TRUE,
											verbose=TRUE){

	#libraries
	library(GenomicRanges,quietly=T)
	library("VariantAnnotation")
	library(data.table,quietly=T)
	library(dplyr,quietly=T)
	library(rtracklayer,quietly=T)

	minentityposcasespct = entitycuttoff 
	corrOutliersThres = outliersThreshold
	corrLB = FALSE
	
	#tweak cut-off for signatures with lower median sample mutation load
	if (signature_test == "SBS84" || signature_test == "SBS9"|| signature_test == "Signature_EX11"|| signature_test == "SBS-E9"){
		print("lower sample mutation load mode")
		minentityposcasespct = 0.02
		min_testing_bin_vars = 7
		frame_pruned_min_nvar = 7
	}

	#run preparation, small differences with sigdriver main routine
	sigdriver_results = read.table(results_file,sep="\t",header=T,stringsAsFactors=F)
	sigdriver_results = sigdriver_results[which(sigdriver_results$p_adjust_BH < 0.05),]
	somaticvarranges = read_variants_ranges_withGT(variant_file)
	sigexpinfo = read_signature_exposures_matrix(signature_file)
	sampleinfo = read_metadata_matrix(covar_file,covariates)
	sampleinfo = merge_signature_samples(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,signature_test=signature_test,thresholdhypmutation=thresholdhypmutation)
	#merge with signatures for correction
	sampleinfo = merge_allsignature_samples(sampleinfo = sampleinfo,sigexpinfo = sigexpinfo)
	sampleinfo = filter_sample_info_matrix_by_vrange(sampleinfo=sampleinfo,somaticvarranges=somaticvarranges)
	entities_include = get_signature_positive_entities(sampleinfo=sampleinfo,minentityposcasespct=minentityposcasespct,maxentityposcasespct=maxentityposcasespct)
	sampleinfofiltered = filter_sample_info_matrix(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,entities_include=entities_include)

	if (sigProfilerInput){
	 	sampleinfofiltered = correctExposuresByEntity(sampleinfofiltered,threshold=corrOutliersThres,correctLowerBound=corrLB)
  }
	#keep only tumors to be tested in variants table
	gns = read_genomic_bins(sigdriver_results)
	#reduce the bins to test by variant distance
	somaticvarranges_lastrun = merge_GR(readRDS(variant_file))
	somaticvarranges = intersetGRs(somaticvarranges_lastrun,somaticvarranges)
	gns = prefilter_genomic_bins(gns,somaticvarranges,framesize_pruned,frame_pruned_min_nvar)
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
		  cl <- parallel::makeCluster(6,useXDR=TRUE, outfile='/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/logs/sigDriverAnno_info_parallel.log')
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

	if (F){	#test
		write(paste("Region","test_site","subregions","test_sites","marker_snvs","n.samples","Q","p.value","importance","entityvar",sep=" "),file=outfile)	 
		#igene=which(gns$SYMBOL=="chr1:206858001-206860001").
		igene=which(gns$SYMBOL=="chr20:57617001-57619001")	#for SBS2 single variant bug
		#igene=which(gns$SYMBOL=="chr3:183272001-183274001")	#for SBS9 importance calculation
		#igene=which(gns$SYMBOL=="chr2:89184001-89186001")	#for SBS9 importance calculation
		#igene=which(gns$SYMBOL=="chr6:18312001-18314001") #1 case, 6 vars problematic frame
		#igene=which(gns$SYMBOL=="chr9:25285001-25287001") #2 case, still failed 6 vars problematic frame
		#igene=1
		samplemetatablewithentity=sampleinfofiltered
		varianttype=50
		pathfile=""
		a=doassocandwriteSKAThotspot(igene,
			gns=gns,
			somaticvarranges=somaticvarranges,
			outfile=outfile,
			samplemetatablewithentity=sampleinfofiltered,
			sigtest=signature_test,pathfile="",
			varianttype=50) 
			
		
		resultsSKATanno=parSapply(cl, 1:3, doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)
		resultsSKATanno=sapply( 1:length(gns$SYMBOL), doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)
		#failed at SBS2:31
		#do.call(function,resultsSKATanno)
		#do.call(rbind,sapply(resultsSKATanno, function(my.anno)
	  #  do.call(rbind, my.anno)))
	}
}