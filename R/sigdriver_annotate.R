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
#' @import data.table
#' @export
sigDriver_annotate <- function(signature_test,
											variant_file,
											signature_file,
											covar_file,
											out_path,
											results_file,
											annotation_gtf,
											threads){

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
	print("Starting association workers...")
	if (grepl("b06",Sys.info()["nodename"])){
	  cl <- parallel::makeCluster(6,useXDR=TRUE)
	  #cl <- parallel::makeCluster(10,useXDR=FALSE,type="PSOCK")
	}else{
	  cl <- parallel::makeCluster(6,useXDR=FALSE,type="FORK")
	}
	
	clusterExport(cl, list("getregionTopMutatedRanges", "doassocandwriteSKAThotspotPerm","getWindowNVarWithWeight","getregionTopMutatedRanges"))
	
	resultsSKATanno=parSapply(cl, 1:length(gns$SYMBOL), doassocandwriteSKAThotspotPerm, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)

	#annotate then write table
	resultsimportancedf = importance_output_to_table(resultsSKATanno)
	#annotate table
	gtfref = readAnnotationGTF(annotation_gtf)
	resultsimportancedf = annotate_importance(resultsimportancedf,gtfref)
	writefields = c("Region","test_site","annotation","annotationtype","remaining_sites","marker_snvs","p.value","importance","entityvar")
	write.table(resultsimportancedf[,writefields],fulloutpath,sep="\t",quote=F,row.names=F)

	#plot importance
	plot_lolli(resultsSKATanno,out_path)

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