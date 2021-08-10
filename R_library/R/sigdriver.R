#' Driver discovery using mutational signature exposures
#' 
#' @param signature_test Name of the signature to test
#' @param variant_file Path to variant simple file
#' @param signature_file Path to tab separated signatures exposure file
#' @param covar_file Path to the metadata file
#' @param out_path Path for output
#' @param threads Number of threads for running
#' @import data.table
#' @export
sigDriver <- function(signature_test,
											variant_file,
											signature_file,
											covar_file,
											testregions,
											backgroundsigs,
											out_path,
											threads,
											sigProfilerInput=TRUE){
	#Libraries
	require(GenomicRanges)
	require(dplyr)
	require(data.table)
	#require(sigDriver)
	
	#tweak cut-off for signatures with lower median sample mutation load
	if (signature_test == "SBS84" || signature_test == "SBS9" || signature_test == "SBS3"){
		print("lower sample mutation load mode")
		minentityposcasespct = 0.015
		##min_testing_bin_vars = 7
		##frame_pruned_min_nvar = 7
	}
	#minentityposcasespct = 0.01
	#min_testing_bin_vars = 7
	#frame_pruned_min_nvar = 7
	
	#run preparation
	somaticvarranges = read_variants_ranges(variant_file)
	sigexpinfo = read_signature_exposures_matrix(signature_file)
	sampleinfo = read_metadata_matrix(covar_file,covariates)
	sampleinfo = merge_signature_samples(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,signature_test=signature_test,thresholdhypmutation=thresholdhypmutation)
	sampleinfo = merge_allsignature_samples(sampleinfo = sampleinfo,sigexpinfo = sigexpinfo)
	sampleinfo = filter_sample_info_matrix_by_vrange(sampleinfo=sampleinfo,somaticvarranges=somaticvarranges)
	entities_include = get_signature_positive_entities(sampleinfo=sampleinfo,minentityposcasespct=minentityposcasespct,maxentityposcasespct=maxentityposcasespct)
	sampleinfofiltered = filter_sample_info_matrix(sampleinfo=sampleinfo,sigexpinfo=sigexpinfo,entities_include=entities_include)
	
	#test ranking
	sampleinfofiltered$sigRank = rank(sampleinfofiltered$normalized_exposures,ties.method="max")
  sampleinfofiltered$sigRank = sampleinfofiltered$sigRank - min(sampleinfofiltered$sigRank)
  sampleinfofiltered$sigRank = sampleinfofiltered$sigRank / max(sampleinfofiltered$sigRank)
  if (sigProfilerInput){
  	#converge the lower bound
  	#sigMean = mean(sampleinfofiltered$normalized_exposures[sampleinfofiltered$normalized_exposures > 0])
  	#sigSD = sd(sampleinfofiltered$normalized_exposures[sampleinfofiltered$normalized_exposures > 0])
  	#sampleinfofiltered$normalized_exposures[which(sampleinfofiltered$normalized_exposures > sigMean + 4 * sigSD)] = max(sampleinfofiltered$normalized_exposures[which(sampleinfofiltered$normalized_exposures < sigMean + 4 * sigSD)])
  	#sampleinfofiltered$normalized_exposures[sampleinfofiltered$normalized_exposures > 0] = sampleinfofiltered$normalized_exposures[sampleinfofiltered$normalized_exposures > 0] - min(sampleinfofiltered$normalized_exposures[sampleinfofiltered$normalized_exposures > 0]) + 0.001
  	sampleinfofiltered$normalized_exposures = correctExposures(sampleinfofiltered$normalized_exposures)
  }
  
	#termination criterias
	if (dim(sampleinfofiltered)[1] == 0){
		print("no samples left, stop")
		quit()
	}


	#keep only tumors to be tested in variants table
	gns = generate_testing_unit_genomic_bins(somaticvarranges)
	#filter if whitelist exists
  if (!is.na(testregions)) gns=genomic_bins_filter_whitelist(testregions,gns)
	#first filt variant range by sample list
	somaticvarranges = samplefilter_somatic_vranges(somaticvarranges,sampleinfofiltered)
	#reduce the bins to test by variant distance
	gns = prefilter_genomic_bins(gns,somaticvarranges,framesize_pruned,frame_pruned_min_nvar)
	#in turn shrink the variant array by bins to test
	somaticvarranges = prefilter_somatic_vranges(gns,somaticvarranges)
	somaticvarranges = split_variants_GR_by_chr(somaticvarranges) #acceleration by splitting chr, only after whole variant file operations finished
	#sigexpinfo = filter_exposures_matrix(sigexpinfo=sigexpinfo,sampleinfo=sampleinfo)
	
	#generate other signature null distribution in bin
	somaticvarranges = signatureRepresentationAdjustment(gns=gns,signature_test=signature_test,sigexpinfo=sigexpinfo,backgroundsigs=backgroundsigs,somaticvarranges=somaticvarranges,samplemetatablewithentity=sampleinfofiltered,threads=threads)
	somaticvarranges = somaticVariantsProbabalisticSubsampling(somaticvarranges)

	outfile = paste(out_path,"/",signature_test,"_intermediate_results.tsv",sep="")
	write(paste("Gene","regions","markers","marker_snvs","n.samples","Q","p.value","entity","entityvar",sep=" "),file=outfile)	

	#create null model
	require(SKAT)
	print("Generating null model...")
	if (length(unique(sampleinfofiltered$entity)) == 1){
    #nullSKATmodel<-SKAT_Null_Model( sampleinfofiltered$normalized_exposures ~ sampleinfofiltered$total_variants + log2(sampleinfofiltered$total_variants + 1) + sampleinfofiltered$gender, out_type="C",n.Resampling=200,type.Resampling="bootstrap")
    nullSKATmodel<-SKAT_Null_Model( sampleinfofiltered$normalized_exposures ~ log2(sampleinfofiltered$total_variants + 1)  + sampleinfofiltered$gender, out_type="C")#,n.Resampling=2000,type.Resampling="bootstrap")
  }else{
    #nullSKATmodel<-SKAT_Null_Model( sampleinfofiltered$normalized_exposures ~ sampleinfofiltered$entity + sampleinfofiltered$total_variants + log2(sampleinfofiltered$total_variants + 1) + sampleinfofiltered$gender, out_type="C",n.Resampling=200,type.Resampling="bootstrap")
    nullSKATmodel<-SKAT_Null_Model( sampleinfofiltered$normalized_exposures ~ sampleinfofiltered$entity +   log2(sampleinfofiltered$total_variants + 1) +  sampleinfofiltered$gender, out_type="C")#,n.Resampling=2000,type.Resampling="bootstrap")
  }
  
	#init parallelization
	gc()
	if (threads > 1){
		print("Starting association workers...")
		if (grepl("b06",Sys.info()["nodename"])){
		  cl <- parallel::makeCluster(threads,useXDR=TRUE)
		  #cl <- parallel::makeCluster(10,useXDR=FALSE,type="PSOCK")
		}else{
		  cl <- parallel::makeCluster(threads,useXDR=FALSE,type="FORK")
		}
	  
		resultsSKAT=parSapply(cl, 1:length(gns$SYMBOL), doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50,nullSKATmodel=nullSKATmodel)
		#resultsSKAT=parSapply(cl, 1:16, doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50)
		stopCluster(cl)
	}else{
		print("Starting without parallelization...")
		resultsSKAT=lapply(1:length(gns$SYMBOL), doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50,nullSKATmodel=nullSKATmodel)

	}
	
	resultsSKATdf = data.frame(do.call(rbind,resultsSKAT),stringsAsFactors=F)
	colnames(resultsSKATdf) = c("Region","Tested_regions","Sites","Variants","n_cases","Q","p_value","p_value_each","Entities")[1:dim(resultsSKATdf)[2]]
	fulloutpath=paste(out_path,"/",signature_test,"_results_FULL.tsv",sep="")
	Collapsedoutpath=paste(out_path,"/",signature_test,"_results.tsv",sep="")

	if (F){	#test
		igene=which(gns$SYMBOL=="chr1:206858001-206860001")
		#igene=which(gns$SYMBOL=="chr6:18312001-18314001") #1 case, 6 vars problematic frame
		#igene=which(gns$SYMBOL=="chr9:25285001-25287001") #2 case, still failed 6 vars problematic frame
		#igene=1
		samplemetatablewithentity=sampleinfofiltered
		varianttype=50
		pathfile=""
		doassocandwriteSKAThotspot(igene,
			gns=gns,
			somaticvarranges=somaticvarranges,
			outfile=outfile,
			samplemetatablewithentity=sampleinfofiltered,
			sigtest=signature_test,pathfile="",
			varianttype=50,nullSKATmodel=nullSKATmodel) 
			
		
		igene=which(gns$SYMBOL=="chr1:206858001-206860001")
		resultsSKAT=parSapply(cl, (igene-100):(igene+100), doassocandwriteSKAThotspot, gns=gns,somaticvarranges=somaticvarranges,outfile=outfile,samplemetatablewithentity=sampleinfofiltered,sigtest=signature_test,pathfile="",varianttype=50,nullSKATmodel=nullSKATmodel)
		

	}

	#collapse output region
	write.table(resultsSKATdf,fulloutpath,col.names = T,quote = F,row.names=F,sep="\t")
	resultsSKATcollapseddf = collapse_hits_by_window(resultsSKATdf)

	#filter hits
	level1filterVariants = 10
	level1filterRecurrence = 2
	level2filterRecurrence = 4
	resultsSKATcollapseddf = 
		resultsSKATcollapseddf[
			which((resultsSKATcollapseddf$Variants >= level1filterVariants & resultsSKATcollapseddf$Variants - resultsSKATcollapseddf$Sites >= level1filterRecurrence) | 
			      (resultsSKATcollapseddf$Variants - resultsSKATcollapseddf$Sites >= level2filterRecurrence)),]

	#p-adjustment calculation
	resultsSKATcollapseddf$p_adjust_BH = p.adjust(resultsSKATcollapseddf$p_value, method="BH")
	resultsSKATcollapseddf$p_adjust_bonferroni = p.adjust(resultsSKATcollapseddf$p_value, method="bonferroni")
	resultsSKATcollapseddf = resultsSKATcollapseddf[,c("Region",	"Tested_regions",	"Sites",	"Variants",	"n_cases",	"Q",	"p_value",	"p_adjust_BH",	"p_adjust_bonferroni",	"p_value_each",	"Entities")]

	#write hits
	write.table(resultsSKATcollapseddf[,1:11],Collapsedoutpath,col.names = T,quote = F,row.names=F,sep="\t")

	#Summary generator
	plot_qq(resultsSKATcollapseddf,paste(out_path,signature_test,"_",calculate_lambda(resultsSKATcollapseddf),"_qq.png",sep=""))
	print("Done")
}
