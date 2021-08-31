#!/usr/bin/env Rscript


read_metadata_matrix <- function(covar_file,covariates){
	sampleinfo = read.table(covar_file,header = T,stringsAsFactors = F,sep="\t")
	
	#check covariates exists
	if (length(which(covariates %in% colnames(sampleinfo))) != length(covariates)){
		stop(paste("Covariates missing in covar file:",colnames(sampleinfo)))
	}
	
	return(sampleinfo)
}

read_signature_exposures_matrix <- function(signature_file){
	sigexpinfo =  read.table(signature_file,header=T,check.names = F,row.names=1,sep="\t")

	return(sigexpinfo)
}

read_context_frequency_matrix <- function(context_file){
	contextinfo =  read.table(context_file,header=T,check.names = F,row.names=1,sep="\t")

	return(contextinfo)
}

merge_context_samples <- function(sampleinfo,contextinfo){
	sampleinfo = merge(sampleinfo,t(contextinfo),by.x="ID",by.y="row.names")
	colnames(sampleinfo) = gsub("(\\[|\\]|>)","_",colnames(sampleinfo))
	return(sampleinfo)
}

merge_allsignature_samples <- function(sampleinfo,sigexpinfo){
	normsig = (t(sigexpinfo)/colSums(sigexpinfo))
	sampleinfo = merge(sampleinfo,normsig,by.x="ID",by.y="row.names")
	colnames(sigexpinfo) = gsub("(\\[|\\]|>)","_",colnames(sigexpinfo))
	return(sampleinfo)
}
merge_allsignature_byrank_samples <- function(sampleinfo,sigexpinfo){
	normsig = (t(sigexpinfo)/colSums(sigexpinfo))
	normsigrank = normsig
	for (i in 1:dim(normsig)[2]){
		normsigrank[,i] = rank(normsig[,i],ties.method = "min")
	}
	sampleinfo = merge(sampleinfo,normsigrank,by.x="ID",by.y="row.names")
	colnames(sigexpinfo) = gsub("(\\[|\\]|>)","_",colnames(sigexpinfo))
	return(sampleinfo)
}

merge_signature_samples <- function(sampleinfo,sigexpinfo,signature_test,thresholdhypmutation){
	signorminfo = melt(colSums(sigexpinfo))
	signorminfo$ID = rownames(signorminfo)
	colnames(signorminfo)[1] = "total_variants"
	
	if (!grepl(",",signature_test)){
		sigtestexpinfo = melt(sigexpinfo[signature_test,])
	}else{
	  sigsincl = strsplit(signature_test,",")
	  sigsincl = sigsincl[[1]]
	  #print(sigexpinfo[sigsincl,])
	  
	  sigtestexpinfo = data.frame(signature_variants=colSums(sigexpinfo[sigsincl,]),stringsAsFactors=F,check.names=F)
	  sigtestexpinfo = data.frame(ID=rownames(sigtestexpinfo),signature_variants=sigtestexpinfo$signature_variants)
	  #print(sigtestexpinfo)
	}
	colnames(sigtestexpinfo) = c("ID","signature_variants")
	sigtestexpinfo = merge(sigtestexpinfo,signorminfo,by="ID")
	
	sigtestexpinfo$normalized_exposures = sigtestexpinfo$signature_variants / sigtestexpinfo$total_variants
	
	sampleinfomerge = merge(sampleinfo,sigtestexpinfo)
	
	#pre-datafile blacklist can vary results of this part
	#filter hypermutation
	print(paste("Hypermutated samples removed:", length(which(sampleinfomerge$total_variants >= thresholdhypmutation))))
	sampleinfomerge = sampleinfomerge[sampleinfomerge$total_variants < thresholdhypmutation,] #~10mut/Mb of 3,234.83 Mb


	#parameters
	minentityposvalue = 0.10
	
	entitiesin = unique(sampleinfomerge$entity)
	
	signature_variants = sampleinfomerge$signature_variants
	signature_variants = signature_variants[signature_variants > 0]
	
	minsigposAbsExpValueLQ = quantile(signature_variants,0.35,na.rm=T)
	minsigposAbsExpValueMEDIAN = quantile(signature_variants,0.50,na.rm=T)

	sampleinfomerge$ispos = unname(sampleinfomerge$normalized_exposures > minentityposvalue | sampleinfomerge$signature_variants > minsigposAbsExpValueLQ)
	sampleinfomerge$isposMEDIAN = unname(sampleinfomerge$signature_variants > minsigposAbsExpValueMEDIAN)

	
	return(sampleinfomerge)
}

#reader for simple variants file
read_variants_ranges <- function(variant_file){
	require(data.table)
	print("reading variants")
	#Example simple file headers
	#<Entity>	<case_ID> <Cohort_ID> <Genome build> <Variant type> <Chr> <start> <end>	<REF> <ALT>
	#ENTITY  tumor_4177987   PEDPANCAN       GRCh37  SNV     X       102346819       102346819       A  G	1       ICGC
	#reading
  fixvcfinfo = fread(file=variant_file,select=c(1,2,3,4,5,6,7,8 ),header=F,sep="\t")
  
  colnames(fixvcfinfo) = c("Entity","case_ID","Cohort_ID","build","type","chr","start","end")
  
	#make sure chr starts with chr
	if (!grepl("^chr",fixvcfinfo$chr[2])){
	  fixvcfinfo$chr = paste("chr",fixvcfinfo$chr,sep="")
	}
	
	#test chromosomes
	testchr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
	fixvcfinfo = fixvcfinfo[fixvcfinfo$chr %in% testchr,]
	
	#convert to granges
  print("Creating variant Granges")
	vcfcolumnsforranges=c("chr","start","end","case_ID")
	fixvcfinfo$start = as.numeric(fixvcfinfo$start)
	fixvcfinfo$start = as.numeric(fixvcfinfo$end)


  somaticvarranges=makeGRangesFromDataFrame(
  			fixvcfinfo[,..vcfcolumnsforranges],
  			seqnames.field="chr",
  			start.field="start",
  			end.field="end",keep.extra.columns=TRUE)

	#done
	return(somaticvarranges)
}


#reader for simple variants file
read_variants_ranges_withGT <- function(variant_file){
	require(data.table)
	print("reading variants")
	#Example simple file headers
	#<Entity>	<case_ID> <Cohort_ID> <Genome build> <Variant type> <Chr> <start> <end>	<REF> <ALT>
	#ENTITY  tumor_4177987   PEDPANCAN       GRCh37  SNV     X       102346819       102346819       A  G	1       ICGC
	#reading
  fixvcfinfo = fread(file=variant_file,select=c(1,2,3,4,5,6,7,8,9,10 ),header=F,sep="\t")
  
  colnames(fixvcfinfo) = c("Entity","case_ID","Cohort_ID","build","type","chr","start","end","REF","ALT")
  
	#make sure chr starts with chr
	if (!grepl("^chr",fixvcfinfo$chr[2])){
	  fixvcfinfo$chr = paste("chr",fixvcfinfo$chr,sep="")
	}
	
	#test chromosomes
	testchr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
	fixvcfinfo = fixvcfinfo[fixvcfinfo$chr %in% testchr,]
	
	#convert to granges
  print("Creating variant Granges")
	vcfcolumnsforranges=c("chr","start","end","case_ID","REF","ALT")
	fixvcfinfo$start = as.numeric(fixvcfinfo$start)
	fixvcfinfo$start = as.numeric(fixvcfinfo$end)
  somaticvarranges=makeGRangesFromDataFrame(
  			fixvcfinfo[,..vcfcolumnsforranges],
  			seqnames.field="chr",
  			start.field="start",
  			end.field="end",keep.extra.columns=TRUE)

	#done
	return(somaticvarranges)
}


merge_GR <- function(GRList){
	GR=GRanges()
	for (i in 1:length(GRList)){
		GR=c(GR,GRList[[i]])
	}
	return(GR)
}


split_variants_GR_by_chr <- function(somaticvarranges){
	
  print("split somatic vars by chr")
  somaticvarrangesbychr=list()
  lookuparray = data.frame(cbind(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
                                 c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)),stringsAsFactors = FALSE)
  rownames(lookuparray) = lookuparray$X1
  lookuparray$X2 = as.numeric(lookuparray$X2)
  for (i in 1:dim(lookuparray)[1]){
  	inclindexies = which(as.character(seqnames(somaticvarranges)) == lookuparray$X1[i])
  	if (length(inclindexies) > 0){
  		somaticvarrangesbychr[[lookuparray$X2[i]]] = somaticvarranges[which(as.character(seqnames(somaticvarranges)) == lookuparray$X1[i])]
  	}else{
  		somaticvarrangesbychr[[lookuparray$X2[i]]] = c()
  	}
   }
  return(somaticvarrangesbychr)
}

get_signature_positive_entities <- function(sampleinfo,minentityposcasespct,maxentityposcasespct){

	entitysig = table( sampleinfo %>% dplyr::select("entity","ispos")) #samplemetatablewithentity[,c(6,10)])
	entitysigpct = data.frame(entitysig[,2]/(entitysig[,1]+entitysig[,2]),row.names=rownames(entitysig))
	print("Entity signature positivity stats")
	print(entitysigpct)
	entitysigpctnumeric=data.frame(entitysigpct)
	entitysigin = rownames(entitysigpct)[which(entitysigpct[,1]>minentityposcasespct & entitysigpct[,1] < maxentityposcasespct) ]

	return(entitysigin)
}

generate_association_matrix <- function(){


	entitysig = table( samplemetatablewithentity %>% dplyr::select("tumor_type","ispos")) #samplemetatablewithentity[,c(6,10)])
	entitysigpct = data.frame(entitysig[,2]/(entitysig[,1]+entitysig[,2]))
	#decide if suspressinflation should be forced
	entitysigpctnumeric=data.frame(entitysigpct)
	if (varianttype != 31 && varianttype != 30 && varianttype == 10 &&
	    suspressinflation == 0 && length(entitysigpctnumeric[entitysigpctnumeric[,1]<0.02,]) / dim(entitysigpctnumeric)[1] > 0.6){
	  print("Exposure across entities too low, forcing suspress entity bias mode")
	  suspressinflation = 1
	  minentityposcasespct = 0.05
	}
	entities_include = rownames(entitysigpct)[which(entitysigpct[,1]>minentityposcasespct & entitysigpct[,1] < maxentityposcasespct) ]
	
	return(entities_include)
}

filter_sample_info_matrix_by_vrange <- function(sampleinfo,somaticvarranges){
	print(paste("Samples:",dim(sampleinfo)[1]))
	sampleslist=unique(somaticvarranges$case_ID)
	sampleinfo=sampleinfo[which(sampleinfo$ID %in% sampleslist),]
	print(paste("Samples remain:",dim(sampleinfo)[1]))
	return(sampleinfo)
}


filter_sample_info_matrix <- function(sampleinfo,sigexpinfo,entities_include){
	#overlap expinfo with sampleinfo
	sampleinfo = sampleinfo[sampleinfo$ID %in% colnames(sigexpinfo),]
	sampleinfo = sampleinfo[sampleinfo$entity %in% entities_include,]
	
	return(sampleinfo)
}

filter_exposures_matrix <- function(sigexpinfo,sampleinfo){
	
}

#perform data checking
data_validation <- function(){
	validated = T
	return ()
}

read_genomic_bins <- function(sigdriver_results){
	region_split = strsplit(sigdriver_results$Region,":")
	region_splitdf = do.call(rbind,region_split)
	position_split = strsplit(region_splitdf[,2],"-")
	position_splitdf = do.call(rbind,position_split)
	
	gns = makeGRangesFromDataFrame(
					data.frame(seqnames=region_splitdf[,1],start=as.numeric(position_splitdf[,1]),end=as.numeric(position_splitdf[,2]),SYMBOL=sigdriver_results$Region,stringsAsFactors=F),
					seqnames.field="seqnames",start.field="start",end.field="end",keep.extra.columns=TRUE)
  
  return(gns)
}

generate_testing_unit_genomic_bins <- function(somaticvarranges){
  gnsnew = getChrCut(somaticvarranges)
  #ONLY for ACCELERATION OF FIXED WINDOW SIZE=2KB on PCAWG
  if (F){
    blacklistcoordinates = fread("./reference/blacklistcoords_hg19",header=F)
    gnsnew = filter(gnsnew, !SYMBOL %in%  blacklistcoordinates$V1)
  }
  
  gns = makeGRangesFromDataFrame(gnsnew,seqnames.field="seqnames",start.field="start",end.field="end",keep.extra.columns=TRUE)
  
  return(gns)
}


genomic_bins_filter_whitelist <- function(testregions,gns){
	whitelist_regions = read.table(testregions,stringsAsFactors=F)
	whitelist_regions_df = data.frame(do.call(rbind,strsplit(whitelist_regions[,1],"[:|-]")),stringsAsFactors=F)
	colnames(whitelist_regions_df) = c("Chr","start","end")
	whitelist_regions_GR = makeGRangesFromDataFrame(whitelist_regions_df)
	gns = subsetByOverlaps(gns, whitelist_regions_GR, type="equal")
	print(paste("Filtering by whitelist",length(gns)))
	return(gns)
}

samplefilter_somatic_vranges <- function(somaticvarranges,sampleinfo){
	print(paste("Filtering",length(somaticvarranges),"variants by sample"))
	somaticvarranges = somaticvarranges[somaticvarranges$case_ID %in% sampleinfo$ID]
	print(paste("Remaining",length(somaticvarranges),"variants by sample"))
	return(somaticvarranges)
}

prefilter_somatic_vranges <- function(gns,somaticvarranges){
	print(paste("Filtering",length(somaticvarranges),"variants"))
	overlapgns = findOverlaps(gns,somaticvarranges)
	keepsomatic = unique(subjectHits(overlapgns))
	somaticvarranges = somaticvarranges[keepsomatic]
	print(paste("Remaining",length(somaticvarranges),"variants"))
	return(somaticvarranges)
}

prefilter_genomic_bins <- function(gns,somaticvarranges,framesize_pruned,frame_pruned_min_nvar){
	print(paste("Filtering",length(gns),"testing bins"))
	somaticvr_nearest = distanceToNearest(somaticvarranges)
	somaticvr_nearestdf = as.data.table(somaticvr_nearest,stringsAsFactors=F)
	somaticvr_nearestdf = somaticvr_nearestdf[order(somaticvr_nearestdf$queryHits),]
	
	if (length(somaticvarranges) != dim(somaticvr_nearestdf)[1]){
		stop("Distance measurement failed, check variants per chromosome >1")
	}
	
	somaticvarranges$distance = somaticvr_nearestdf$distance
	#min way to put 5 variants inside <framesize> is <framesize> / 4
	somaticvarranges_subset = somaticvarranges[which(somaticvarranges$distance < round(framesize_pruned / frame_pruned_min_nvar))]


	#OLD: filter by recurrence
	if (F){
		recurrentsites = data.frame(table(start(somaticvarranges)),check.names=F,stringsAsFactors=F)
		recurrentsitesstart = recurrentsites[which(recurrentsites$Freq>1),]$Var1
		somaticvarranges_recurrent = somaticvarranges[start(somaticvarranges) %in% recurrentsitesstart]
		
		#gns=gns[ingnstable[which(ingnstable$Freq >= min_testing_bin_vars),]$Var1]
		overlapgns = findOverlaps(gns,somaticvarranges_recurrent)
	}
	
	overlapgns = findOverlaps(gns,somaticvarranges_subset)
	ingnstable = as.data.table(table(queryHits(overlapgns)),stringsAsFactors=F)
	gns=gns[as.numeric(ingnstable[which(ingnstable$N >= (frame_pruned_min_nvar - 1)),]$V1)]
	print(paste("Remaining",length(gns),"testing bins"))
	return(gns)
}

run_sigdriver_association <- function(signature_test,somaticvarranges,sigexpinfo,sampleinfo,gns,out_path,test_mode,write_intermediate){
	
	outfile_intemediate=paste(out_path,"/",pathsuffix,gsub("Signature ","",signature_test),"_intermediate.tsv",sep="")
	outfile_full=paste(out_path,"/",pathsuffix,gsub("Signature ","",signature_test),"_results.tsv",sep="")
	
	#execute association function in parallel
	results=parSapply(cl, 1:length(gns$SYMBOL), 
		doassocandwriteSKAThotspot, 
		gns=gns,
		somaticvarranges=somaticvarranges,
		outfile=outfile_intemediate,
		samplemetatablewithentity=sampleinfo,
		samplevariantscnttable=samplevariantscnttable,
		sigexpinfo=sigexpinfo,
		sigtest=signature_test,
		pathfile="",
		varianttype=test_mode)
	    
	write.table(results,file=outfile_full,col.names = F,quote = F,sep="\t")

}


plot_qq <- function(resultsSKATdf,outpath){
	require(qqman)
	resultsSKATp = as.numeric(resultsSKATdf[resultsSKATdf$p_value != "NA",]$p_value)
	png(outpath)
	qq(resultsSKATp)
	dev.off()
}

calculate_lambda <- function(resultsSKATdf){
	resultsSKATp = as.numeric(resultsSKATdf[resultsSKATdf$p_value != "NA",]$p_value)
	resultsSKATp = resultsSKATp[resultsSKATp < 1]
	z = qnorm(resultsSKATp / 2)
	lambda = round(median(z^2) / 0.454, 3)
	return(lambda)
}




correctExposuresByEntity <- function(sampleinfofiltered,threshold=100,correctLowerBound=FALSE){
	listEntities = unique(sampleinfofiltered$entity)
	for (i in 1:length(listEntities)){
		sampleinfofiltered[which(sampleinfofiltered$entity == listEntities[i]),]$normalized_exposures = correctExposures(sampleinfofiltered[which(sampleinfofiltered$entity == listEntities[i]),]$normalized_exposures,threshold=threshold)
	}
	sampleinfofiltered$normalized_exposures = correctExposures(sampleinfofiltered$normalized_exposures,threshold=threshold*1.4,correctLowerBound=correctLowerBound)
	return (sampleinfofiltered)
}


correctExposures <- function(values,threshold = 100,correctLowerBound=FALSE){
	#min 5 non-zero values to do correction
	nonZeroValues = length(which(values > 0))
	if (nonZeroValues < 6){
		print("too few positives")
		return(values)
	}

	# values=sampleinfofiltered$SBS10a
	require(DDoutlier)
	kMax = 10
	if (nonZeroValues <= 10){
		kMax = nonZeroValues - 1
	}
	#nonneg_exp=data.frame(values,stringAsFactors=F)#[sampleinfofiltered$SBS10a>0],stringAsFactors=F)
	nonneg_exp=data.frame(values,stringAsFactors=F)
	knn_scores <- sort(DDoutlier::KNN_SUM(nonneg_exp))#',k_max=kMax))
	knn_scores = knn_scores[knn_scores>0]
	#threshold = median(knn_scores) + sd(knn_scores[1:round(length(knn_scores)*0.75)]) * threshold
	LBknnscore <- min(knn_scores[knn_scores>0])
	threshold <-  ((median(knn_scores[knn_scores>0])-LBknnscore) * threshold) + LBknnscore
	UQOutliersCnt <- length(which(knn_scores[round(length(knn_scores)/2):length(knn_scores)] > threshold))
	#lof_scores <- DDoutlier::LDF(nonneg_exp)
	#cof_scores <- DDoutlier::COF(nonneg_exp)
	print(paste("N outliers",threshold,UQOutliersCnt))

	valuesRank = rank(values,ties.method="min")
	values = sort(values[values > 0],decreasing=F)
	dist = list()
	for (i in 1:length(values)-1){
		dist[[length(dist)+1]] = c(x=i,y=abs(values[i] - values[i+1]))
	}
	distdf = data.frame(do.call(rbind,dist),stringsAsFactors=F)
	
	
	#########Correction#############
	#quartile definition
	LQloweridx = 1
	LQupperidx = floor(length(values) / 8)
	#UQloweridx = ceiling(length(values)  - length(values) / 8)
	UQloweridx = length(values) - UQOutliersCnt + 1
	UQupperidx = length(values) 
	SDDistMax = 2
	
	#distance based correction
	meanDist = median(abs(distdf$y))
	DistUQ = quantile(abs(distdf$y),0.75)
	SDDist = sd(abs(distdf$y[LQupperidx:UQloweridx]))
	#correct lower 12.5 % of data
#	LQmeanDist = mean(abs(distdf$y[LQloweridx:LQupperidx]))
#	LQCorrFactor = LQmeanDist / meanDist
#	LQCorrector = distdf$y[LQloweridx:LQupperidx]
#	LQMedian = median(LQCorrector)
#	LQCorrector[LQCorrector > meanDist + SDDist * SDDistMax ] = LQMedian
#	distdf$y[LQloweridx:LQupperidx] = LQCorrector
	#correct Upper 12.5 % of data
	#UQmeanDist = mean(abs(distdf$y[UQloweridx:UQupperidx]))
	#UQCorrFactor = UQmeanDist / meanDist
	#UQCorrector = distdf$y[UQloweridx:UQupperidx]
	#UQMedian = median(UQCorrector)
	#UQCorrector[UQCorrector > meanDist + SDDist * SDDistMax ] = UQMedian
	if (UQloweridx <= UQupperidx){
		distdf$y[UQloweridx:UQupperidx] = DistUQ #UQCorrector
	} 
	
	#########rebuild numbers############
	valuescorrlist = list()
	if (correctLowerBound){
		nowvalue = meanDist * 3
	}else{
		nowvalue = min(values)#meanDist
	}
	for (i in 1:dim(distdf)[1]){
		nowvalue = nowvalue + distdf$y[i]
		valuescorrlist[[length(valuescorrlist)+1]] = nowvalue
	}
	valuescorr = do.call(c,valuescorrlist)
	#fill zeros back
	valuescorr = c(valuescorr, rep(0,length(valuesRank) - length(valuescorr)))
	valuescorr = sort(valuescorr)
	valuescorr = valuescorr[valuesRank]
	
	return(valuescorr)
}