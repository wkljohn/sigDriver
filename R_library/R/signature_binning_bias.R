somaticVariantsProbabalisticSubsampling <- function(somaticvarranges,rand_seed=1){
	set.seed(rand_seed)
	#probablistic subsampling of variants associated to signatures
	for (i in 1:length(somaticvarranges)){
		print(paste(i,length(somaticvarranges[[i]])))
		somaticvarranges[[i]]$rand=runif(length(somaticvarranges[[i]]))
		somaticvarranges[[i]] = somaticvarranges[[i]][which(somaticvarranges[[i]]$rand < somaticvarranges[[i]]$WEIGHT)]
		if (length(which(somaticvarranges[[i]]$WEIGHT < 1)) > 0){
			somaticvarranges[[i]][which(somaticvarranges[[i]]$WEIGHT < 1)]$WEIGHT = 1
		}
		print(paste(i,length(somaticvarranges[[i]])))
	}
	return(somaticvarranges)
}


signatureRepresentationAdjustment <- function(gns,
                             signature_test,
                             sigexpinfo,
                             backgroundsigs,
                             excludeSigs,
                             somaticvarranges,
                             samplemetatablewithentity,
                             threads,
                             variantFactor=1.5,
                             entityFactor=1,
                             binStatsList=NULL,
                             verbose=TRUE){
   print("Calculating signature-bin distribution")
   print(paste(variantFactor,entityFactor))
		if (is.null(binStatsList)){
			cl <- parallel::makeCluster(threads,useXDR=TRUE)#, outfile='/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/logs/bin_bias_info_parallel.log')
			binStatsList=parLapply(cl, 1:length(gns$SYMBOL), getBinSignatureAverage, gns=gns,somaticvarranges=somaticvarranges,samplemetatablewithentity=samplemetatablewithentity,sigweight=FALSE,sigNames=rownames(sigexpinfo))
			#getBinSignatureAverage(1,gns=gns,somaticvarranges=somaticvarranges,samplemetatablewithentity=samplemetatablewithentity,sigweight=FALSE,sigNames=rownames(sigexpinfo))
			stopCluster(cl)
		}
		#binStatsDF = do.call(rbind,binStatsList)
		binStatsDF = do.call(rbind,lapply(binStatsList,function(x){if (!is.null(x)){return(x[[1]])}}))
		binStatsDFAVG = do.call(rbind,lapply(binStatsList,function(x){if (!is.null(x)){return(x[[2]])}}))
		binStatsENTITY = do.call(rbind,lapply(binStatsList,function(x){if (!is.null(x) ){return(x[[3]])}}))
		binEntityOccr = table(binStatsENTITY) / length(binStatsENTITY)
		
		#Compute non-zero averages
		#sampleSignatureSigExpAvg = samplemetatablewithentity[,rownames(sigexpinfo)]
		#sampleSignatureSigExpAvg[sampleSignatureSigExpAvg == 0] = NA
		#sigExpPositivityInSamples = colMeans(sampleSignatureSigExpAvg,na.rm=T)
		
#		#Compute signature exposures in samples
		sigExpPositivityInSamples = colSums(samplemetatablewithentity[,rownames(sigexpinfo)] > 0.05) / dim(samplemetatablewithentity)[1]

		#rank samples by exposures
		samplemetatablewithentityRank = samplemetatablewithentity[,!colnames(samplemetatablewithentity) %in% rownames(sigexpinfo)]
		samplemetatablewithentityRank = merge_allsignature_byrank_samples(sampleinfo = samplemetatablewithentityRank, sigexpinfo = sigexpinfo)
		
		#entity correction
#		if (length(unique(samplemetatablewithentity$entity)) > 0){
#			sampleEntityOccurence = table(samplemetatablewithentity$entity) / length(samplemetatablewithentity$entity)
#			sampleEntityOccurenceDiff = cbind(sampleEntityOccurence,round(binEntityOccr[names(sampleEntityOccurence)],10))
#			sampleEntityOccurenceCorr = data.frame(((sampleEntityOccurenceDiff[,1]+0.001)/(sampleEntityOccurenceDiff[,2]+0.001)) ^ entityFactor)	
#			sampleEntityOccurenceCorr[is.na(sampleEntityOccurenceCorr)] = 1
#			print("Entity corr")
#			print(sampleEntityOccurenceCorr)
#			somaticvarranges = variantsAddEntity(somaticvarranges,samplemetatablewithentity,sampleEntityOccurenceCorr)
#		}else{
#			sampleEntityOccurenceCorr = 1
#		}
		
		#Computer exposures in bins and create weight
		binStatsMeans = data.frame(colMeans(binStatsDF,na.rm=T),stringsAsFactors=F)
		#method 1, weight only by prevalence
		#binSignatureWeights = data.frame((1-binStatsMeans)^2)
		#method 2, expectation diff
		#binSignatureWeights =  data.frame(1/(100^(binStatsMeans-sigExpPositivityInSamples)))
		#print("sig pos")
		#print(sigExpPositivityInSamples)
		#print("bin pos")
		#print(binStatsMeans)
		#binSignatureWeights =  data.frame(((sigExpPositivityInSamples+0.001)/(binStatsMeans+0.001)) ^1.5)	#0.001 is the error
		binSignatureWeights =  data.frame(((sigExpPositivityInSamples+0.0001)/(binStatsMeans+0.0001)) ^ variantFactor)	
		#center weights, converge mode
#		for (i in 1:10){
#			binNormFactor = mean(log2(binSignatureWeights[,1])*sigExpPositivityInSamples)
#			binSignatureWeights = 2^(log2(binSignatureWeights) - binNormFactor)
#			print(mean(log2(binSignatureWeights[,1])*sigExpPositivityInSamples))
#		}
		if (verbose){
			print("sigExpPositivityInSamples:")
			print(sigExpPositivityInSamples)
			print("binStatMeans:")
			print(binStatsMeans)
		}
		#for average
		#binSignatureWeights =  data.frame(((sigExpPositivityInSamples+0.001)/(binStatsMeans+0.001)) ^ 3)	#0.001 
		#binSignatureWeights[binSignatureWeights > 3] = 3
		#get list of backgrounds, do not weight on background
		backgroundsigslist=c()
		if (nchar(backgroundsigs) > 0){
			backgroundsigslist = strsplit(gsub("\\s","",backgroundsigs),",")[[1]]
			#weighting of background
			binSignatureWeights[which(rownames(binSignatureWeights) %in% backgroundsigslist),] = 1
		}


		#overweight disabled
		if (verbose){
			print("remove overweight:")
			print(binSignatureWeights[which(binSignatureWeights > 1),])
		}
		binSignatureWeights[which(binSignatureWeights > 1),]   = 1
		#weighting of signature
		#upweight underrepresented signatures
		ThresUnderPositivity = 0.5 ^ variantFactor	#disabled when FALSE
		rankscaler=20	#disabled when FALSE
		underRepresentationThreshold = 0.5
		if (FALSE){ #sigExpPositivityInSamples[signature_test] < ThresUnderPositivity){
			ThBasline = ThresUnderPositivity + 0.001
			#adjustment for under-and-over represented signatures
			baselineUnderRepWeight = 10^(log10(1.5) + 
				log10(binSignatureWeights[which(rownames(binSignatureWeights) %in% signature_test),]) + 
				log10(( ThBasline - sigExpPositivityInSamples[signature_test]) / ThresUnderPositivity * 1.7))
			binSignatureWeights[which(rownames(binSignatureWeights) %in% signature_test),] = baselineUnderRepWeight #+  ( ThresUnderPositivity - sigExpPositivityInSamples[signature_test]) / ThresUnderPositivity * 1.7
			#uprank
			samplemetatablewithentityRank[,signature_test] = round(samplemetatablewithentityRank[,signature_test] *  ( 1+10^(log10(1.5) + log10(ThBasline - sigExpPositivityInSamples[signature_test]))*rankscaler))
			#* (1.5 +  ( ThresUnderPositivity - sigExpPositivityInSamples[signature_test]) / ThresUnderPositivity * 1.2 ))
			samplemetatablewithentityRank[,signature_test] = samplemetatablewithentityRank[,signature_test] - min(samplemetatablewithentityRank[,signature_test]) + 1
		}else if (signature_test %in% rownames(binSignatureWeights)  &&
						 (binSignatureWeights[which(rownames(binSignatureWeights) %in% signature_test),] > underRepresentationThreshold || binStatsMeans[signature_test,] < 0.1)){
			#weighting of tested-signature, by conditioned on weighting 
			binSignatureWeights[which(rownames(binSignatureWeights) %in% signature_test),] = -1
		}
		
		#exclude signatures
		if (nchar(excludeSigs) > 0){
			excludeSigslist = strsplit(gsub("\\s","",excludeSigs),",")[[1]]
			#weighting of background
			print("Exclude signatures from correction:")
			print(rownames(binSignatureWeights)[which(rownames(binSignatureWeights) %in% excludeSigslist)])
			binSignatureWeights[which(rownames(binSignatureWeights) %in% excludeSigslist),] = -1
		}
		print("bin weight")
		print(binSignatureWeights)
		
		
		
		#test binning weighting results
		#compute weight on each variant
		sigweight=binSignatureWeights
		
		print("Computing variant-wise weight")
		sigweight=binSignatureWeights
		
		backgroundsigsidx=c()
		if (length(backgroundsigslist) > 0){
			backgroundsigsidx = which(rownames(sigexpinfo) %in%  c(signature_test,backgroundsigslist))
		}
		#preprocessing on matricies
		#sigSampleInfoMatrix = samplemetatablewithentity[,rownames(sigexpinfo)]
		#signature positivity matrix
		#sigSampleInfoMatrix[sigSampleInfoMatrix>0.05] = 1
		#sigSampleInfoMatrix[sigSampleInfoMatrix<=0.05] = 0
		#######
		#Rank
		sigSampleInfoMatrix = samplemetatablewithentityRank[,rownames(sigexpinfo)]
		rownames(sigSampleInfoMatrix) =  samplemetatablewithentity$ID
		
		cl <- parallel::makeCluster(threads,useXDR=TRUE)#, outfile='/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/logs/bin_bias_info_parallel.log')
		#somaticvarranges=parSapply(cl, 1:length(somaticvarranges), signatureRepresentationWeight, somaticvarranges=somaticvarranges,sigweight=sigweight,backgroundsigsidx=backgroundsigsidx,sigSampleInfoMatrix=sigSampleInfoMatrix,sigexpinfo=sigexpinfo)
		somaticvarranges=lapply(1:length(somaticvarranges), signatureRepresentationWeight, somaticvarranges=somaticvarranges,sigweight=sigweight,backgroundsigsidx=backgroundsigsidx,sigSampleInfoMatrix=sigSampleInfoMatrix,sigexpinfo=sigexpinfo)
		stopCluster(cl)
		
		return (somaticvarranges)
}

variantsAddEntity <- function(somaticvarranges,samplemetatablewithentity,EntityWeight){
	for (i in 1:length(somaticvarranges)){
		rownames(samplemetatablewithentity) = samplemetatablewithentity$ID
		somaticvarranges[[i]]$entityWeight = EntityWeight[samplemetatablewithentity[somaticvarranges[[i]]$case_ID,]$entity,1]
	}
	return(somaticvarranges)
}

signatureRepresentationWeight <- function(i,somaticvarranges,sigweight,
																			backgroundsigsidx,
																			sigSampleInfoMatrix,
																			sigexpinfo){
		require(dplyr)
		print(paste("weighting chr",i))
		somaticvarrangesDF = as.data.frame(somaticvarranges[[i]])
		if (dim(somaticvarrangesDF) == 0){
			return(somaticvarranges[[i]])
		}
		somaticvarrangesDF$WEIGHT = 1
		#somaticvarranges[[i]]$WEIGHT = 1
		pos_unique_sites = unique(start(somaticvarranges[[i]]))
		
		#add exposures
		somaticvarrangesDF = merge(somaticvarrangesDF, sigSampleInfoMatrix, by.x="case_ID", by.y="row.names")
		vrange_summary_site = somaticvarrangesDF %>%
	  group_by(start) %>%
	  summarise_at(vars(rownames(sigexpinfo)), mean)
	  vrange_summary_site_sig_pos = vrange_summary_site
	  #vrange_summary_site_sig_pos[,backgroundsigsidx+1] = 0	#not in rank
	  vrange_summary_site$max_exp = do.call(pmax,  vrange_summary_site[,rownames(sigexpinfo)])
	  
	  #for rank, max to binary
	  vrange_summary_site[,rownames(sigexpinfo)] = vrange_summary_site[,rownames(sigexpinfo)] - vrange_summary_site$max_exp + 1
	  vrange_summary_site[vrange_summary_site < 0] = NA
	  vrange_summary_site[,rownames(sigexpinfo)] = (1-vrange_summary_site[,rownames(sigexpinfo)]) + (t(t(vrange_summary_site[,rownames(sigexpinfo)]) * sigweight[,1]))
	  # (1-vrange_summary_site[756,rownames(sigexpinfo)]) + (vrange_summary_site[756,rownames(sigexpinfo)] * sigweight[,1])
	  # vrange_summary_site[vrange_summary_site$start == 191784318,]
	  #OLD
	  #vrange_summary_site[vrange_summary_site<0.0001] = 1
	  #vrange_summary_site[,backgroundsigsidx+1] = 1
	  #NEW
	  #vrange_summary_site[0.9999 < vrange_summary_site  & vrange_summary_site < 1.0001] = NA
	  #vrange_summary_site[,backgroundsigsidx+1] = NA

	  #weighting method
	  #method 1: minimum
	  #vrange_summary_site$useweight =   do.call(pmin, c( vrange_summary_site[,rownames(sigexpinfo)], list( na.rm=T)))
	  #method 2: mean
	  #vrange_summary_site$useweight =   rowMeans(vrange_summary_site[,rownames(sigexpinfo)],na.rm=T)
	  #vrange_summary_site$useweight[is.nan(vrange_summary_site$useweight)] = 1
	  #method 3: mean of exposure max
#	  vrange_summary_site_sig_pos$max_exp = do.call(pmax,  vrange_summary_site_sig_pos[,rownames(sigexpinfo)])
#	  vrange_summary_site_sig_pos = vrange_summary_site_sig_pos - vrange_summary_site_sig_pos$max_exp
#	  vrange_summary_site_sig_pos[vrange_summary_site_sig_pos>=0] = 1
#	  vrange_summary_site_sig_pos[vrange_summary_site_sig_pos<0] = NA
#	  vrange_summary_site = vrange_summary_site * vrange_summary_site_sig_pos	#keep only max
#	  vrange_summary_site$maxWeight = do.call(pmax,c(vrange_summary_site[,rownames(sigexpinfo)],list(na.rm=T)))
#	  vrange_summary_site$minWeight = do.call(pmin,c(vrange_summary_site[,rownames(sigexpinfo)],list(na.rm=T)))
#	  vrange_summary_site$maxWeight[vrange_summary_site$maxWeight < 1] = 1
#	  vrange_summary_site$minWeight[vrange_summary_site$minWeight > 1] = 1
#	  vrange_summary_site$useweight = vrange_summary_site$maxWeight * vrange_summary_site$minWeight
#	  #vrange_summary_site$useweight =   rowMeans(vrange_summary_site[,rownames(sigexpinfo)],na.rm=T)
#	  vrange_summary_site$useweight[is.nan(vrange_summary_site$useweight)] = 1
#	  #vrange_summary_site[vrange_summary_site$isTestsigMax,]$useweight = 1
	  #method 4: mean of weighting max
	  #vrange_summary_site$maxWeight = do.call(pmax,c(vrange_summary_site[,rownames(sigexpinfo)],list(na.rm=T)))
	  #vrange_summary_site$minWeight = do.call(pmin,c(vrange_summary_site[,rownames(sigexpinfo)],list(na.rm=T)))
	  #vrange_summary_site$maxWeight[vrange_summary_site$maxWeight < 1] = 1
	  #vrange_summary_site$minWeight[vrange_summary_site$minWeight > 1] = 1
	  #vrange_summary_site$useweight = vrange_summary_site$maxWeight * vrange_summary_site$minWeight
	  #vrange_summary_site$useweight =   rowMeans(vrange_summary_site[,rownames(sigexpinfo)],na.rm=T)
	  #vrange_summary_site$useweight[is.nan(vrange_summary_site$useweight)] = 1
	  #method 5: mean rank of exposures * weight
	  vrange_summary_site$useweight = do.call(pmin,c(vrange_summary_site[,rownames(sigexpinfo)],list(na.rm=T)))
	  
	  
	 	
	  #entity weight
	  if ("entityWeight" %in% colnames(somaticvarrangesDF)){
	  	vrange_summary_site_entity = somaticvarrangesDF %>%
				  group_by(start) %>%
				  summarise_at("entityWeight", mean)
			vrange_summary_site$useweight = vrange_summary_site$useweight * vrange_summary_site_entity$entityWeight
	  }
	  
	  vrange_summary_site$useweight[vrange_summary_site$useweight < 0] = 1
	  vrange_summary_site = data.frame(vrange_summary_site,stringsAsFactors=F,row.names=vrange_summary_site$start)
	 	
	  somaticvarranges[[i]]$WEIGHT = vrange_summary_site[as.character(start(somaticvarranges[[i]])),]$useweight
	  return(somaticvarranges[[i]])
}

getBinSignatureRepresentation <- function (igene,
                             gns,
                             somaticvarranges,
                             samplemetatablewithentity,
                             sigweight=FALSE,sigNames){
                             
	  require(GenomicRanges)
	  require(dplyr)
	  require(reshape2)
	  require(matrixStats)
	  require(sigDriver)
	  out <- tryCatch({
		  splitByOverlaptolist <-
		    function(query, subject, column="ENTREZID", ...)
		    {
		      olaps <- findOverlaps(query, subject, ...)
		      return(subjectHits(olaps))
		    }
	    
			#annotation variables
		  restricted = TRUE
			lookuparray = data.frame(cbind(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
		                                 c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)),stringsAsFactors = FALSE)
		  rownames(lookuparray) = lookuparray$X1
		  lookuparray$X2 = as.numeric(lookuparray$X2)
		  genetest = gns$SYMBOL[igene]
	    testgns1gene = gns[which(gns$SYMBOL == genetest)]
	    
	    idxchr = lookuparray[as.character(seqnames(testgns1gene)[1]),]$X2
	    variantsassociated = somaticvarranges[[idxchr]][somaticvarranges[[idxchr]]$case_ID %in% samplemetatablewithentity$ID, ]
	    listvarinregion = unique(splitByOverlaptolist(testgns1gene, variantsassociated, "SYMBOL"))
	    topregions <- getregionTopMutatedRanges(testgns1gene,
	                                            variantsassociated[listvarinregion,], #somaticvarranges[[idxchr]],
	                                            samplemetatablewithentity$ID,
	                                            samplemetatablewithentity,
	                                            pctin,restricted,useSigWeight=sigweight)
			if (dim(topregions)[1] != 0){
	      print(paste("extract variants for top:", dim(topregions)[1]))
		    #get what samples are in the bin       
		    topsomaticvrangesdf = cbind(as.character(seqnames(testgns1gene)),topregions$START,topregions$END)
		    topsomaticvrangesdf = data.frame(topsomaticvrangesdf,stringsAsFactors=FALSE)
		    colnames(topsomaticvrangesdf)  = c("CHROM","POS","END")                           
				topsomaticvranges = makeGRangesFromDataFrame(topsomaticvrangesdf,seqnames.field="CHROM",start.field="POS",end.field="END",keep.extra.columns=TRUE)
		    listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, variantsassociated, "SYMBOL"))
		    varframe = data.frame( variantsassociated[listvarinregion,])
		    #print(varframe$case_ID)
		    #print(head(samplemetatablewithentity[varframe$case_ID,sigNames]))
		    #calculate prevalence of signatures in a bin
		   	#print(head(samplemetatablewithentity))
		    #print(head(varframe))
		    sample_sig_table_pos5 = colSums(samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),sigNames]>0.05) / length(varframe$case_ID)
		    #print(sample_sig_table_pos5)#samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),sigNames]>0.05)
		    return(sample_sig_table_pos5)
	    }else{
	    	return(NULL)
	    }
	  },
  error=function(cond) {
    print(paste("Error:",cond))
    return(NULL)
    #return(cond)
  })
}


getBinSignatureAverage <- function (igene,
                             gns,
                             somaticvarranges,
                             samplemetatablewithentity,
                             sigweight=FALSE,sigNames){
                             
	  require(GenomicRanges)
	  require(dplyr)
	  require(reshape2)
	  require(matrixStats)
	  require(sigDriver)
	  out <- tryCatch({
		  splitByOverlaptolist <-
		    function(query, subject, column="ENTREZID", ...)
		    {
		      olaps <- findOverlaps(query, subject, ...)
		      return(subjectHits(olaps))
		    }
	    
			#annotation variables
		  restricted = TRUE
			lookuparray = data.frame(cbind(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
		                                 c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)),stringsAsFactors = FALSE)
		  rownames(lookuparray) = lookuparray$X1
		  lookuparray$X2 = as.numeric(lookuparray$X2)
		  genetest = gns$SYMBOL[igene]
	    testgns1gene = gns[which(gns$SYMBOL == genetest)]
	    
	    idxchr = lookuparray[as.character(seqnames(testgns1gene)[1]),]$X2
	    variantsassociated = somaticvarranges[[idxchr]][somaticvarranges[[idxchr]]$case_ID %in% samplemetatablewithentity$ID, ]
	    listvarinregion = unique(splitByOverlaptolist(testgns1gene, variantsassociated, "SYMBOL"))
	    topregions <- sigDriver:::getregionTopMutatedRanges(testgns1gene,
	                                            variantsassociated[listvarinregion,], #somaticvarranges[[idxchr]],
	                                            samplemetatablewithentity$ID,
	                                            samplemetatablewithentity,
	                                            pctin,restricted,useSigWeight=sigweight)
			if (dim(topregions)[1] != 0){
	      print(paste("extract variants for top:", dim(topregions)[1]))
		    #get what samples are in the bin       
		    topsomaticvrangesdf = cbind(as.character(seqnames(testgns1gene)),topregions$START,topregions$END)
		    topsomaticvrangesdf = data.frame(topsomaticvrangesdf,stringsAsFactors=FALSE)
		    colnames(topsomaticvrangesdf)  = c("CHROM","POS","END")                           
				topsomaticvranges = makeGRangesFromDataFrame(topsomaticvrangesdf,seqnames.field="CHROM",start.field="POS",end.field="END",keep.extra.columns=TRUE)
		    listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, variantsassociated, "SYMBOL"))
		    varframe = data.frame( variantsassociated[listvarinregion,])
		    #print(varframe$case_ID)
		    #print(head(samplemetatablewithentity[varframe$case_ID,sigNames]))
		    #calculate prevalence of signatures in a bin
		   	#print(head(samplemetatablewithentity))
		    #print(head(varframe))
		    
		    sample_sig_table_pos5 = colSums(samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),sigNames]>0.05) / length(varframe$case_ID)
				
				sampleEntityOccurence = samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),]$entity
				
				sampleSignatureSigExpAvg = samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),sigNames]
				sampleSignatureSigExpAvg[sampleSignatureSigExpAvg == 0] = NA
				sigExpPositivityInSamples = colMeans(sampleSignatureSigExpAvg,na.rm=T)
				
		    #print(sample_sig_table_pos5)#samplemetatablewithentity[which(samplemetatablewithentity$ID %in% varframe$case_ID),sigNames]>0.05)
		    return(list(sample_sig_table_pos5,sigExpPositivityInSamples,sampleEntityOccurence))
	    }else{
	    	return(NULL)
	    }
	  },
  error=function(cond) {
    print(paste("Error:",cond))
    return(NULL)
    #return(cond)
  })
}

