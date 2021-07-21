signatureRepresentationAdjustment <- function(gns,
                             signature_test,
                             sigexpinfo,
                             backgroundsigs,
                             somaticvarranges,
                             samplemetatablewithentity,
                             threads){
   print("Calculating signature-bin distribution")
		cl <- parallel::makeCluster(threads,useXDR=TRUE)#, outfile='/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/logs/bin_bias_info_parallel.log')
		binStatsList=parSapply(cl, 1:length(gns$SYMBOL), getBinSignatureRepresentation, gns=gns,somaticvarranges=somaticvarranges,samplemetatablewithentity=sampleinfofiltered,sigweight=FALSE,sigNames=rownames(sigexpinfo))
		stopCluster(cl)
		binStatsDF = do.call(rbind,binStatsList)
		binStatsMeans = data.frame(colMeans(binStatsDF),stringsAsFactors=F)
		binSignatureWeights = data.frame((1-binStatsMeans)^2)
		#get list of backgrounds, do not weight on background
		backgroundsigslist = strsplit(gsub("\\s","",backgroundsigs),",")[[1]]
		binSignatureWeights[which(rownames(binSignatureWeights) %in% c(signature_test,backgroundsigslist)),] = 0
		
		#test binning weighting results
		#compute weight on each variant
		sigweight=binSignatureWeights
		#samplemetatablewithentity=sampleinfofiltered
		
		print("Computing variant-wise weight")
		sigweight=binSignatureWeights
		backgroundsigsidx = which(rownames(sigexpinfo) %in%  c(signature_test,backgroundsigslist))
		#preprocessing on matricies
		samplemetatablewithentity = sampleinfofiltered
		sigSampleInfoMatrix = samplemetatablewithentity[,rownames(sigexpinfo)]
		rownames(sigSampleInfoMatrix) =  sampleinfofiltered$ID
		#signature positivity matrix
		sigSampleInfoMatrix[sigSampleInfoMatrix>0.05] = 1
		sigSampleInfoMatrix[sigSampleInfoMatrix<=0.05] = 0
		
		cl <- parallel::makeCluster(threads,useXDR=TRUE)#, outfile='/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/results/logs/bin_bias_info_parallel.log')
		somaticvarranges=parSapply(cl, 1:length(somaticvarranges), signatureRepresentationWeight, somaticvarranges=somaticvarranges,sigweight=sigweight,backgroundsigsidx=backgroundsigsidx,sigSampleInfoMatrix=sigSampleInfoMatrix,sigexpinfo=sigexpinfo)
		stopCluster(cl)
		
		return (somaticvarranges)
}


signatureRepresentationWeight <- function(i,somaticvarranges,sigweight,
																			backgroundsigsidx,
																			sigSampleInfoMatrix,
																			sigexpinfo){
		require(dplyr)
		print(paste("weighting chr",i))
		somaticvarrangesDF = as.data.frame(somaticvarranges[[i]])
		somaticvarrangesDF$WEIGHT = 1
		#somaticvarranges[[i]]$WEIGHT = 1
		pos_unique_sites = unique(start(somaticvarranges[[i]]))
		
		#add exposures
		somaticvarrangesDF = merge(somaticvarrangesDF, sigSampleInfoMatrix, by.x="case_ID", by.y="row.names")
		vrange_summary_site = somaticvarrangesDF %>%
	  group_by(start) %>%
	  summarise_at(vars(rownames(sigexpinfo)), mean)
	  vrange_summary_site$max_exp = do.call(pmax,  vrange_summary_site[,rownames(sigexpinfo)])
	  #vrange_summary_site$isTestsigMax = vrange_summary_site[,signature_test] >= vrange_summary_site$max_exp
	  
	  vrange_summary_site[,rownames(sigexpinfo)] = (1-vrange_summary_site[,rownames(sigexpinfo)]) + (t(t(vrange_summary_site[,rownames(sigexpinfo)]) * sigweight[,1]))
	  # (1-vrange_summary_site[756,rownames(sigexpinfo)]) + (vrange_summary_site[756,rownames(sigexpinfo)] * sigweight[,1])
	  # vrange_summary_site[vrange_summary_site$start == 191784318,]
	  vrange_summary_site[vrange_summary_site<0.0001] = 1
	  vrange_summary_site[,backgroundsigsidx+1] = 1
	  vrange_summary_site$min =   do.call(pmin,  vrange_summary_site[,rownames(sigexpinfo)])
	  #vrange_summary_site[vrange_summary_site$isTestsigMax,]$min = 1
	  
	  vrange_summary_site = data.frame(vrange_summary_site,stringsAsFactors=F,row.names=vrange_summary_site$start)
	 
	  somaticvarranges[[i]]$WEIGHT = vrange_summary_site[as.character(start(somaticvarranges[[i]])),]$min
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
  })
}
