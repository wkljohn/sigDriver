#cat EX1KBrun.cmd | sed 's/ 10 / 30 /g' | sed 's/merged.avout.hg19_multianno.nofilt.txt.gz/merged.avout.hg19_multianno.txt.gz/g'
#cat SBS_qsub.cmd | sed 's/ 10 / 30 /g' | sed 's/merged.avout.hg19_multianno.nofilt.txt.gz/merged.avout.hg19_multianno.txt.gz/g'

doassocandwriteSKAThotspot <- function (igene,
                             bigtabledesc,gns,somaticvarranges,outfile,
                             samplemetatablewithentity,sigtest,pathfile,varianttype,
                             verbose=0){
  #require(bigmemory)
  require(GenomicRanges)
  require(dplyr)
  require(reshape2)
  require(pscl)
  require(data.table)
  require(SKAT)
  require(matrixStats)
  require(sigDriver)
  
  
  
  splitByOverlaptolist <-
    function(query, subject, column="ENTREZID", ...)
    {
      olaps <- findOverlaps(query, subject, ...)
      return(subjectHits(olaps))
      #f1 <- factor(subjectHits(olaps),
      #             levels=seq_len(subjectLength(olaps)))
      #splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
    }
  
  
  out <- tryCatch({
	  lookuparray = data.frame(cbind(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
	                                 c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)),stringsAsFactors = FALSE)
	  rownames(lookuparray) = lookuparray$X1
	  lookuparray$X2 = as.numeric(lookuparray$X2)
	  testedRegionCoordinates = ""
	  
	  
	  #gene based test
	  if (length(pathfile) == 1){
	  	#igene=which(gns$SYMBOL=="RAD51B")
	    #igene=which(gns$SYMBOL=="ADGRG6")
	    restricted = TRUE
	    pctin = 0.005
	    genetest = gns$SYMBOL[igene]
	    
	    testgns1gene = gns[which(gns$SYMBOL == genetest)]
	    idxchr = lookuparray[as.character(seqnames(testgns1gene)[1]),]$X2
	    variantsassociated = somaticvarranges[[idxchr]][somaticvarranges[[idxchr]]$case_ID %in% samplemetatablewithentity$ID, ]
	    listvarinregion = unique(splitByOverlaptolist(testgns1gene, variantsassociated, "SYMBOL"))
	    

	    #direct GT table or most mutated region(nc + coding)
	    if (varianttype == 10 || varianttype == 11 || varianttype == 20 || varianttype == 50){
	      topregions <- getregionTopMutatedRanges(testgns1gene,
	                                              variantsassociated[listvarinregion,], #somaticvarranges[[idxchr]],
	                                              samplemetatablewithentity$ID,
	                                              samplemetatablewithentity,
	                                              pctin,restricted)
	      #toc()
	      #NOT ANY MORE on 10kb windows:50 only should contain 1 window
	      #if (varianttype == 50) { topregions = topregions[1,]}
	      if (verbose==1) print(topregions)
	      if (dim(topregions)[1] != 0){
	        print("extract variants for top")
	        #topregions[order(topregions[,3],decreasing = TRUE),]
	        topsomaticvrangesdf = cbind(as.character(seqnames(testgns1gene)),topregions$START,topregions$END)
	        topsomaticvrangesdf = data.frame(topsomaticvrangesdf,stringsAsFactors=FALSE)
	        colnames(topsomaticvrangesdf)  = c("CHROM","POS","END")
	        
	        topsomaticvranges = makeGRangesFromDataFrame(topsomaticvrangesdf,seqnames.field="CHROM",start.field="POS",end.field="END",keep.extra.columns=TRUE)
	        testedRegionCoordinates = paste(t(cbind(as.character(seqnames(topsomaticvranges)),start(topsomaticvranges),end(topsomaticvranges))),sep=":",collapse=":")

	        #listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, somaticvarranges[[idxchr]], "SYMBOL"))
	        listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, variantsassociated, "SYMBOL"))
	        if (verbose==1) print(topsomaticvranges)
	        if (verbose==1) print(length(listvarinregion))
	      }else{
	        listvarinregion = c()
	      }
	    #CNV type handing
	    }else  if (varianttype == 51){
	      #igene = 1664
	      #igene = 1273
	      #source("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/scripts/R_association/hotspot_routine.R")
	      topregions <- getregionTopMutatedRanges(testgns1gene,
	                                              variantsassociated[listvarinregion,], #somaticvarranges[[idxchr]],
	                                              samplemetatablewithentity$ID,
	                                              samplemetatablewithentity,
	                                              pctin,restricted,isCNV=TRUE)
	      #NOT ANY MORE on 10kb windows:50 only should contain 1 window
	      #if (varianttype == 50) { topregions = topregions[1,]}
	      if (verbose==1) print(topregions)
	      if (dim(topregions)[1] != 0){
	        if (verbose==1) print("extract variants for top")
	        #topregions[order(topregions[,3],decreasing = TRUE),]
	        topsomaticvrangesdf = cbind(as.character(seqnames(testgns1gene)),topregions$START,topregions$END)
	        topsomaticvrangesdf = data.frame(topsomaticvrangesdf,stringsAsFactors=FALSE)
	        colnames(topsomaticvrangesdf)  = c("CHROM","POS","END")
	        
	        topsomaticvranges = makeGRangesFromDataFrame(topsomaticvrangesdf,seqnames.field="CHROM",start.field="POS",end.field="END",keep.extra.columns=TRUE)
	        testedRegionCoordinates = paste(t(cbind(as.character(seqnames(topsomaticvranges)),start(topsomaticvranges),end(topsomaticvranges))),sep=":",collapse=":")
	        
	        #listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, somaticvarranges[[idxchr]], "SYMBOL"))
	        listvarinregion = unique(splitByOverlaptolist(topsomaticvranges, variantsassociated, "SYMBOL"))
	        if (verbose==1) print(topsomaticvranges)
	        if (verbose==1) print(length(listvarinregion))
	      }else{
	        listvarinregion = c()
	      }
	    }else{-
	      testedRegionCoordinates = paste(t(cbind(as.character(seqnames(testgns1gene)),start(testgns1gene),end(testgns1gene))),sep=":",collapse=":")
	      
	    }
	    
	    
	    #split into information tables
	   # library("tictoc")
	   # tic("DT")
	    varframe = data.frame( variantsassociated[listvarinregion,])
	    coordinatesframe = varframe[,c(1,2,3)] %>% distinct()
			
			#remove redundant variants if variants are 1) singletons 2) from the same case 3) >=2 variants per case
			max_case_variants = 2
			min_tested_variants = 6
	    if (dim(varframe)[1] >= min_tested_variants){
	      if (verbose==1) print("filt excess singleton 3")
				variants_per_case = as.data.table(table(varframe$case_ID))
				variants_per_case_idx = which(variants_per_case$N >= max_case_variants)
				if (length(variants_per_case_idx) > 0){
					variants_per_case_ids = variants_per_case[variants_per_case_idx,]$case_ID
					#print(variants_per_case_ids)
					variantstable = table(varframe$start)
					variantstable_singleton_idx = which(variantstable == 1)
					if (length(variantstable_singleton_idx) > 0){
						#print(names(variantstable)[variantstable_singleton_idx])
						variantstable_singletons = as.numeric(names(variantstable)[variantstable_singleton_idx])
						
						#process each excessive case
						for (x in 1:length(variants_per_case_idx)){
						  if (verbose==1) print(paste("Filt",variants_per_case[variants_per_case_idx[x],1,with=FALSE]))
							#print(as.character(varframe$case_ID))
							#print(variants_per_case[variants_per_case_idx[x],1,with=FALSE])
							case_remove = which(as.character(varframe$case_ID) %in% variants_per_case[variants_per_case_idx[x],1,with=FALSE])
							#print(case_remove)
							case_singletons = variantstable_singletons[which(variantstable_singletons %in% varframe$start[case_remove])]
							nvar_remove = as.numeric(variants_per_case[variants_per_case_idx[x],2,with=FALSE]) - 2
							#print(case_singletons)
							varframe = varframe[!(varframe$start %in% case_singletons[1:nvar_remove]),]
							#print(paste("removing",case_singletons[1:nvar_remove]))

						}
					}
					
				}
			}

	    
	    #case-genotype matrix, at least two variants
	    if (dim(varframe)[1] >= min_tested_variants){
	      #print("prepare var frame by case")
	      #print(head(varframe))
	      varframebycase = dcast(varframe,start~case_ID,fun.aggregate = function(x){return(as.integer(length(x)));})
	      #print(varframebycase)
	      #print("done var frame by case")
	      varsamplesabsent = samplemetatablewithentity$ID[!(samplemetatablewithentity$ID %in% colnames(varframebycase))]
	      varsamplespresent = colnames(varframebycase)
	      
	      #append non-positive cases with positive cases
	      zeromatrix = matrix(0, nrow = dim(varframebycase)[1], ncol = length(varsamplesabsent) )
	      colnames(zeromatrix) = varsamplesabsent
	      varframebycase = cbind(varframebycase,zeromatrix)
	      
	      rownames(varframebycase) = varframebycase$start
	      varframebycase = as.matrix(varframebycase[,samplemetatablewithentity$ID])
	    
	      #print("compute summary")
	      #tumor load X varframe
	      varframebycasewithtumorload = t(as.matrix(t(varframebycase) * (1 / (samplemetatablewithentity$total_variants))))
	      varframebycasewithtumorload[varframebycasewithtumorload==0] = NA

	      varframebycasewithtumorloadmedians = rowMedians((varframebycasewithtumorload),na.rm = T)

	      weightframenc <- data.frame(cbind(varframe$start,rep(1,dim(varframe)[1])),stringsAsFactors = F)
	      colnames(weightframenc) = c("start","weight")
	      
	      #print("compute weight")
	      weightframencLIST = (as.vector(rowSums(varframebycase))) * (varframebycasewithtumorloadmedians  )
	      if (verbose==1) print(weightframencLIST)
	      
	      gttable = as.matrix(table(variantsassociated[listvarinregion,]$case_ID))

	      samplemetatablewithentitygt=(samplemetatablewithentity)
	      samplemetatablewithentitygt$gt = merge(x=samplemetatablewithentitygt,y=gttable,all.x=TRUE,all.y=FALSE,by.x="ID",by.y="row.names")$V1
	      samplemetatablewithentitygt$gt[is.na(samplemetatablewithentitygt$tmp)] = 0
	      samplemetatablewithentitygt$gt[is.na(samplemetatablewithentitygt$gt)] = 0
	    }else{
	      varframebycase = 0
	    }
	    
	  }else{
		  #NOT IMPLEMENTED
	    
	  }
	  
	  #head(samplemetatablewithentitygt)
	  #print(sum(varframebycase))
    if ( sum(varframebycase,na.rm = T) >= 6){
      #Calculate entity information
      if("MASS" %in% (.packages())){detach("package:MASS",character.only = TRUE)}
      #print("DOING SPLIT")
      
      tmpdcast = data.frame(samplemetatablewithentitygt) %>% dplyr::select("entity","gt")
      tmpdcast$gt[tmpdcast$gt>1] = 1
      tmpdcast=dcast(tmpdcast,entity~gt,fun.aggregate=sum)
      tmpdcast=tmpdcast[,c(1,3)]
      tmpdcast$entity=as.character(tmpdcast$entity)
      tmpdcast$entity=gsub(" ","",tmpdcast$entity)
      #Entity count correction does not work
      #fullentitygttable = tmpdcast
      #rownames(fullentitygttable) = fullentitygttable$entity
      #samplemetatablewithentity$entityMutationCnt = fullentitygttable[as.character(samplemetatablewithentity$entity),2]
      tmpdcast=head(tmpdcast[order(tmpdcast[,2],decreasing=TRUE),],n=5)
      
      
      
      #the y ~ confounders declaration
      #log2(samplemetatablewithentity$sigvar + 1)
      #samplemetatablewithentity$normalized_exposures
      #samplemetatablewithentity$rank   ,produce bad p for background EX1,SBS1
      #
      if (verbose==1) print("DOING SKAT")
			#obj<-SKAT_Null_Model( samplemetatablewithentity$normalized_exposures ~ samplemetatablewithentity$entity + samplemetatablewithentity$total_variants + log2(samplemetatablewithentity$total_variants + 1)+ samplemetatablewithentity$gender, out_type="C")#, type.Resampling="bootstrap", Adjustment=TRUE,n.Resampling=10000)
      
      #handle sex chromosomes
      samplemetatablewithentity$gender = as.factor(samplemetatablewithentity$gender)
      samplemetatablewithentity$gender = 3- as.numeric(samplemetatablewithentity$gender)
      if (verbose==1) print(table(samplemetatablewithentity$gender))
      #if (idxchr == 23){
      #  attach(samplemetatablewithentity)
      #  obj<-SKAT_Null_Model_ChrX( normsig ~ Entity + nvar + log2(nvar + 1)+ gender, SexVar="gender",out_type="C")
      #}else{
      if (length(unique(samplemetatablewithentity$entity)) == 1){
        obj<-SKAT_Null_Model( samplemetatablewithentity$sigRank ~ samplemetatablewithentity$total_variants + log2(samplemetatablewithentity$total_variants + 1)+ samplemetatablewithentity$gender, out_type="C",n.Resampling=1000)
      }else{
        obj<-SKAT_Null_Model( samplemetatablewithentity$sigRank ~ samplemetatablewithentity$entity + samplemetatablewithentity$total_variants + log2(samplemetatablewithentity$total_variants + 1)+ samplemetatablewithentity$gender, out_type="C",n.Resampling=1000)
      }
        #}
      
			#samplemetatablewithentity$entitsiglevel 
			#samplemetatablewithentity$entity
			
			 obj$id_include = samplemetatablewithentity$ID
			varframebycasemat = as.matrix(t(varframebycase))
			#SKAT test on the genotype matrix given outcome-confounders obj
			#print("DOING SKATO")
			print("DOING SKAT multi rho")
			
			#DEBUGGER
			if (verbose==1){ 
			  print("weights")
			  print(head(weightframencLIST,n = 50))
			  print("varframe")
			  print(head(varframebycasemat[head(order(rowSums(varframebycasemat),decreasing = T)),]))
			}
			
			out = SKAT(varframebycasemat, obj, weights=weightframencLIST, method="SKATO")
			#out = SKAT(varframebycasemat, obj, weights=weightframencLIST, method="Burden")
			#out = SKAT(varframebycasemat, obj, weights=weightframencLIST, method="SKAT", r.corr=seq(0,1,0.1))
			  
			if (verbose==1 && !is.null(out$p.value.resampling)){
			  print(min(out$p.value.resampling))
			  print(median(out$p.value.resampling))
			  print(max(out$p.value.resampling))
			  print(Get_Resampling_Pvalue(out))
			}
			
			if (verbose==1) print(out)
			
			if (is.na(out$Q)){
			  qvaluescast=paste(out$param$q.val.each,collapse=",")
			  pvaluescast=paste(out$param$p.val.each,collapse=",")
			  rtnvar=c(genetest,testedRegionCoordinates,out$param$n.marker.test,sum(varframebycase),dim(samplemetatablewithentity)[1],qvaluescast,out$p.value,pvaluescast,paste(t(tmpdcast),collapse=" "))
			  rtnline=rtnline=paste(rtnvar,collapse=" ",sep=" ")
			}else{
			  rtnvar=c(genetest,testedRegionCoordinates,out$param$n.marker.test,sum(varframebycase),dim(samplemetatablewithentity)[1],out$Q,out$p.value,"NA",paste(t(tmpdcast),collapse=" "))
			  rtnline=paste(rtnvar,collapse=" ",sep=" ")
			}
			
      #print(paste(genetest,coef(summary(g))[2,1],coef(summary(g))[2,4],adjsandwitch[2,1],paste(adjsandwitch[,4],collapse = " "),sep=" "))
      #write.table(paste(genetest,coef(summary(g))[2,1],coef(summary(g))[2,4],adjsandwitch[2,1],paste(adjsandwitch[,4],collapse = " "),sep=" "),outfile,sep="\t",append = T,col.names = NA,quote = F)
			if (verbose==1) print(rtnline)
      write(rtnline,file=outfile,append=TRUE)
      return(rtnvar)
    }else{
      rtnvar=c(genetest,"NA","NA","NA")
      rtnline=paste(rtnvar,collapse=" ",sep=" ")
      if (verbose==1) print(rtnline)
      write(rtnline,file=outfile,append=TRUE)
      return(rtnvar)
      
    }
   
  },
  error=function(cond) {
    rtnvar=c(genetest,"ERR","ERR","ERR")
    rtnline=paste(rtnvar,collapse=" ",sep=" ")
    print(rtnline)
    write(rtnline,file=outfile,append=TRUE)
    return(rtnvar)
  })
}
