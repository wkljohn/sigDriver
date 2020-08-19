

getWindowNVarWithWeight <- function(startcoord,endcoord,framesize,variantRanges){
  
  #bin by matrix elements per generated  interval
  bincoord = seq(startcoord,endcoord,framesize)
  variantRangesdt = variantRanges[,bin:=findInterval(variantRanges$start, bincoord)]
  chrbins <- variantRangesdt[,list(start=bincoord[bin],end=bincoord[bin]+framesize,.N,sum(total_variants)), by=bin]
  return(chrbins)
}

#get max 5% mutated region for given genomic grange
getregionTopMutatedRanges <- function(gregion,variantRanges,tumorsincluded,samplemetatablewithentity,topProportion = 0.05,restricted=FALSE,minvariantinframe = 4,minimumframeinclude = 2,isCNV=FALSE){
  if (restricted){
    if (isCNV) {
      framesize = 10000
      framesizeweightsize = 0 #no weighting frame allowed
      SDsmaxthres = 4
      divfrompeak = 0.99
    }else{
      #refineframesize=5
      framesize = 30
      framesizeweightsize = 12
      
      
      #minvariantinframe = 4
      #minimumframeinclude = 2
      
      SDsmaxthres = 4
      divfrompeak = 0.5
    }
  }else{	      
    framesize = 200
    framesizeweightsize = 12
    minvariantinframe = 5
    SDsmaxthres = 4
  }
  print("INIT")
  nvarmedian=median(samplemetatablewithentity$total_variants)
  #add nvar column to variant ranges
  
  
  #extract tumors included only
  variantRanges = variantRanges[variantRanges$case_ID %in% tumorsincluded]
  variantRanges = sort(variantRanges,ignore.strand=TRUE)
  variantRanges = as.data.table(variantRanges)
  rownames(samplemetatablewithentity) = samplemetatablewithentity$ID 
  variantRanges$total_variants = log2(nvarmedian + 1)/log2(samplemetatablewithentity[as.character(variantRanges$case_ID),]$total_variants + 1)
  #variantRanges$total_variants = (samplemetatablewithentity[as.character(variantRanges$case_ID),]$total_variants - nvarmedian) / sd(samplemetatablewithentity[as.character(variantRanges$case_ID),]$total_variants)
  #variantRanges$total_variants[which(variantRanges$total_variants >= 0)] = 1 / (variantRanges$total_variants[which(variantRanges$total_variants >= 0)] + 1)
  #variantRanges$total_variants[which(variantRanges$total_variants < 0)] = 2^(-variantRanges[variantRanges$total_variants < 0,]$total_variants + 1)
  #variantRanges$total_variants = (nvarmedian)/(samplemetatablewithentity[as.character(variantRanges$case_ID),]$total_variants)
  
  #mutationfreqlist = list()
  startcoord = start(gregion)
  endcoord = end(gregion)
  print("S1")
  mutationfreqlist=getWindowNVarWithWeight(startcoord,endcoord,framesize,variantRanges)
  #50 frame
  startcoord = start(gregion) + framesize / 2
  endcoord = endcoord - (framesize / 2)
  mutationfreqlist=rbind(mutationfreqlist,getWindowNVarWithWeight(startcoord,endcoord,framesize,variantRanges))
  
  print("S2")
  #print((mutationfreqlist))
  mutationfreqlistframe=data.frame(mutationfreqlist,stringsAsFactors=FALSE)
  colnames(mutationfreqlistframe) = c("BIN","START","END","COUNT","WEIGHTED")
 # colnames(mutationfreqlistframe) = c("START","END","COUNT","WEIGHTED")
  #print(mutationfreqlistframe)
  print(sum(mutationfreqlistframe$COUNT))
  #print(quantile(mutationfreqlistframe$COUNT,0.9))
  #Exclude zeros for calculation of stats
  mutationfreqlistframe = mutationfreqlistframe[which(mutationfreqlistframe$COUNT > 0),]
  mutationfreqlistframe = mutationfreqlistframe[which(!is.na(mutationfreqlistframe$START)),]
  
  totallen = dim(mutationfreqlistframe)[1]
  medianweight = median(mutationfreqlistframe$WEIGHTED)
  sdmut = sd(mutationfreqlistframe$WEIGHTED)
  maxmut = max(mutationfreqlistframe$WEIGHTED)
  maxmutcoord = mutationfreqlistframe$START[which(mutationfreqlistframe$WEIGHTED == maxmut)]
  print(medianweight)
  print(sdmut * SDsmaxthres)
  #print(mutationfreqlistframe)
  mutationfreqlistframe = mutationfreqlistframe[which(mutationfreqlistframe$COUNT >= minvariantinframe),]
  
  #print(mutationfreqlistframe[mutationfreqlistframe$WEIGHTED > 5,])
  print(head(mutationfreqlistframe[order(mutationfreqlistframe$WEIGHTED,decreasing = T),],n=10))
  
  if (!is.null(dim(mutationfreqlistframe)[1]) && dim(mutationfreqlistframe)[1] > 0){
    mutationfreqlistframe <- mutationfreqlistframe[order(mutationfreqlistframe$COUNT,decreasing = TRUE),] 
    mutationfreqlistframe$REGIONWEIGHTED = FALSE
    
    
    #enrich all the peaks which are not neighbours
    curridx = 1
    while (curridx != -1){
      if (dim(mutationfreqlistframe)[1] >= curridx &&
          mutationfreqlistframe[curridx,]$WEIGHTED > minvariantinframe && 
          mutationfreqlistframe[curridx,]$WEIGHTED > (maxmut * divfrompeak)){
        #mutationfreqlistframe[curridx,]$WEIGHTED > (maxmut / 2)){
        #add to other peaks their neighbours
        if (curridx == 1 || abs(mutationfreqlistframe[curridx,]$START - mutationfreqlistframe[(curridx-1),]$START) > 500){
          nowstart = mutationfreqlistframe[curridx,]$START
          mutationfreqlistframe$REGIONWEIGHTED[curridx] = TRUE
          print(paste("weight",nowstart))
          nowweightranges = which(abs(mutationfreqlistframe$START - nowstart) <= (framesize * framesizeweightsize) & mutationfreqlistframe$REGIONWEIGHTED == FALSE)
          print(paste("on range:",nowweightranges))
          print(paste("org:",mutationfreqlistframe[nowweightranges,]$WEIGHTED))
          mutationfreqlistframe$WEIGHTED[nowweightranges] = mutationfreqlistframe$WEIGHTED[nowweightranges] * 2
          mutationfreqlistframe$REGIONWEIGHTED[nowweightranges] = TRUE
        }
        curridx = curridx + 1
      }else{
        curridx = -1  #stop at <peak/2 or when too small
      }
    }
    
    print(paste("thres cutoff:",maxmut*divfrompeak,restricted,sdmut,SDsmaxthres,medianweight,medianweight+ sdmut * SDsmaxthres))
    #restricted frame vs jumbo frame : adopt to the expected motif size of Mutational signatures
    if (restricted){
      if (dim(mutationfreqlistframe[which( mutationfreqlistframe$WEIGHTED > (maxmut*divfrompeak) & 
                                           mutationfreqlistframe$WEIGHTED > (sdmut * SDsmaxthres)),])[1] > minimumframeinclude){
        mutationfreqlistframe = mutationfreqlistframe[which( mutationfreqlistframe$WEIGHTED > (maxmut*divfrompeak) & 
                                                               mutationfreqlistframe$WEIGHTED > (sdmut * SDsmaxthres)),]
      }else if (dim(mutationfreqlistframe)[1] >= minimumframeinclude){
        mutationfreqlistframe = mutationfreqlistframe[c(1:minimumframeinclude),]

      }
      
      #refine windows further if possible
      # USEextraSharpPeak = FALSE
      # for (i in 1:dim(mutationfreqlistframe)[1]){
      #   mutationfreqlist=getWindowNVarWithWeight(mutationfreqlistframe$START[i],mutationfreqlistframe$END[i],refineframesize,variantRanges)
      #   mutationfreqlist=c(mutationfreqlist,getWindowNVarWithWeight(mutationfreqlistframe$START[i] + as.integer(refineframesize/2),mutationfreqlistframe$END[i],refineframesize,variantRanges))
      # 
      #   mutationfreqlistframe2=do.call(rbind, mutationfreqlist)
      #   mutationfreqlistframe2=data.frame(mutationfreqlistframe2,stringsAsFactors=FALSE)
      #   mutationfreqlistframe2 <- mutationfreqlistframe2[order(mutationfreqlistframe2[,4],decreasing = TRUE),]
      #   colnames(mutationfreqlistframe2) = c("START","END","COUNT","WEIGHTED")
      #   sumvar = sum(mutationfreqlistframe2$COUNT) / 2
      # 
      #   #print(paste("frame",i,"totalvar",sumvar))
      #   #print(mutationfreqlistframe2)
      #   if (mutationfreqlistframe2$COUNT[1] > sumvar * 0.8 || USEextraSharpPeak) {
      #     USEextraSharpPeak = TRUE
      #     #very condensed peak type
      #     print(paste("replace",i,"totalvar",sumvar))
      #     mutationfreqlistframe$START[i] = mutationfreqlistframe2$START[1]
      #     mutationfreqlistframe$END[i] = mutationfreqlistframe2$END[1]
      #     #do not replace weight as they might be filtered
      #   }
      # }
    }else{

      if (dim(mutationfreqlistframe[mutationfreqlistframe$WEIGHTED > 
                                    medianweight+(sdmut * SDsmaxthres),])[1] > ceiling(totallen * topProportion)){
        mutationfreqlistframe = mutationfreqlistframe[mutationfreqlistframe$WEIGHTED > medianweight+(sdmut * SDsmaxthres),]
      }else if (dim(mutationfreqlistframe)[1] >= totallen * topProportion){
        mutationfreqlistframe = mutationfreqlistframe[c(1:ceiling(totallen * topProportion)),]
      }
      
      #peaks only
      mutationfreqlistframe = mutationfreqlistframe[mutationfreqlistframe$WEIGHTED >= (medianweight+(sdmut * SDsmaxthres)),]
      
    }
  }
  
  
  return(mutationfreqlistframe)
  #PENDING, toppct filter
}
