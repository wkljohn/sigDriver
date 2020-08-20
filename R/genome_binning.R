getChrCut = function(somaticvarranges){
  #binsize=100000#20000
  binsize=2000
  #overlapsize=1000#2000
  overlapsize=binsize/2#2000
  minvarperwindow=6 #minimum is the hard cut off in the end
  chromosomearray = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  gnsnew = c()
  fixvcfinfo=as.data.table(somaticvarranges)
  for (i in 1:length(chromosomearray)){
    print(paste("generating variant bins:",chromosomearray[i]))
    perchrdtR = fixvcfinfo[seqnames==chromosomearray[i],]
    #perchrdt = perchrdt[1:100000,]
    #print(paste(i,min(perchrdtR$start)))
    chrmin = 1#min(perchrdtR$POS)
    chrmax = max(perchrdtR$end)
    
    #bin by matrix elements by first window offset
    bincoord = seq(chrmin,chrmax,binsize)
    perchrdt = perchrdtR[,bin:=findInterval(perchrdtR$start, bincoord)]
    
    #filtering
    chrbins <- perchrdt[,list(seqnames=chromosomearray[i],start=bincoord[bin],end=bincoord[bin]+binsize,.N), by=bin]
    #chrbins = chrbins[which(chrbins$N >= minvarperwindow),c("seqnames","start","end")]
    #chrbins$SYMBOL = paste(chrbins$seqnames,":",chrbins$start,"-",chrbins$end,sep="")
    
    #half overlaping bin
    #bin by matrix elements per generated 25bp interval
    bincoord = seq(chrmin+overlapsize,chrmax-overlapsize,binsize)
    perchrdt = perchrdtR[,bin:=findInterval(perchrdtR$start, bincoord)]
    chrbins <- rbind(chrbins,perchrdt[,list(seqnames=chromosomearray[i],start=bincoord[bin],end=bincoord[bin]+binsize,.N), by=bin])
    
    #filtering
    #chrbins = rbind(chrbins,perchrdt[,list(seqnames=chromosomearray[i],start=bincoord[bin],end=bincoord[bin]+binsize,.N), by=bin])
    #print(median(chrbins$N))
    chrbins = chrbins[which(chrbins$N >= minvarperwindow),c("seqnames","start","end")]
    chrbins$SYMBOL = paste(chrbins$seqnames,":",chrbins$start,"-",chrbins$end,sep="")
    chrbins = chrbins[order(chrbins$start),]
    
    gnsnew = rbind(gnsnew,chrbins)
  }
  
  gnsnew = gnsnew[!is.na(gnsnew$start),]
  
  return(gnsnew)
}

getChrCutNonOverlapping = function(){
  binsize=25
  chromosomearray = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  gnsnew = c()
  for (i in 1:length(chromosomearray)){
    print(paste("generating variant bins:",chromosomearray[i]))
    perchrdt = fixvcfinfo[seqnames==chromosomearray[i],]
    #perchrdt = perchrdt[1:100000,]
    chrmin = min(perchrdt$POS)
    chrmax = max(perchrdt$POS)
    
    #bin by matrix elements per generated 25bp interval
    bincoord = seq(chrmin,chrmax,binsize)
    perchrdt = perchrdt[,bin:=findInterval(perchrdt$POS, bincoord)]
    
    #filtering
    chrbins <- perchrdt[,list(seqnames=chromosomearray[i],start=bincoord[bin],end=bincoord[bin]+binsize,.N), by=bin]
    chrbins = chrbins[which(chrbins$N >= 5),c("seqnames","start","end")]
    chrbins$SYMBOL = paste(chrbins$seqnames,":",chrbins$start,"-",chrbins$end,sep="")
    
    
    
    gnsnew = rbind(gnsnew,chrbins)
  }
  
  return(gnsnew)
}
