
readAnnotationGTF <- function(annotation_gtf){

	if (file.exists(annotation_gtf)){
		print("Reading GTF")
		gtfref = rtracklayer::import(annotation_gtf)
		gtfref$level = as.numeric(gtfref$level)
		return(gtfref)
		
	}else{
		
		print("Warning: gtf not found, skipping annotation")
		return(NA)
	}
}

plot_lolli <- function(resultsSKATanno,out_path){

	require(rtracklayer)
	require(trackViewer)
	resultsimportancelist = list()
	for (i in 1:length(resultsSKATanno)){
		siteimportance = data.frame(do.call(rbind,resultsSKATanno[[i]]),stringsAsFactors=F,check.names=F)
		colnames(siteimportance) = c("Region","test_site","subregions","remaining_sites","marker_snvs","n.samples","Q","p.value","importance","entityvar")
		
		siteimportance$p.value = as.numeric(siteimportance$p.value)
	  minimportance = min(siteimportance$p.value)
	  importantcasescnt = 0
	  siteimportance$p.value[which(siteimportance$p.value < 0)] = 0
	  
	  #get sites to use
	  sitestested = siteimportance$test_site
	  sitestested = sitestested[which(sitestested != "NA")]
	  if(length(sitestested) == 1){
		siteimportance$entityvar[1] = siteimportance$entityvar[2]
	  }
	  #FIXFIX
	  sitestesteddf = t(as.data.table(strsplit(sitestested,":")))
	  colnames(sitestesteddf) = c("chr","start")
	  chrarrayindex = as.numeric(sitestesteddf[1,1])
	  sitestesteddf = as.data.frame(sitestesteddf,stringsAsFactors=F,check.names=F)
	  sitestesteddf$chr = paste("chr",sitestesteddf$chr,sep="")
	  sitestesteddf$chr = gsub("chr23","chrX",sitestesteddf$chr)
	  sitestesteddf$chr = gsub("chr24","chrY",sitestesteddf$chr)
	  sitestestedGR = makeGRangesFromDataFrame(sitestesteddf,
	  			seqnames.field="chr",
	  			start.field="start",
	  			end.field="start",keep.extra.columns=TRUE)
	  			
	  
	  #get the variants
	  overlapvariants = findOverlaps(sitestestedGR,somaticvarranges[[chrarrayindex]])
	  snvsigtable = data.frame(somaticvarranges[[chrarrayindex]][subjectHits(overlapvariants)],stringsAsFactors=F)
	  
	  siteimportance$test_site = gsub(".*:","",siteimportance$test_site)
	  uniqsites = unique(snvsigtable[,c("seqnames","start","REF","ALT")])
	  sitestats = data.table(table(snvsigtable$start,snvsigtable$REF,snvsigtable$ALT))
	  uniqsites = merge(x=uniqsites,y=sitestats,by.x=c("start","REF","ALT"),by.y=c("V1","V2","V3"))
	  uniqsites = merge(x=uniqsites,y=siteimportance,by.x="start",by.y="test_site")
	  uniqsites$importance = as.numeric(uniqsites$importance)
	  minimportanceperturbed = min(uniqsites$p.value )
	  uniqsites$importance = -log10(uniqsites$importance)
	  
	  
	  #get site importance stats
	  siteimportancestats=c()
	  ##not useful anymore after perterb improvements
	  siteimportancestats$min = min(uniqsites$importance)

	  siteimportancestats$average = mean(uniqsites$importance)
	  siteimportancestats$max = as.numeric(max(uniqsites$importance))
	  siteimportancestats$sd = sd(uniqsites$importance)
	  
          #make entity info string
	  majorsiteentities = c()
	  countedpval = c()
	  uniqsites$EntityString = ""
	  uniqsites$entityvar = gsub("  "," ",uniqsites$entityvar)
	  for (h in 1:dim(uniqsites)[1]){
	    entitystring = ""
	    tmpsite1 = data.frame(t(matrix(nrow=2,strsplit(uniqsites$entityvar[h]," ")[[1]])),stringsAsFactors=F)
	    colnames(tmpsite1) = c("Entity","Count")
	    tmpsite1$Count = as.numeric(tmpsite1$Count)
	    if (grepl("-",tmpsite1$Entity)[1]){
	    	tmpsite1$Subtype = do.call(rbind,strsplit(tmpsite1$Entity,"-"))[,2]
	    }else{
	    	tmpsite1$Subtype = cbind(tmpsite1$Entity,"N")
	    }
	    
	    tmpsite1 = tmpsite1[which(tmpsite1$Count > 0),]
	    uniqsites$EntityString[h] = paste(substr(tmpsite1$Entity,1,2),substr(tmpsite1$Subtype,1,1),":",tmpsite1$Count,sep="",collapse=",")
	    
	    #importance by n-entities
	    if (uniqsites$importance[h] >= siteimportancestats$max * 0.33 && tmpsite1$Count[1] >= 1){
	      majorsiteentities = c(majorsiteentities,tmpsite1$Entity[1])
	    }
	    #most supported site recurrence
	    if (uniqsites$importance[h] >= siteimportancestats$max * 0.33 && !(uniqsites$p.value[h] %in% countedpval)){
	      countedpval = c(countedpval,uniqsites$p.value[h])
	      importantcasescnt = importantcasescnt + sum(tmpsite1$Count)
	    }
	  }
	  
	  sitescountonly = uniqsites[,c("start","N")]
	  sitescountonly = aggregate(sitescountonly$N~sitescountonly$start,FUN=sum)
	  
	  #filter singletons
	  uniqsitesplot=""
	  if (length(which(uniqsites$N<=2)) > 0){
	    sitescountonly = sitescountonly[which(sitescountonly[,2] == 1),]
	    excludestarts = uniqsites[which(uniqsites$importance < 1 & uniqsites$start %in% sitescountonly[,1]),]$start
	    uniqsitesplot = uniqsites[!(uniqsites$start %in% excludestarts),]
	  }else{
	    uniqsitesplot = uniqsites
          }
	  
	  #check output directories
	  if (!dir.exists(paste(out_path,"/plots/",sep=""))){
	  	dir.create(paste(out_path,"/plots/",sep=""), showWarnings = FALSE)
	  }
	  #check output directories
	  outfolderplots=paste(out_path,"/plots/",signature_test,"/",sep="")
	  if (!dir.exists(outfolderplots)){
	  	dir.create(outfolderplots, showWarnings = FALSE)
	  }


	  if (typeof(uniqsitesplot) == "list" && dim(uniqsitesplot)[1] > 0){
	    
	    sample.gr <- GRanges(uniqsitesplot$seqnames, IRanges( uniqsitesplot$start, width=1, 
	                                                          names=paste(uniqsitesplot$N,paste(uniqsitesplot$seqnames,uniqsitesplot$start,sep=":"),uniqsitesplot$EntityString,sep=" / ")))
	    sample.gr$score <- uniqsitesplot$importance
	    
	  
	    sample.gr$color <- as.character(paste(uniqsitesplot$REF,uniqsitesplot$ALT,sep="_"))
	    muttypes=c("C_A","C_G","C_T","T_A","T_C","T_G","OTHER")
	    sample.gr$color = gsub("A_C","T_G",sample.gr$color) #
	    sample.gr$color = gsub("A_G","T_C",sample.gr$color) #[sample.gr$color == "A_G"] = "T_C"
	    sample.gr$color = gsub("A_T","T_A",sample.gr$color) #[sample.gr$color == "A_T"] = "T_A"
	    sample.gr$color = gsub("G_A","C_T",sample.gr$color) #[sample.gr$color == "G_A"] = "C_T"
	    sample.gr$color = gsub("G_C","C_G",sample.gr$color) #[sample.gr$color == "G_A"] = "C_T"
	    sample.gr$color = gsub("G_T","C_A",sample.gr$color) #[sample.gr$color == "G_A"] = "C_T"
	    sample.gr$color[!(sample.gr$color %in% muttypes)] = "OTHER"
	    sample.gr$color <- factor(sample.gr$color, levels=muttypes)
	    legend <- as.integer(factor(muttypes,levels=muttypes))
	    names(legend) <-muttypes

	  }else{
	    print("Plot not possible")
	    dolollipopplot = F
	    sample.gr = c()
	  }
	  
	  #match the current results with annotation
	  annotationrecord = resultsimportancedf[which(resultsimportancedf$Region == siteimportance$Region[1] & resultsimportancedf$test_site == "ALL"),]
	  if (is.na(annotationrecord$annotation[1])){
	  	print(paste("No annotation for region",annotationrecord$Region[1]))
	  	next()
	  }
	  #search GTF
	  gtfmatch_transcript = gtfref[which(gtfref$gene_id == annotationrecord$annotationID[1] & gtfref$type == "transcript")]
	  gtfmatch_transcript = gtfmatch_transcript[order(gtfmatch_transcript$level),]
	  gtfmatch_exons = gtfref[which(gtfref$transcript_id == gtfmatch_transcript$transcript_id[1] & gtfref$type == "exon")]
	  #gtfmatch_transcript = gtfmatch_transcript[order(gtfmatch_transcript$transcript_type)]
	  
	  #get exons
	  exonstart = start(gtfmatch_exons)
	  exonend = end(gtfmatch_exons)
	  exonstrand = "+"
	  if (as.character(strand(gtfmatch_exons)[1]) == "+"){
	    print("+")
	    exonstrand = "+"
	    if (!is.na(exonstart)){
	      features <- GRanges(seqnames(gtfmatch_exons), IRanges(exonstart, 
	                                                          width=exonend-exonstart,
	                                                          names=paste0("exon", 1:length(exonstart))))
	      features$height <- list(rep(unit(12, "points"),length(exonend)))
	      
	      features=  c(GRanges(seqnames(gtfmatch_exons)[1], IRanges(exonstart[1]-1500, 
	                                                           width=1500,
	                                                           names="upstream")),features)
	      features$height <- as.list(c(unit(0.01,"points"),rep(unit(0.05, "points"),length(exonend))))
	      
	      
	      features$fill <- c("gray",rainbow(length(exonend)))
	    }else{
	#      names(testgns1gene)="unannotated gene"
	#      features <- testgns1gene
	#      features$height <- list(rep(unit(12, "points"),length(testgns1gene)))
	#      
	#      features$height <- as.list(rep(unit(0.05, "points"),length(testgns1gene)))
	#      features$fill <- c(rainbow(length(testgns1gene)))
	    }
	  }else{
	    
	    
	    if (!is.na(exonstart)){
	      print("-")
	      exonstrand = "-"
	      exonstart=rev(exonstart)
	      exonend=rev(exonend)
	      features <- GRanges(seqnames(gtfmatch_exons), IRanges(exonstart, 
	                                                          width=exonend-exonstart,
	                                                          names=paste0("exon", length(exonstart):1)))
	      features=  c(features,GRanges(seqnames(gtfmatch_exons)[1], IRanges(exonend[ length(exonend)], 
	                                                                    width=1500,
	                                                                    names="upstream")))
	      features$height <- as.list(c(rep(unit(0.05, "points"),length(exonend)),unit(0.01,"points")))
	      features$fill <- c(rainbow(length(exonend)),"gray")
	      
	    }else{
	#      names(testgns1gene)="unannotated gene"
	#      features <- testgns1gene
	#      features$height <- list(rep(unit(12, "points"),length(testgns1gene)))
	#      
	#      features$height <- as.list(rep(unit(0.05, "points"),length(testgns1gene)))
	#      features$fill <- c(rainbow(length(testgns1gene)))
	    }
	  }
	  
	  #features$fill <- brewer.pal(n = length(exonend), name = 'Spectral')
	  #do exonic annotation here
	  if (length(features) >= 2 && length(which(uniqsites$importance > 1 & uniqsites$importance > (max(uniqsites$importance) * 0.33 ))) > 0){
	    sitegranges <- GRanges(seqnames(gtfmatch_exons)[1], IRanges(uniqsites$start[which(uniqsites$importance > 1 & uniqsites$importance > (max(uniqsites$importance) * 0.33 ))],   width=1))
	    overlapcoding <- queryHits(findOverlaps(sitegranges,features))
	  }else{
	    overlapcoding <- c()
	  }
	  
	  if (length(somaticvarranges) > 1){
	    features$featureLayerID = "exons"
	   # topsomaticvranges$height = unit(0.03, "points")
	   # topsomaticvranges$fill = "red"
	   # topsomaticvranges$featureLayerID = "hotspot"
	    #names(topsomaticvranges) = rep("Top 1% hotspots",length(topsomaticvranges))
	    #features=  c(features,topsomaticvranges)
	  
	  	outdirplot = paste(outfolderplots,"/",signature_test,"_",annotationrecord$annotation[1],"_",gsub(":","-",annotationrecord$Region),".png",sep="")
	  	outdirplotsites = paste(outfolderplots,"/",signature_test,"_",annotationrecord$annotation[1],"_",gsub(":","-",annotationrecord$Region),".tsv",sep="")
	    print(paste("Plotting",outdirplot))
	   print(sample.gr)
	    png(outdirplot,width = 2000,height = 700)
	    lolliplot(sample.gr, features,ylab=paste(signature_test, " -log(p) change"), yaxis=c(0,ceiling(max(sample.gr$score))),legend=legend)
	    dev.off()
	      
	    outtablesnv = cbind(snvsigtable$seqnames,snvsigtable$start,snvsigtable$end,exonstrand)
	    
	    write(paste(annotationrecord$annotation[1],
	                annotationrecord$Region[1],
	                dim(uniqsites)[1],
	                minimportance,
	                minimportanceperturbed,
	                siteimportancestats$min,
	                siteimportancestats$max,
	                siteimportancestats$average,
	                siteimportancestats$sd,
	                length(unique(majorsiteentities)),
	                importantcasescnt,
	                length(overlapcoding)),
	          outdirplotsites,append = T)
	    
	    convertToComplement<-function(x){
	      bases=c("A","C","G","T")
	      xx<-unlist(strsplit(toupper(x),NULL))
	      paste(unlist(lapply(xx,function(bbb){
	        if(bbb=="A") compString<-"T"
	        if(bbb=="C") compString<-"G"
	        if(bbb=="G") compString<-"C"
	        if(bbb=="T") compString<-"A"
	        if(!bbb %in% bases) compString<-"N"
	        return(compString)
	      })),collapse="")
	    }
	#    
	#    #get motif per top hit
	#    uniqsites = uniqsites[order(uniqsites$importance,decreasing = T),]
	#    #totalsitesannotated = 0
	#    totalmotifswritten = 0
	#    lastsite = 0
	#    caseEntityCountingTable = c()
	#    for (k in 1:dim(uniqsites)[1]){
	#      #only get top 25bp window
	#      if (k == 1){
	#        topwindowstart = uniqsites$start[k]
	#      }
	#      
	#      #print(head(snvsigtable))
	#      tmpdcast = data.frame(table(snvsigtable[which(snvsigtable$start == uniqsites$start[k] & 
	#                                                      snvsigtable$REF == uniqsites$REF[k] & 
	#                                                      snvsigtable$ALT == uniqsites$ALT[k]),]$Entity),stringsAsFactors = F)
	#      casesfulltable = snvsigtable[which(snvsigtable$start == uniqsites$start[k] & 
	#                                                      snvsigtable$REF == uniqsites$REF[k] & 
	#                                                      snvsigtable$ALT == uniqsites$ALT[k]),c("samplesinvcf","Entity")]
	#      
	#      if (uniqsites$importance[k] > siteimportancestats$max / 3 && 
	#          uniqsites$importance[k] >= 1 ){
	#        #totalsitesannotated = totalsitesannotated + 1
	#        motifseqstart = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k]-20, end = uniqsites$start[k]-1))
	#        motifseqend = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k]+1, end = uniqsites$start[k]+20))
	#        motifseqsite = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k], end = uniqsites$start[k]))
	#        
	#        dowrite=1
	#        #if (motifseqsite == "G" || motifseqsite == "A"){
	#        #  #motifseqsitemut = gsub("-","",convertToComplement(uniqsites$ALT[k]))
	#        #  motifseqsitemut = gsub("-","",uniqsites$ALT[k])
	#        #}else 
	#        if (uniqsites$REF[k] == "-"){
	#          INSlength = length(uniqsites$ALT[k])
	#          
	#          motifseqsitemut = paste(motifseqsite,uniqsites$ALT[k],sep="")
	#        }else if (uniqsites$ALT[k] == "-"){
	#          dellength = nchar(uniqsites$REF[k])
	#          motifseqstart = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k]-20, end = uniqsites$start[k]-1))
	#          motifseqend = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k]+dellength, end = uniqsites$start[k]+20 + dellength - 1))
	#          motifseqsite = as.character(getSeq(Hsapiens, uniqsites$seqnames[k], start = uniqsites$start[k], end = uniqsites$start[k] + dellength - 1))
	#          
	#          motifseqsitemut = ""#paste(motifseqsite,uniqsites$ALT[k])
	#        }else if (length(uniqsites$ALT[k]) > 1 || length(uniqsites$REF[k]) > 1 || uniqsites$ALT[k] == "-"){
	#          dowrite = 0
	#        }else{
	#          motifseqsitemut = gsub("-","",uniqsites$ALT[k])
	#        }
	#        
	#        if (dowrite == 1) {
	#          if (totalmotifswritten < 2 &&
	#              uniqsites$start[k] != lastsite){
	#            totalmotifswritten = totalmotifswritten + 1
	#            lastsite = uniqsites$start[k]
	#            write(paste(genename,region.annotation,uniqsites$seqnames[k],uniqsites$start[k]-20,uniqsites$importance[k],
	#                        uniqsites$REF[k],
	#                        uniqsites$ALT[k],
	#                        paste(motifseqstart,motifseqsite,motifseqend,sep=""),
	#                        paste(motifseqstart,motifseqsitemut,motifseqend,sep=""),
	#                        paste(t(tmpdcast),collapse=" ")),paste("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/visualization/motif/regional/",sigtest,"_motif_gene",sep=""),append = T)
	#          }
	#          
	#         # write(paste(genename,region.annotation,uniqsites$seqnames[k],uniqsites$start[k],uniqsites$importance[k],
	#        #              uniqsites$REF[k],
	#        #              uniqsites$ALT[k],
	#        #              paste(t(tmpdcast),collapse=" ")),paste("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/visualization/entities/regional/",sigtest,"_entities_summary",sep=""),append = T)
	#          caseEntityCountingTable = rbind(caseEntityCountingTable,casesfulltable)
	#        }
	#      }
	#    }
	#    
	#    tmpdcast = data.frame(table(unique(caseEntityCountingTable)$Entity),stringsAsFactors = F)
	#    write(paste(genename,region.annotation,uniqsites$seqnames[k],uniqsites$start[k],
	#                              paste(t(tmpdcast),collapse=" ")),paste("/b06x-isilon/b06x-c/chromothripsis/results/icgc/stratton_breast/mutSig/Publication_Master/Association/visualization/entities/regional/",sigtest,"_entities_summary",sep=""),append = T)
	#                
	  }
	}
}

annotate_importance <- function(resultsimportancedf,gtfref){
		
		
		annotation_elements_priority = c("exon","UTR","gene")
		#search for gencode tags and use only basic
		if ("tag" %in% colnames(mcols(gtfref))){
			gtfref = gtfref[which(gtfref$tag == "basic" | (gtfref$type == "gene" & gtfref$level <= 2) )]
		}
		
		#create coordinates
	  sitestested = resultsimportancedf$test_site
	  sitestested = sitestested[which(sitestested != "NA")]
	  sitestesteddf = t(as.data.table(strsplit(sitestested,":")))
	  colnames(sitestesteddf) = c("chr","start")
	  chrarrayindex = as.numeric(sitestesteddf[1,1])
	  sitestesteddf = as.data.frame(sitestesteddf,stringsAsFactors=F,check.names=F)
	  sitestesteddf$test_site = paste(sitestesteddf$chr,":",sitestesteddf$start,sep="")
	  sitestesteddf$chr = paste("chr",sitestesteddf$chr,sep="")
	  sitestesteddf$chr = gsub("chr23","chrX",sitestesteddf$chr)
	  sitestesteddf$chr = gsub("chr24","chrY",sitestesteddf$chr)
	  sitestestedGR = makeGRangesFromDataFrame(sitestesteddf,
	  			seqnames.field="chr",
	  			start.field="start",
	  			end.field="start",keep.extra.columns=TRUE)
		sitestestedGR$annotation = NA
		sitestestedGR$annotationID = NA
		sitestestedGR$annotationtype = NA
		
		for (j in length(annotation_elements_priority):1){
			subsetgtf = gtfref[which(gtfref$type == annotation_elements_priority[j])]
			overlapanno = findOverlaps(sitestestedGR,subsetgtf)
			if (length(overlapanno) > 0){
				overlapannodf = as.data.table(overlapanno)
				overlapannodf = overlapannodf[!duplicated(overlapannodf$queryHits),]
				sitestestedGR[overlapannodf$queryHits]$annotation = subsetgtf[overlapannodf$subjectHits]$transcript_name
				sitestestedGR[overlapannodf$queryHits]$annotationID = subsetgtf[overlapannodf$subjectHits]$gene_id
				sitestestedGR[overlapannodf$queryHits]$annotationtype = annotation_elements_priority[j]
			}
		}
		#add missing ones
		overlapannodfNAidx = which(is.na(sitestestedGR$annotation))
		if (length(overlapannodfNAidx) > 0){
			missingAnnoGR = sitestestedGR[overlapannodfNAidx]
			missingAnnoGR$idx = overlapannodfNAidx
			#gene only
			subsetgtf = gtfref[which(gtfref$type == "gene")]
			hit_nearest = distanceToNearest(missingAnnoGR,subsetgtf)
			hit_nearestdf = as.data.table(hit_nearest,stringsAsFactors=F)
			hit_nearestdf = hit_nearestdf[order(hit_nearestdf$queryHits),]
			hit_nearestdf = data.table(hit_nearestdf,annotation=NA,annotationID=NA,annotationtype="gene")
			hit_nearestdf$annotation = subsetgtf[hit_nearestdf$subjectHits]$transcript_name
			hit_nearestdf$annotationID = subsetgtf[hit_nearestdf$subjectHits]$gene_id
			missingAnnoGR[hit_nearestdf$queryHits]$annotation = hit_nearestdf$annotation
			missingAnnoGR[hit_nearestdf$queryHits]$annotationID = hit_nearestdf$annotationID
			missingAnnoGR[hit_nearestdf$queryHits]$annotationtype = paste(hit_nearestdf$annotationtype,",dist=",hit_nearestdf$distance,sep="")
			sitestestedGR[missingAnnoGR$idx]$annotation = missingAnnoGR$annotation
			sitestestedGR[missingAnnoGR$idx]$annotationID = missingAnnoGR$annotationID
			sitestestedGR[missingAnnoGR$idx]$annotationtype = missingAnnoGR$annotationtype
			
		}
		
		sitestesteddf = as.data.table(sitestestedGR)
		
		resultsimportancedfannotated = merge(resultsimportancedf,sitestesteddf,by="test_site",all.x=T)
		resultsimportancetophits = resultsimportancedfannotated[order(resultsimportancedfannotated$importance,decreasing=T),]
		resultsimportancetophits = resultsimportancetophits[!duplicated(resultsimportancetophits$Region),]
		#put most important hit to annotation
		baslineidx = which(resultsimportancedfannotated$test_site == "NA")
		for ( y in 1:length(baslineidx)){
			resultsimportancedfannotated[baslineidx[y],]$annotation = resultsimportancetophits$annotation[which(resultsimportancetophits$Region == resultsimportancedfannotated[baslineidx[y],]$Region)]
			resultsimportancedfannotated[baslineidx[y],]$annotationID = resultsimportancetophits$annotationID[which(resultsimportancetophits$Region == resultsimportancedfannotated[baslineidx[y],]$Region)] 
			resultsimportancedfannotated[baslineidx[y],]$annotationtype = resultsimportancetophits$annotationtype[which(resultsimportancetophits$Region == resultsimportancedfannotated[baslineidx[y],]$Region)]
		}
		resultsimportancedfannotated$test_site[baslineidx] = "ALL"
		
		resultsimportancedfannotated = resultsimportancedfannotated[order(resultsimportancedfannotated$Region),]
		return(resultsimportancedfannotated)
}

importance_output_to_table <- function(resultsSKATanno){
	resultsimportancelist = list()
	for (i in 1:length(resultsSKATanno)){
		sitetable = data.frame(do.call(rbind,resultsSKATanno[[i]]),stringsAsFactors=F,check.names=F)
		resultsimportancelist[[length(resultsimportancelist) + 1]] = sitetable
	}
	resultsimportancedf = do.call(rbind,resultsimportancelist)
	colnames(resultsimportancedf) = c("Region","test_site","subregions","remaining_sites","marker_snvs","n.samples","Q","p.value","importance","entityvar")
	resultsimportancedf = resultsimportancedf[,c("Region","test_site","remaining_sites","marker_snvs","p.value","importance","entityvar")]
	resultsimportancedf$importance = as.numeric(resultsimportancedf$importance)
	resultsimportancedf$importance = -log10(resultsimportancedf$importance)
	
	return(resultsimportancedf)
}
