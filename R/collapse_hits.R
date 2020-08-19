

collapse_hits_by_window <- function(intable){
	require(data.table)

	#sort table by coordinates and remove NA
	intable = intable[which(grepl(":",intable$Region)),]
	#print(which(!grepl(":",intable$Region)))
	intable = intable[!is.na(intable$Tested_regions),]
	intable = intable[which(intable$Tested_regions != "NA"),]
	intable$p_value = as.numeric(intable$p_value)
	intable$Variants = as.numeric(intable$Variants)
	intable$Sites = as.numeric(intable$Sites)
	
	
	
	chrname = do.call(rbind,strsplit(intable$Region,":"))
	chrpos = do.call(rbind,strsplit(chrname[,2],"-"))

	intable$chr = chrname[,1]
	intable$start = as.numeric(chrpos[,1])
	intable$end = as.numeric(chrpos[,2])

	intable=intable[order(intable$chr,intable$start,intable$end),]
	intable$REMOVE = FALSE

	intablebycount=intable[order(intable$Variants,decreasing = T),]
	intablebycount$INDEX = match(intable$Region,intablebycount$Region)
	intablebycount = data.frame(intablebycount)
	#print(head(intablebycount))
	#setkey(intablebycount, INDEX)
	#intablebycount = as.data.table(intablebycount)


	print(paste("#Table size original:",dim(intable)[1]))
	#next biggest eats smaller tiplets
	#print(head(intablebycount))
	for (i in 1:(dim(intablebycount)[1]-1)){
	  #idx = which(intable$Region == intablebycount[i]$Region)
	  idx = intablebycount$INDEX[i]
	  
	  if (intable$REMOVE[idx] == FALSE && idx > 1 && idx < dim(intable)[1]){
	    #purn only when both +1,-1 window overlap with current idx
	    if ((intable$end[idx-1] > intable$start[idx] && intable$end[idx-1] < intable$end[idx] ) && 
	        (intable$start[idx+1] < intable$end[idx] && intable$start[idx+1] > intable$start[idx] )){

				#Check exceptions
				#1. significance(V8) drifts by log2 >= 3 AND n.var(V5) > 20, REASON: complex region cannot use sites only to determine segmentation
				if (((intable$p_value[idx] / intable$p_value[idx+1] > 100 && intable$Variants[idx+1] > 20) || 
						 (intable$p_value[idx] / intable$p_value[idx-1] > 100 && intable$Variants[idx-1] > 20)) &&
						 intable$Variants[idx] > 20){
					print(paste("Exceptions:",intable$chr[idx],intable$start[idx]))
					next
				}
				
	      intable$REMOVE[idx+1]=TRUE
	      intable$REMOVE[idx-1]=TRUE
	    }
	  }
	}
	
	#intable = intable[which(intable$REMOVE == FALSE),]
	print(paste("#Table size stage 1:",dim( intable[which(intable$REMOVE == FALSE),])[1]))
	#print(paste("#Table size stage 1:",dim(intable)[1]))
	intablebycount=intable[order(intable$Variants,decreasing = T),]
	intablebycount$INDEX = match(intable$Region,intablebycount$Region)
	#next biggest eats smaller doublets
	for (i in 1:(dim(intablebycount)[1]-1)){
	  #idx = which(intable$Region == intablebycount[i]$Region)
	  idx = intablebycount$INDEX[i]
	  #if (length(idx) > 0 && idx > 1 && idx < dim(intable)[1]){
	  if (intable$REMOVE[idx] == FALSE && idx > 1 && idx < dim(intable)[1]){
	    #purn only when both +1,-1 window overlap with current idx
	    #but do not leave empty spaces
	    if ((intable$end[idx-1] > intable$start[idx] && intable$end[idx-1] < intable$end[idx] ) && 
	       (idx - 2 < 1 || intable$REMOVE[idx-2] == FALSE || intable$end[idx-2] != intable$start[idx] )){
	       
	       
	       #Check exceptions
				#1. significance(V8) drifts by log2 >= 3 AND n.var(V5) > 20
				if ((intable$p_value[idx] / intable$p_value[idx-1] > 1000 &&
						 intable$Variants[idx-1] > 20) &&
						 intable$Variants[idx] > 20){
					print(paste("Exceptions:",intable$chr[idx],intable$start[idx]))
					next
				}
				
	      intable$REMOVE[idx-1] = TRUE
	    }else if((intable$start[idx+1] < intable$end[idx] && intable$start[idx+1] > intable$start[idx] ) &&
	             (idx + 2 > dim(intable)[1] || intable$REMOVE[idx+2] == FALSE  || intable$start[idx+2] != intable$end[idx])){
	             #(idx + 2 > dim(intable)[1] || intable$start[idx+2] != intable$end[idx+1])){

		      #Check exceptions
					#1. significance(V8) drifts by log2 >= 3 AND n.var(V5) > 20
					if ((intable$p_value[idx] / intable$p_value[idx+1] > 1000 && 
					     intable$Variants[idx+1] > 20)  &&
							 intable$Variants[idx] > 20){
						print(paste("Exceptions:",intable$chr[idx],intable$start[idx]))
						next
					}
		
		      intable$REMOVE[idx+1] = TRUE
	    }
	  }
	}


	intable = intable[which(intable$REMOVE == FALSE),]
	print(paste("#Table size stage 2:",dim(intable)[1]))
	  

	#write.table(intable,outfile,sep=" ",quote = F,row.names = F)
	return(intable)
}