# wigfilename: the name of the wiggle file
# annotationfilename: the names of annotation file hg18
# width: the window size upstream
## smooth: Logic value, if smooth=T perform smoothing for the TF binding profile;  

wig2weight<-function(wigfilename, annotationfilename, width=10000, smooth=T, outfilename)
{
	chr.len = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
					146274826,140273252,135374737,134452384,132349534,114142980,106368585,
					100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,
					154913754,57772954,16571)
	chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
	chr.len = chr.len[1:(length(chr.len)-2)]
	chr.nam = chr.nam[1:(length(chr.nam)-2)]

	# read in gene info
	mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
	tag = substr(mygene[,1], 1, 3)
	mygene = mygene[tag=="NM_", 1:5]
	#care about mRNA, from tts to tts..
	colnames(mygene) = c("name", "chr", "str", "sta", "end")
	tmp = mygene$sta
	tmp1 = tmp-width
	tmp2 = tmp+width
	tag = mygene[,3]=="-"
	tmp[tag==1] = mygene[tag==1,5]
	tmp1[tag==1] = tmp[tag==1]+width 
	tmp2[tag==1] = tmp[tag==1]-width 
	mygene$sta = tmp1
	mygene$end = tmp2
	#at the end, consider only -10000 to 10000 around TSS for all NM,
	#make sure st and ed won't go beyond the boundaries
	for(k in 1:length(chr.nam))
	{
	   tag = mygene$chr==chr.nam[k]
	   tmp = mygene$sta[tag==1]
	   tmp[tmp<1] = 1
	   tmp[tmp>chr.len[k]] = chr.len[k]	
	   mygene$sta[tag==1] = tmp	
	   tmp = mygene$end[tag==1]
	   tmp[tmp<1] = 1
	   tmp[tmp>chr.len[k]] = chr.len[k]	
	   mygene$end[tag==1] = tmp	
	}
	
	# read in wiggle file
	conIn = file(wigfilename, "r")
	data = readLines(conIn, -1)
	close(conIn)
	
	###
	pos = grep("fixedStep", data)
	info = data[pos]
	pos = c(pos, length(data)+1)
	info = unlist(strsplit(info, " "))
	cnum = length(info)/4
	mychr = info[(1:cnum)*4-2]
	mysta = info[(1:cnum)*4-1]
	mychr = gsub("chrom=", "", mychr)
	mysta = as.numeric(gsub("start=", "", mysta))
	pos.sta = pos[1:(length(pos)-1)]+1
	pos.end = pos[2:length(pos)]-1
	sig.tk = data.frame(mychr, mysta, pos.sta, pos.end)

	## create weight file
	myw = rep(0, width*2+1)
	for(k in 1:length(chr.nam))
	{
	   	cat("\r", chr.nam[k])
	   	# signal vector
	   	read.cov = rep(0, chr.len[k])
	   	tag = sig.tk[,1]==chr.nam[k]
	   	if(sum(tag)==0) next
	   	mysig = sig.tk[sig.tk[,1]==chr.nam[k], ]
 	  	for(i in 1:nrow(mysig))
		{
			tmp1 = mysig[i,3]   # start line
			tmp2 = mysig[i,4]   # end line
			tmp.sta = mysig[i,2]  # start position in chr
			read.cov[(tmp.sta+1):(tmp.sta+tmp2-tmp1+1)] = as.numeric(data[tmp1:tmp2])
		}
   	# refseq
   		curgene = mygene[mygene[,2]==chr.nam[k],]
   		for(i in 1:nrow(curgene))
   		{
   	  	 	myw = myw+read.cov[curgene[i,4]:curgene[i,5]]   	      	   
   		}
 	}
	myw = myw/nrow(mygene)
	tmp = myw
	if(smooth==T)
	{
		for(i in 1:length(myw))
		{
		    myw[i] = mean(tmp[max(i-250,1):min(i+250, length(myw))])	
		}
	}
	myw = myw/sum(myw)
	write.table(myw, outfilename, sep="\t", row.names=F, col.names=F, quote=F)
}

#####

peak2weight<-function(peakfilename, annotationfilename, width=10000, smooth=T, outfilename)
{
	chr.len = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
					146274826,140273252,135374737,134452384,132349534,114142980,106368585,
					100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,
					154913754,57772954,16571)
	chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
	chr.len = chr.len[1:(length(chr.len)-2)]
	chr.nam = chr.nam[1:(length(chr.nam)-2)]

	# read in gene info
	mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
	tag = substr(mygene[,1], 1, 3)
	mygene = mygene[tag=="NM_", 1:5]
	colnames(mygene) = c("name", "chr", "str", "sta", "end")
	tmp = mygene$sta
	tmp1 = tmp-width
	tmp2 = tmp+width
	tag = mygene[,3]=="-"
	tmp[tag==1] = mygene[tag==1,5]
	tmp1[tag==1] = tmp[tag==1]+width 
	tmp2[tag==1] = tmp[tag==1]-width 
	mygene$sta = tmp1
	mygene$end = tmp2
	for(k in 1:length(chr.nam))
	{
	   tag = mygene$chr==chr.nam[k]
	   tmp = mygene$sta[tag==1]
	   tmp[tmp<1] = 1
	   tmp[tmp>chr.len[k]] = chr.len[k]	
	   mygene$sta[tag==1] = tmp	
	   tmp = mygene$end[tag==1]
	   tmp[tmp<1] = 1
	   tmp[tmp>chr.len[k]] = chr.len[k]	
	   mygene$end[tag==1] = tmp	
	}
	
	# read in peak file
	data = read.table(peakfilename, sep="\t", header=T)
	data$chr = paste("chr", data$chr, sep="")
	
	## create weight file
	myw = rep(0, width*2+1)
	for(k in 1:length(chr.nam))
	{
	   	cat("\r", chr.nam[k])
	   	# signal vector
	   	read.cov = rep(0, chr.len[k])
	   	tag = data$chr==chr.nam[k]
	   	if(sum(tag)==0) next
	   	mysig = data[tag==1,]
 	  	for(i in 1:nrow(mysig))
		{
			tmp1 = as.numeric(mysig$start[i])+1   # start line
			tmp2 = as.numeric(mysig$end[i])   # end line
			read.cov[tmp1:tmp2] = read.cov[tmp1:tmp2] + as.numeric(mysig$height[i])
		}
   	# refseq
   		curgene = mygene[mygene[,2]==chr.nam[k],]
   		for(i in 1:nrow(curgene))
   		{
   	  	 	myw = myw+read.cov[curgene[i,4]:curgene[i,5]]   	      	   
   		}
 	}
	myw = myw/nrow(mygene)
	tmp = myw
	if(smooth==T)
	{
		for(i in 1:length(myw))
		{
		    myw[i] = mean(tmp[max(i-250,1):min(i+250, length(myw))])	
		}
	}
	myw = myw/sum(myw)
	write.table(myw, outfilename, sep="\t", row.names=F, col.names=F, quote=F)
}
#####

calscoreonwig<-function(wigfilename, annotationfilename, weightfilename, outfilename)
{
	chr.len = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
					146274826,140273252,135374737,134452384,132349534,114142980,106368585,
					100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,
					154913754,57772954,16571)
	chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
	chr.len = chr.len[1:(length(chr.len)-2)]
	chr.nam = chr.nam[1:(length(chr.nam)-2)]

	# read in weight
	conIn = file(weightfilename, "r")
	myw = readLines(conIn)
	close(conIn)
	myw = as.numeric(myw)
	width = (length(myw)-1)/2
	
	# read in gene info
	mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
	tag = substr(mygene[,1], 1, 3)
	mygene = mygene[tag=="NM_", 1:5]
	colnames(mygene) = c("name", "chr", "str", "sta", "end")
	
	# read in wiggle file
	conIn = file(wigfilename, "r")
	data = readLines(conIn, -1)
	close(conIn)
	
	###
	pos = grep("fixedStep", data)
	info = data[pos]
	pos = c(pos, length(data)+1)
	info = unlist(strsplit(info, " "))
	cnum = length(info)/4
	mychr = info[(1:cnum)*4-2]
	mysta = info[(1:cnum)*4-1]
	mychr = gsub("chrom=", "", mychr)
	mysta = as.numeric(gsub("start=", "", mysta))
	pos.sta = pos[1:(length(pos)-1)]+1
	pos.end = pos[2:length(pos)]-1
	sig.tk = data.frame(mychr, mysta, pos.sta, pos.end)

	## calculate score
	mysco = rep(0, nrow(mygene))	#TSS
	gname = rep("",nrow(mygene))
	count = 0
	for(k in 1:length(chr.nam))
	{
	   	cat("\r", chr.nam[k])
	   	# signal vector
	   	read.cov = rep(0, chr.len[k])
	   	tag = sig.tk[,1]==chr.nam[k]
	   	if(sum(tag)==0) next
	   	mysig = sig.tk[sig.tk[,1]==chr.nam[k], ]
 	  	for(i in 1:nrow(mysig))
		{
			tmp1 = mysig[i,3]   # start line
			tmp2 = mysig[i,4]   # end line
			tmp.sta = mysig[i,2]  # start position in chr
			read.cov[(tmp.sta+1):(tmp.sta+tmp2-tmp1+1)] = as.numeric(data[tmp1:tmp2])
		}
   		curgene = mygene[mygene[,2]==chr.nam[k],]
	   	for(i in 1:nrow(curgene))
	   	{
	   	   count = count + 1
	   	   if(curgene$str[i]=="+")
	   	   {
	   	   		tmp= read.cov[max(1,curgene$sta[i]-width):min(curgene$sta[i]+width, chr.len[k])]
	   	   		b1 = 1-(curgene$sta[i]-width)
	   	   		b2 = curgene$sta[i]+width-chr.len[k]
	   	   		if(b1>0)
	   	   		{
	   	   		   tmp=c(rep(0, b1), tmp)	
	   	   		}
	   	   		if(b2>0)
	   	   		{
	   	   		   tmp=c(tmp, rep(0, b2))	
	   	   		}
	   	   		mysco[count] = sum(tmp*myw)
	   	   }
	   	   if(curgene$str[i]=="-")
	   	   {
	   	   		tmp= read.cov[min(chr.len[k],curgene$end[i]+width):max(curgene$end[i]-width, 1)]
	   	   		b1 = curgene$end[i]+width-chr.len[k]
	   	   		b2 = 1-(curgene$end[i]-width)
	   	   		if(b1>0)
	   	   		{
	   	   		   tmp=c(rep(0, b1), tmp)	
	   	   		}
	   	   		if(b2>0)
	   	   		{
	   	   		   tmp=c(tmp, rep(0, b2))	
	   	   		}
	   	   		mysco[count] = sum(tmp*myw)
	   	   }
	   	   gname[count] = curgene[i,1]   	      	   
	   	}
	}
	zscore = (mysco-mean(mysco))/sd(mysco)	#Z-score based onTSS
	pvalue = pnorm(-zscore)
    res = cbind(gname, mysco, zscore, pvalue)
	colnames(res) = c("name", "raw.score", "zscore", "p.value")
	res = res[order(res[,1], decreasing=T), ]
	write.table(res, outfilename, sep="\t", row.names=F, quote=F)
}


calscoreonpeak<-function(peakfilename, annotationfilename, weightfilename, outfilename)
{
	chr.len = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
					146274826,140273252,135374737,134452384,132349534,114142980,106368585,
					100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,
					154913754,57772954,16571)
	chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
	chr.len = chr.len[1:(length(chr.len)-2)]
	chr.nam = chr.nam[1:(length(chr.nam)-2)]

	# read in weight
	conIn = file(weightfilename, "r")
	myw = readLines(conIn)
	close(conIn)
	myw = as.numeric(myw)
	width = (length(myw)-1)/2
	
	# read in gene info
	mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
	tag = substr(mygene[,1], 1, 3)
	mygene = mygene[tag=="NM_", 1:5]
	colnames(mygene) = c("name", "chr", "str", "sta", "end")
	
	# read in peak file
	data = read.table(peakfilename, sep="\t", header=T)
	data$chr = paste("chr", data$chr, sep="")
	
	## calculate score
	mysco = rep(0, nrow(mygene))	#TSS
	gname = rep("",nrow(mygene))
	count = 0
	for(k in 1:length(chr.nam))
	{
	   	cat("\r", chr.nam[k])
	   	# signal vector
	   	read.cov = rep(0, chr.len[k])
	   	tag = data$chr==chr.nam[k]
	   	if(sum(tag)==0) next
	   	mysig = data[tag==1,]
 	  	for(i in 1:nrow(mysig))
		{
			tmp1 = as.numeric(mysig$start[i])+1   # start line
			tmp2 = as.numeric(mysig$end[i])   # end line
			read.cov[tmp1:tmp2] = read.cov[tmp1:tmp2] + as.numeric(data$height[i])
		}
   	# refseq
   		curgene = mygene[mygene[,2]==chr.nam[k],]
	   	for(i in 1:nrow(curgene))
	   	{
	   	   count = count + 1
	   	   if(curgene$str[i]=="+")
	   	   {
	   	   		tmp= read.cov[max(1,curgene$sta[i]-width):min(curgene$sta[i]+width, chr.len[k])]
	   	   		b1 = 1-(curgene$sta[i]-width)
	   	   		b2 = curgene$sta[i]+width-chr.len[k]
	   	   		if(b1>0)
	   	   		{
	   	   		   tmp=c(rep(0, b1), tmp)	
	   	   		}
	   	   		if(b2>0)
	   	   		{
	   	   		   tmp=c(tmp, rep(0, b2))	
	   	   		}
	   	   		mysco[count] = sum(tmp*myw)
	   	   }
	   	   if(curgene$str[i]=="-")
	   	   {
	   	   		tmp= read.cov[min(chr.len[k],curgene$end[i]+width):max(curgene$end[i]-width, 1)]
	   	   		b1 = curgene$end[i]+width-chr.len[k]
	   	   		b2 = 1-(curgene$end[i]-width)
	   	   		if(b1>0)
	   	   		{
	   	   		   tmp=c(rep(0, b1), tmp)	
	   	   		}
	   	   		if(b2>0)
	   	   		{
	   	   		   tmp=c(tmp, rep(0, b2))	
	   	   		}
	   	   		mysco[count] = sum(tmp*myw)
	   	   }
	   	   gname[count] = curgene[i,1]   	      	   
	   	}
	}
	zscore = (mysco-mean(mysco))/sd(mysco)	#Z-score based onTSS
	pvalue = pnorm(-zscore)
    res = cbind(gname, mysco, zscore, pvalue)
	colnames(res) = c("name", "raw.score", "zscore", "p.value")
	res = res[order(res[,1], decreasing=T), ]
	write.table(res, outfilename, sep="\t", row.names=F, quote=F)
}


