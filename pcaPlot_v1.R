#=========Thoughts==================
# As for the samples, we need to elborate directly from the graphs that where each sample belongs to.
# Plan for now would be the followings:
# separate GWAS1 and 2 by thickness (cex)
# separate Cohorts by symbols (col)
# separate platforms by colors (pch)
# Ideas are the same as the previous version (150510)
# Just need to distinguish the sample dots as above
# So some function can still be useable from the previous version (150510)
######UPDATED: MAY 13th######
# Separation was done differently since cex is not visible and no pch option
# So GWAS1/2 was shrinked to platform (610 for GWAS1, and remaining for GWAS2)
# use letter to represent cohorts
# different colors for platforms

#=========Functions=================
boxArea <- function(dataPoints) {
  # PRE: given 3 vectors of data points like 3 PCs in data.frame format
  # POST: returns the xrange, yrange and zrange of the 3D box as a vector 
  #     in the following format: c(xmin, xmax, ymin, ymax, zmin, zmax)
  xmin <- min(dataPoints[,1])
  xmax <- max(dataPoints[,1])
  ymin <- min(dataPoints[,2])
  ymax <- max(dataPoints[,2])
  zmin <- min(dataPoints[,3])
  zmax <- max(dataPoints[,3])
  return(c(xmin, xmax, ymin, ymax, zmin, zmax))
}
isInRange <- function(range, point) (point >= range[1] & point <= range[2])
isInBox <- function(hap, sample) {
  # PRE: get the 'box' from hapmap data for specific ethnic group in hap and 
  #       the sample individuals to be tested for belonging in the
  #       cluster under sample
  # POST: returns boolean vector on whether individual is within ethnic group
  inBox <- logical(length=dim(sample)[1])
  box <- boxArea(hap[,7:9])
  xrange <- c(box[1], box[2])
  yrange <- c(box[3], box[4])
  zrange <- c(box[5], box[6])
  for (i in 1:dim(sample)[1]) {
    if (isInRange(xrange, sample[i,7]) & 
        isInRange(yrange, sample[i,8]) &
        isInRange(zrange, sample[i,9]))
    {
      inBox[i] <- TRUE
    } else {
      inBox[i] <- FALSE
    }
  }
  return(inBox)
}
sepeth= function(hapmap_vec,relationship,populations){
  # read in  the hapmap individuals
  # Return a list with each population in one element
  # ordering is: 
  sep_list = list(); 
  for (i in 1:length(populations)){
    pop = subset(relationship, relationship$population %in% populations[i])
    eigen = subset(hapmap_vec, hapmap_vec$V2 %in% pop$IID)
    sep_list = c(sep_list, list(eigen))
  }
  return(sep_list)
}
inWhichBox = function(sample, hap_list, populations){
  # sample -> inds to be tested
  # hap_list -> list of different ethnicity groups calculated from hapmap samples
  # This function calls the above isInBox function
  # return the list of inds that falls in other ethnicity boxes
  rmv_list=NULL;
  name = populations[!(populations %in% c('CEU','TSI'))]
  for (i in 1:length(hap_list)){
    rmv = sample[isInBox(hap_list[[i]],sample),]
    rmv_list = rbind(rmv_list, data.frame(FID=rmv[,1],IID =rmv[,2],ETH =rep(name[i],dim(rmv)[1]),stringsAsFactors=F))
  }
  return(rmv_list)
}
getGroup = function(ind_info){
  # helper function prepare group stratification for plot_sample
  # Cohorts are represented by their first letter
  # Up tp 6 Platforms are saved as color, in particular:
  # 370K -> 'magenta', 610-> 'black', 660K -> 'limegreen',
  # 660W -> 'blue1', GWAS3 -> 'peru', Omni5 -> 'darkmagenta'
  txt=rep(NA, dim(ind_info)[1])
  col = rep('red', dim(ind_info)[1])
  if ('cohort' %in% tolower(names(ind_info))){
    txt = substring(levels(as.factor(ind_info$Cohort)),1,1) #take the first letter as the representation
  } 
  if ('platform' %in% tolower(names(ind_info))){
    col = c('magenta','black','limegreen','blue1','peru','darkmagenta')[as.numeric(as.factor(ind_info$Platform))]
  }
  
  return(data.frame(FID=ind_info$FID,IID=ind_info$IID,TXT=txt,COL=col))
}
getOutlier = function(hap, sample,pcs=3) {
  # hap map are the people with selected ethnicity  
  # samples are the people not in hapmap
  # pcs -> numbers of PCs wanna use, likely 3, could be 2
  # output the column of those people that are outside the range
  
  hap_sd = sapply(7:(7+pcs-1), function(i) sd(as.numeric(hap[,i]),na.rm=T))
  hap_mean = sapply(7:(7+pcs-1), function(i) mean(as.numeric(hap[,i]),na.rm=T))
  left_limit = hap_mean - 6*hap_sd
  right_limit= hap_mean + 6*hap_sd
  outlier = NULL
  for(i in 7:(7+pcs-1)) {
    out = subset(sample,(as.numeric(sample[,i])<left_limit[i-6] | as.numeric(sample[,i])>right_limit[i-6])) 
    outlier = rbind(outlier,out)
  }
  return(unique(outlier))
}
plot_sample3d = function(sample_vec,ind_info,group){
  # This is the undernearth function that serves the main function plot_eth
  # Plot samples by groups
  # The groups are categorized by the combination of cohort and platform
  if (dim(sample_vec)[1] != dim(ind_info)[1]){
    ind_info = subset(ind_info, IID %in% sample_vec[,2])
  }

  group$GRU = as.factor(paste(group$TXT,'-',group$COL,sep=''))  
  for (i in 1:length(levels(group$GRU))){
    groupsample = sample_vec[group$GRU == levels(group$GRU)[i],] 
    m = dim(groupsample)[1]
    grouptxt = as.character(strsplit(levels(group$GRU)[i],'-')[[1]][1])
    groupcol = as.character(strsplit(levels(group$GRU)[i],'-')[[1]][2])
    if (prod(grouptxt == 'NA')){
      points3d(groupsample[,7],groupsample[,8],groupsample[,9],col=rep(groupcol,m),size=6)
    } else{
    text3d(groupsample[,7],groupsample[,8],groupsample[,9],text=rep(grouptxt,m),
           col=rep(groupcol,m))
    }
    legend3d('topright', legend = levels(as.factor(ind_info$Platform)),
                         col=c('magenta','black','limegreen','blue1','peru','red')[1:length(levels(as.factor(ind_info$Platform)))],
                         pch=16)
    
}}
plot_sample2d = function(sample_vec,ind_info,group,v1,v2){
  # This is the undernearth function that serves the main function plot_eth
  # if the sample_vec and ind_info don't agree, pick the subset of ind_info
  # Plot samples by groups
  # The groups are categorized by the combination of cohort and platform
  if (dim(sample_vec)[1] != dim(ind_info)[1]){
    ind_info = subset(ind_info, IID %in% sample_vec[,2])
  }
  
  group$GRU = as.factor(paste(group$TXT,'-',group$COL,sep=''))  
  for (i in 1:length(levels(group$GRU))){
    groupsample = sample_vec[group$GRU == levels(group$GRU)[i],] 
    m = dim(groupsample)[1]
    grouptxt = strsplit(levels(group$GRU)[i],'-')[[1]][1]
    grouppch = as.numeric(sum((c('F','J','T','U') %in% grouptxt)*c(0,1,2,3)))
    groupcol = as.character(strsplit(levels(group$GRU)[i],'-')[[1]][2])
    points(groupsample[,v1],groupsample[,v2], col=rep(groupcol,m),pch=rep(grouppch,m))
  }
}

plot_eth3d= function(pc_vec, ind_info,group=getGroup(ind_info),relationship, populations  
                     ,colours){
  # This is the main function, the first part built the axis and the hapmap clusters
  # and call plot_sample to plot samples onto it with different symbols representing different cohort/platforms
  # Legend created as well
  # This plot function is compatiable without hapmap individuals (Although the legend stays there)
  library('rgl')
  hapmap_vec =pc_vec[pc_vec[,2] %in% relationship$IID,]
  sample_vec = pc_vec[!(pc_vec[,2] %in% relationship$IID),]
  with(pc_vec, plot3d(V7,V8,V9, type="n",xlab='PC1',ylab='PC2',zlab='PC3'))
  for (i in 1:length(populations)){
    pop = subset(relationship, relationship$population %in% populations[i])
    eigen = subset(hapmap_vec, hapmap_vec$V2 %in% pop$IID)
    with(eigen, points3d(V7,V8,V9, col=colours[i],size=6))
  }
  plot_sample3d(sample_vec, ind_info,group)
}

plot_eth2d = function(pc_vec, ind_info,group=getGroup(ind_info),relationship, coordinate,
                      populations,colours){
  # This is the main function, the first part built the axis and the hapmap clusters
  # and call plot_sample to plot samples onto it with different symbols representing different cohort/platforms
  # input coordinates 'X','Y','Z' to represent PC1-3
  # Legend created as well
  # This plot function is compatiable without hapmap individuals (Although the legend stays there)
  c1 = coordinate[1];c2=coordinate[2];
  v1 = sum((c('X','Y','Z') %in% toupper(c1))*c(1,2,3))+6 # The sum thing returned 1,2,3 for being x,y,z;
  v2 = sum((c('X','Y','Z') %in% toupper(c2))*c(1,2,3))+6 # Corresponding to 7th, 8th, 9th Columns
  hapmap_vec =pc_vec[pc_vec[,2] %in% relationship$IID,]
  sample_vec = pc_vec[!(pc_vec[,2] %in% relationship$IID),]
  plot(pc_vec[,v1],pc_vec[,v2], type="n",xlab=paste0('PC',sum((c('X','Y','Z') %in% toupper(c1)*c(1,2,3))))
       ,ylab=paste0('PC',sum((c('X','Y','Z') %in% toupper(c2))*c(1,2,3))))
  for (i in 1:length(populations)){
    pop = subset(relationship, relationship$population %in% populations[i])
    eigen = subset(hapmap_vec, hapmap_vec$V2 %in% pop$IID)
    points(eigen[,v1],eigen[,v2], col=colours[i],pch=16)
  }
  legend('topright',legend=populations ,col=colours ,pch=rep(16,length(populations), cex=0.7))
  
  #   legend('topleft',legend = c(levels(as.factor(ind_info$Platform)),c('FR','JHU','TOR','UNC')),
  #          col= c('magenta','black','limegreen','blue1','peru','red',rep('black',4)),
  #          pch = c(rep(16,length(levels(as.factor(ind_info$Platform)))),0:3),cex=0.7)}
  plot_sample2d(sample_vec, ind_info,group,v1=v1,v2=v2)
}

plotLegend=function(populations,color){
  n=length(color)
  par(mfrow=c(1,1))
  plot(c(0,1),c(0,12),type='n',xlab='',ylab='')
  x1=rep(0,n)
  x2=rep(0.2,n)
  y1=0:(n-1)
  y2=1:n
  rect(x1,y1,x2,y2,col=color)
  #x.text=rep(0.4,12)
  #y.text=1:12
  #text(x.text,y.text,express
  ref = c('ASW-African','CEU-Caucasian','CHB-Chinese-Beijing','CHD-Chinese-Denver','GIH-Indian'
          ,'JPT-Japanese','LWK-Luhya','MEX-Mexican','MKK-Maasia','TSI-Italian','YRI-Yoruba')
  for (i in 1:n){
    exp =ref[grep(paste0('^',populations[i]),ref)] 
    text(0.4,0.5+i-1,exp)
  }
}
