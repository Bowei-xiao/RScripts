# This function created the forest plot for particular SNPs
# THe graph contains 4 cohorts: FR JHU TOR UNC
# It allowed show Gene Name on the title
# Choice of weight and CI
# Method was set to random effect but could be changed
# Input file should be in the form of chrx.blockxx_{cohort}
# Data should contain 4 columns: snp, beta, se, p-value without header
forest_plot = function(chr,blk,snp,weight=T,ci=F,method='REML',gene='Unknown'){
  library('data.table')
  library('metafor')  
  file_list = list.files(pattern = paste0('^chr',chr,'.block',blk,'_+'))
  frst=NULL;
  for (i in (1:length(file_list))){
    file= fread(file_list[i],header=F, sep=' ',data.table = F)
    inx = which(file$V1 == snp)
    frst = rbind(frst, data.frame(file=file_list[i], mean = file[inx,2]
                                  ,se = file[inx,3],p=file[inx,4],stringsAsFactors = F))
  }
  x =rma(yi=as.numeric(mean), sei=as.numeric(se), data=frst,method=method,slab = c('FR(1215)','JHU(1514)','TOR(1569)','UNC(1547)'))
  png(paste0('forPlot_',snp,'On_',gene,'.png'), width=1920, height=1440, pointsize=24)
  forest(x, showweights=weight ,addcred=ci,
         transf=exp,xlab=paste0('Meta Analysis of Odds Ratio for SNP ',snp,' on Gene ',gene),refline=1
         ,order=order(frst$file))
  mtext('Cohorts',line=-7,col='black',adj=0.03)
  mtext('Weights',line=-7,col='black',adj=0.83)
  mtext('Odds Ratio [95% CI]',line=-7,col='black',adj=1)
  dev.off()
  print('Plot ouputed')
  return(x)
}
if(FALSE){
  forest_plot(1,21,'rs7419153',gene='SLC26A9') #I=0
  forest_plot(13,3,'rs61948108',gene='ATP12A')  #I=0
  forest_plot(23,12,'rs3788766',gene='SLC6A14') #I=0
  forest_plot(2,5,'rs72796743',gene='ABCG8')   #I=68.79%
  forest_plot(7,15,'rs9969188', gene='PRSS1')  #I=0
  
}