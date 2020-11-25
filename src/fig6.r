#' Fig. 6: Alignment with later stage human data: eye development as an example.
#' Retinal cell types at late stage: Hu, Y. et al. Dissecting the transcriptome landscape of the human fetal neural retina and retinal pigment epithelium by single-cell RNA-seq analysis. PLoS biology 17, e3000365, doi:10.1371/journal.pbio.3000365 (2019).

# Load both datasets
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
knitr::opts_chunk$set(warning = FALSE)
load('our_mx.rdata')
load('../../../cross_embryo/code/ge2an.rdata')
ficlu<-read.table('~/BaoLab/xuy/human/cross_embryo/result/reclustering3/TJ_clu/cell_annotate.csv',sep=',',as.is=T,header=T)
# clean Hu's data
get_id<-function(x){ sapply(x, function(y){names(ge2an)[ge2an==y][1]})}
tang_cell<-read.table('../data/s5.txt',sep='\t',as.is=T,header=T) #
tang_mx<- read.table('../data/GSE107618_Merge.TPM.csv',sep=',',as.is=T,header=T) #
sum( tang_cell[,8] %in% colnames(tang_mx) ) # all cells with ID have UMI data
rownames(tang_cell)<- tang_cell[,8]
tang_deg<- read.table('../data/s4_deg.txt',sep='\t',as.is=T,header=T)
length(late_cell<- rownames(tang_cell)[tang_cell[,6] %in% c('5W','6W') & (!tang_cell[,7] %in% c('Blood','Fibroblast')) ]) #
sum(ge2an %in% tang_mx[,1] ) #
late_mx<- data.matrix(tang_mx[ tang_mx[,1] %in% ge2an, late_cell ])
rownames(late_mx)<- tang_mx[tang_mx[,1] %in% ge2an,1]
late_mx<- late_mx/100 # TPM to per 10000
rownames(late_mx)<- get_id( rownames(late_mx) )


# Merge dataset and calculate umap
mg_mx<- cbind( our_mx[ rownames(late_mx) , unlist(neu1_ref_clu[c(6,10,13)]) ], late_mx)

# compile an annotation for merged data
tmp<-unique(read.table('~/baolab/xuy/ciona/cluster2/product/code/type_all_col3.txt',sep='\t',as.is=T,comment.char='')[,2])
col6<-c('lightgrey','red',tmp[3],'purple',tmp[2],tmp[5],'blue','green4','orchid','turquoise','sienna',tmp[9],'yellow2','hotpink','navy','steelblue','skyblue','pink','black',tmp[4],rainbow(7))
cell_id<-rep('',dim(mg_mx)[2])
names(cell_id)<- colnames(mg_mx)
cell_id[neu1_ref_clu[[6]]]<- 'o-rpe'
cell_id[names(cell_id)%in% ficlu[ficlu[,4]=='retinal progenitor cell',1] ]<- 'o-rpc'
cell_id[names(cell_id)%in% ficlu[ficlu[,4]=='optic vesicle',1]]<- 'o-ov'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('RGC_1','RGC_2')] ]<- 't-rgc'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('RPC_1','RPC_2','RPC_3','RPC_4','RPC_5','RPC_6')] ]<- 't-rpc'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('RPE_1','RPE_3')] ]<- 't-rpe'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('Microglia')] ]<- 't-mg'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('Horizontal_2','Horizontal_1')] ]<- 't-hor'
cell_id[ late_cell[tang_cell[late_cell,7] %in% c('Amacrine_1')] ]<- 't-ama'
# Color theme for cell types and stages
type_col<-col6[c(6,2,3,27,25,24,4,11,12)]
names(type_col)<-names(table(cell_id))
time_col<-col6[c(23,25,25,25,25,22,23,25,25,25,22,23,25,25,25,23,25,25,25,27,27)]
names(time_col)<-c('h0','h9b','h9a','h5','h6','ht7','l0','l9','lv5','lv6','tv7','t0','t9','t5','t6','v0','v9','v6','vl5','RET','RPE')

# Union of signature genes for pseudo-time analysis
our_sig<- unique(unlist(type_diag[c(4,6,7)]))
late_sig<- get_id(unique(unlist(lapply( c('RPE','RPC','RGC'), function(x,nn=30){
  deg<-tang_deg[ tang_deg[,8]==x, ]
  deg<- deg[ order(-deg[,1]), ]
  n<- ifelse( dim(deg)[1]<nn, dim(deg)[1], nn)
  return(deg[1:n,9])
})))) # control number of DEGs
length(mg_sig<- setdiff( unique(c(our_sig,late_sig)), rm_gene)) 

#' Pseudo-time analysis
b2_cell<-  colnames(mg_mx)[ cell_id[colnames(mg_mx)] %in% c('o-rpc','o-ov','o-rpe','t-rpc','t-rpe','t-rgc') ]
b2_anno<- data.frame(cell_id[b2_cell])
library(monocle)
b2_mono<- newCellDataSet( mg_mx[,b2_cell],expressionFamily=negbinomial.size(),phenoData = new("AnnotatedDataFrame", data =b2_anno ))
b2_mono <- estimateSizeFactors(b2_mono)
b2_mono <- estimateDispersions(b2_mono,remove_outliers =F)
b2_mono <- setOrderingFilter(b2_mono,mg_sig)
b2_mono<-reduceDimension(b2_mono, max_components = 2, reduction_method = 'DDRTree', ncenter=51)
b2_mono <- orderCells(b2_mono)
b2_rand<- sample.int(length(b2_cell))
#cstat<- b2_mono@phenoData@data[,4]
#pdf('../result/mono.pdf',width=7,height=7)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3,pch=20,mfrow=c(2,2))
plot( b2_mono@reducedDimS[1,b2_rand],b2_mono@reducedDimS[2,b2_rand],pch=16,xlab='dim-1',ylab='dim-2',frame=F,col=type_col[cell_id[b2_cell[b2_rand]]], cex=.5,main='cell type')
ind<- cell_id[colnames(b2_mono@reducedDimS)] %in% 'o-rpc'
points( b2_mono@reducedDimS[1, ind],b2_mono@reducedDimS[2,ind], col=type_col['o-rpc'], cex=.5, pch=16)
legend('bottomright', bty='n', legend=c('o-rpc','o-ov','o-rpe','t-rpc','t-rpe','t-rgc') , col=type_col[c('o-rpc','o-ov','o-rpe','t-rpc','t-rpe','t-rgc') ], cex=.8,xpd=NA, pch=16)
plot( b2_mono@reducedDimS[1,b2_rand],b2_mono@reducedDimS[2,b2_rand],pch=16,xlab='dim-1',ylab='dim-2',frame=F,col=time_col[sapply(b2_cell[b2_rand],function(x){strsplit(x,split='_')[[1]][1]})], cex=.5,main='stage')
legend( 'bottomright',bty='n', legend=c('CS12','CS13','CS15','5-6W') , col=time_col[c(6,1,2,20) ], cex=.8,xpd=NA,pch=16)
#dev.off()

#' Plot marker expression on pseudo-time
source('../../../AllCode/function.r')
plot_genes_on_tsne<-function( tsne, mx, genes, file_name,is_order=T, num=length(genes),path='',plot_cex=.5,is_ret=F, data_max=0, zerof=T, xx=c(0,0), yy=c(0,0), title_lab=''){
  if(is_order){
     if(length(genes)==1) exp_num<-sum(mx[genes,]>0) else exp_num<- apply(mx[ genes,  ]>0,1,sum) 
     plot_gene<-genes[order(exp_num)] 
  }else plot_gene<-genes
  plot_gene<-plot_gene[1:num]
  if(xx[1]==0&xx[2]==0){
      xx<-range(tsne[,1])
      yy<-range(tsne[,2])
  }
  if(title_lab[1]=='') title_lab<-ge2an[plot_gene]
  if(data_max[1]==0) data_max<-rep(0,length(plot_gene))
#  pdf(paste(path,file_name,'.pdf',sep=''),height=12,width=9)
  par(las=1,mar=c(5,1,1,4),lwd=2, mfrow=c(2,2))
  if(length(genes)==1) res<-plot_mk_in_cluster4(x=plot_gene,plot_max=data_max[1],gene_name=title_lab,title=title_lab,tsne_res=tsne,count_mx=mx[,rownames(tsne)],nl=999,is_0f=zerof,plot_cex1=plot_cex, plot_cex2=plot_cex, xlm=xx, ylm=yy)
  else res<-mapply(plot_mk_in_cluster4,x=plot_gene,plot_max=data_max,gene_name=title_lab,title=title_lab,MoreArgs=list(tsne_res=tsne,count_mx=mx[,rownames(tsne)],nl=999,is_0f=zerof,plot_cex1=plot_cex, plot_cex2=plot_cex, xlm=xx, ylm=yy))   
#  dev.off()
  if(is_ret) return(plot_gene)
}
plot_genes_on_tsne(tsne=t(b2_mono@reducedDimS), mx=mg_mx, genes=c(get_id(c('SOX3','VSX2','MITF','PMEL'))), file_name=paste('tmp',sep=''),path=paste('../result/',sep=''),plot_cex=.6,is_order=F) # marker from paper

#' Plot waves of gene expression in two branches of pseudo-axis
# A function to plot dynamic of gene expression over pseudo-axis
plot_dynamic<-function( gene, cell, file_name='_', sp=0.25, eps=10^(-10),is_ret=T ){
  if(length(sp)==1) sp<- rep(sp,length(gene))
  res<-list('')
  pdf(file_name,width=6,height=8)
  par(cex=2,las=1,mar=c(4,4,1,3),lwd=3,mfrow=c(4,3))
  for( i in 1:length(gene)){
    val<- log2(mg_mx[gene[i],cell]+1)+eps
    sm_curve<- loess.smooth( 1:length(cell), val ,span=sp[i])
    plot( sm_curve$x, sm_curve$y-eps, type='l', xlab='pseudo time',ylab='relative expression',frame=F,main=ge2an[gene[i]] )
#    text( length(cell), sm_curve$y[length(sm_curve$y)], label=ge2an[gene[i]], xpd=NA, pos=4)
   res<-c( res, list(sm_curve$y-eps))
  }
  dev.off()
  if(is_ret) return(res[-1] )
}
# RPC branch
# Cluster on gene expression
library(cluster) # use PAM clustering
b2_ret<- c( b2_cell[ b2_mono@phenoData@data[,4] == 2 ][order(-b2_mono@phenoData@data[b2_mono@phenoData@data[,4]==2,3])],
  b2_cell[ b2_mono@phenoData@data[,4] == 3 ][order(b2_mono@phenoData@data[b2_mono@phenoData@data[,4]==3,3])] )
ret_mx<- log2(mg_mx[ ret_gene , b2_ret ]+1)
ret_pam<- pam( t(apply(ret_mx,1,function(x){ x/max(x) })), k=4)
ret_pam_clu<-ret_pam$clustering
# adjust the cluster of some genes
ret_pam_clu[get_id(c('PAX6','UCHL1','S1PR3','FEZF2'))]<-3
ret_pam_clu[get_id(c('PRSS23'))]<-2
ret_pam_clu[get_id(c('CCL2'))]<-1
ret_wave<- plot_dynamic( get_id(c('SOX3','SIX6','VSX2','SIX3')), cell= b2_ret, file_name=paste('../result/ret_ref_wave.pdf',sep=''), sp=c(.55,.55,.35,.55))
ret_wave<- sapply( ret_wave, function(x){ (x-min(x))/(max(x)-min(x)) })
# Calculate the mean line of wave
ret_wave_res<- list('')
for(i in 1:4){
rm_gene<- c( list(''),list(get_id(c('CKB'))),
list(get_id(c('CCL2','PAX6','UCHL1','ID1','HES1'))),
list(''))[[i]]
raw_wave<-plot_dynamic( setdiff(names(ret_pam_clu)[ret_pam_clu==i], rm_gene), cell= b2_ret, file_name=paste('../result/wave/ret_all_wave',i,'.pdf',sep=''), sp=.3, is_ret=T)
if(i==4 | i==1) raw_wave<-plot_dynamic( setdiff(names(ret_pam_clu)[ret_pam_clu==i], rm_gene)[ sapply(raw_wave, function(x){max(x)-min(x)}) > 0.1 ], cell= b2_ret, file_name=paste('../result/wave/ret_all_wave',i,'.pdf',sep=''), sp=.3, is_ret=T)
raw_wave<-t(sapply(raw_wave, function(x){ (x-min(x))/(max(x)-min(x)) }))
raw_val<- cbind( apply(raw_wave,2,mean), apply(raw_wave,2,function(x){ sd(x)/sqrt(length(x)) } ))
ret_wave_res<- c( ret_wave_res, list(raw_val))
}
ret_wave_res<-ret_wave_res[-1]
# plot four mean+sd waves together
plot_col<- col6[c(5,10,18,7)]
#pdf('../result/wave.pdf',width=8,height=4)
par(cex=2,las=1,mar=c(2,4,2,1),lwd=3)
plot(0,0,xlim=c(1,dim(ret_wave_res[[3]])[1]),ylim=c(0,1),col='white',frame=F,xlab='',ylab='relative expression',main='',xaxt='n')
for(i in 1:4){
  polygon(x = c( 1:(dim(ret_wave_res[[i]])[1]) , rev(1:(dim(ret_wave_res[[i]])[1])) ), y = c( ret_wave_res[[i]][,1]+ret_wave_res[[i]][,2], rev(ret_wave_res[[i]][,1]-ret_wave_res[[i]][,2]) ),col = adjustcolor(plot_col[i], alpha.f = 0.5), border = NA)
  if(i!=3&i!=2) text( which.max( ret_wave[,i]), par('usr')[4]*1.05, label=c('SOX3','SIX6','VSX2','SIX3')[i], col=plot_col[i] ,xpd=NA)
  else if(i==2) text( which.max( ret_wave[,i])-5, par('usr')[4]*1.05, label=c('SOX3','SIX6','VSX2','SIX3')[i], col=plot_col[i] ,xpd=NA)
  else text( which.max( ret_wave[,i])+5, par('usr')[4]*1.05, label=c('SOX3','SIX6','VSX2','SIX3')[i], col=plot_col[i] ,xpd=NA)
}
text( (par('usr')[1]+par('usr')[2])/2, -0.15, label='pseudo time of RPC',xpd=NA)
#dev.off()

# RPE branch
b2_rpe<- c( b2_cell[ b2_mono@phenoData@data[,4] == 2 ][order(-b2_mono@phenoData@data[b2_mono@phenoData@data[,4]==2,3])],
  b2_cell[ b2_mono@phenoData@data[,4] == 1 ][order(-b2_mono@phenoData@data[b2_mono@phenoData@data[,4]==1,3])] )
# Cluster on gene expression
rpe_mx<- log2(mg_mx[ rpe_gene , b2_rpe ]+1)
rpe_pam<- pam( t(apply(rpe_mx,1,function(x){ x/max(x) })), k=4)
rpe_pam_clu<-rpe_pam$clustering
# adjust
rpe_pam_clu[get_id(c('AC009501.4','CPAMD8','TRPM3'))]<-4
rpe_pam_clu[get_id(c('ALDH1A3','LHX2'))]<-3
rpe_pam_clu[rpe_pam_clu==2]<-1
rpe_pam_clu[get_id(c('POMC','CCK','FAM19A5','GPC4','CLYBL','CRMP1','LMO3','ACP5','FABP3','MLLT11','COL2A1','SOX2','CKB','DAPL1','PAX2','HCRT'))]<-2
rpe_wave_gene<-get_id(c('SOX3','CCK','PMEL','OTX1'))
rpe_wave<- plot_dynamic( rpe_wave_gene, cell= b2_rpe, file_name=paste('../result/rpe_ref_wave.pdf',sep=''), sp=seq(.1,.9,0.05)[c(6,6,6,6)])
rpe_wave<- sapply( rpe_wave, function(x){ (x-min(x))/(max(x)-min(x)) })
# calculate the mean line of wave
rpe_wave_res<- list('')
for(i in 1:4){
rm_gene<- c( list(get_id('MAB21L1')),list(get_id(c(''))),
list(get_id(c('LHX2'))),
list(get_id('PLA2G16')))[[i]]
raw_wave<-plot_dynamic( setdiff(names(rpe_pam_clu)[rpe_pam_clu==i], rm_gene), cell= b2_rpe, file_name=paste('../result/wave/rpe_all_wave',i,'.pdf',sep=''), sp=.3, is_ret=T)
if(i==4 | i==1 | i==2) raw_wave<-plot_dynamic( setdiff(names(rpe_pam_clu)[rpe_pam_clu==i], rm_gene)[ sapply(raw_wave, function(x){max(x)-min(x)}) > 0.1 ], cell= b2_rpe, file_name=paste('../result/wave/rpe_all_wave',i,'.pdf',sep=''), sp=.3, is_ret=T)
raw_wave<-t(sapply(raw_wave, function(x){ (x-min(x))/(max(x)-min(x)) }))
raw_val<- cbind( apply(raw_wave,2,mean), apply(raw_wave,2,function(x){ sd(x)/sqrt(length(x)) } ))
rpe_wave_res<- c( rpe_wave_res, list(raw_val))
}
rpe_wave_res<-rpe_wave_res[-1]
# plot four mean+sd waves together
par(cex=2,las=1,mar=c(2,4,2,1),lwd=3)
plot(0,0,xlim=c(1,dim(rpe_wave_res[[3]])[1]),ylim=c(0,1),col='white',frame=F,xlab='',ylab='relative expression',main='',xaxt='n')
for(i in 1:4){
  polygon(x = c( 1:(dim(rpe_wave_res[[i]])[1]) , rev(1:(dim(rpe_wave_res[[i]])[1])) ), y = c( rpe_wave_res[[i]][,1]+rpe_wave_res[[i]][,2], rev(rpe_wave_res[[i]][,1]-rpe_wave_res[[i]][,2]) ),col = adjustcolor(plot_col[i], alpha.f = 0.5), border = NA)
  if(i==1) text( which.max( rpe_wave[,i]), par('usr')[4]*1.05, label=ge2an[rpe_wave_gene[i]], col=plot_col[i] ,xpd=NA)
  else if(i==2) text( which.max( rpe_wave[,i])-5, par('usr')[4]*1.05, label='SOX2', col=plot_col[i] ,xpd=NA)
  else if(i==3) text( which.max( rpe_wave[,i])-8, par('usr')[4]*1.05, label=ge2an[rpe_wave_gene[i]], col=plot_col[i] ,xpd=NA)
  else text( which.max( rpe_wave[,i])+2, par('usr')[4]*1.05, label=ge2an[rpe_wave_gene[i]], col=plot_col[i] ,xpd=NA)
}
text( (par('usr')[1]+par('usr')[2])/2, -0.15, label='pseudo time of RPE',xpd=NA)

#' Heatmap of gene expression on two branches
# RPC branch
source('../../../../distance_mx/heatmap3.r')
plot_mx<- ret_mx[order(ret_pam_clu),]
plot_mx[plot_mx>4]<-4
heatmap3( plot_mx,labRow=ge2an[rownames(plot_mx)] ,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=colorRampPalette(c("blue","white","red"))(499),color_key_label='norm UMI',color_key_label_cex=1,margins=c(0,0),color_key_axis_cex=1,key_mar=c(3, 0, 1, 1),labCol='',labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=.6, key_axis_num=2^seq(0,12,1) ,ColSideColors= time_col[sapply(colnames(plot_mx),function(x){strsplit(x,split='_')[[1]][1]})] ,ColSideColors2=type_col[cell_id[colnames(plot_mx)]], RowSideColor=col6[c(1,3,6,10)][as.factor(sort(ret_pam_clu))] )
# RPE branch
plot_mx<- rpe_mx[order(rpe_pam_clu),]
plot_mx[plot_mx>4]<-4
heatmap3( plot_mx,labRow=ge2an[rownames(plot_mx)] ,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=colorRampPalette(c("blue","white","red"))(499),color_key_label='norm UMI',color_key_label_cex=1,margins=c(0,0),color_key_axis_cex=1,key_mar=c(3, 0, 1, 1),labCol='',labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=.6, key_axis_num=2^seq(0,12,1) ,ColSideColors= time_col[sapply(colnames(plot_mx),function(x){strsplit(x,split='_')[[1]][1]})] ,ColSideColors2=type_col[cell_id[colnames(plot_mx)]], RowSideColor=col6[c(1,3,6,10)][as.factor(sort(rpe_pam_clu))] )


# rmarkdown::render('fig6.r')
