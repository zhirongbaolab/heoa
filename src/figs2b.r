#' Fig. S2b: pairwise correlation between 333 cell types
#' Pearson's correlation based on mean expression of DEGs in each cell type

#' Load gene expression data and cluster information
library(Matrix)
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
load('../../../cross_embryo/code/recluster_deg4_deg_0312.rdata') # DEGs
load(file='../../../cross_embryo/result/all_emb/norm_mx_c333.rdata') # normalized matrix
ficlu<-read.table('~/BaoLab/xuy/human/cross_embryo/result/reclustering3/TJ_clu/cell_annotate.csv',sep=',',as.is=T,header=T) # cell list for each cell type
tj_clu<-read.table( '../../../cross_embryo/code/tj_tables1_sheet1.txt', sep='\t',as.is=T, header=T) # annotation of cell type (Table S1)
c333_od<- tj_clu[,2]
c333_od_chunk<- tj_clu[,5] # developmental systems for each cell type
names(c333_od_chunk)<- c333_od
tmp<-unique(read.table('~/baolab/xuy/ciona/cluster2/product/code/type_all_col3.txt',sep='\t',as.is=T,comment.char='')[,2])
col6<-c('lightgrey','red',tmp[3],'purple',tmp[2],tmp[5],'blue','green4','orchid','turquoise','sienna',tmp[9],'yellow2','hotpink','navy','steelblue','skyblue','pink','black',tmp[4],rainbow(7))
tj_col<-col6[c(3:4,16,6:15,18,17,2,25,20,5)] # color codes of developmental systems
names(tj_col)<- unique(tj_clu[,5])

#' Calculate mean expression of DEGs in each cell type
c333_deg2_mn<- sapply( c333_od, function(x){
  clu<-tj_clu[tj_clu[,2]==x,3]	       
  cell<-ficlu[ficlu[,2]==clu,1]
  res<-apply( norm_mx_c333[all_deg[[3]],cell],1,mean )
  return(res)
})
c333_deg2_corpe<- cor( c333_deg2_mn, method='pearson')

#' Plot correlation as a heatmap
c333_ex<- read.table( 'c333_od_for_expand.txt', sep='\t',as.is=T)[,1] # duplicate cell types in small developmental system to make them more visable
plot_mx<-c333_deg2_corpe[c333_ex,c333_ex]
source('../../../../distance_mx/heatmap3.r')
heatmap3( plot_mx,labRow=NA,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=colorRampPalette(c("blue","white","red"))(499),color_key_label='cor',color_key_label_cex=1,margins=c(0,0),color_key_axis_cex=1,key_mar=c(3, 0, 1, 1),labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=1,cexCol=.3, ColSideColors= tj_col[ c333_od_chunk[rownames(plot_mx)] ], RowSideColors=tj_col[ c333_od_chunk[rownames(plot_mx)] ] ,labCol_pos=c(1,3), labCol=NA)

# rmarkdown::render('figs2b.r')