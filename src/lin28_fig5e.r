#' Punctuated embryonic stages in vertebrates.
#' Hierarchical clustering of similarity between timepoints revealed four stages of gene expression with clear temporal boundaries. Similarity of gene expression between timepoints within each species by calculating the pairwise correlation of transcription factors (TFs) expression between timepoints.

# Zebrafish
source('../../../../distance_mx/heatmap3.r')
library(dendextend)
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
knitr::opts_chunk$set(warning = FALSE)
load(file='../../White/code/dr_all_exp.rdata')
go_tf<- read.table('../../White/ensembl_GO0003700.txt',sep='\t',as.is=T,header=T) # zebrafish TFs
rownames(go_tf)<-go_tf[,1]
tf<- intersect( rownames(dr_all_exp), go_tf[,1])
tf_exp<- dr_all_exp[tf,]
self_tf_cor<- cor(tf_exp[apply(tf_exp,1,max)>10,])
plot_mx<- self_tf_cor
st_anno<- c(paste( c(0,0.75,2.25,3,4.3,5.25,6,8,10.3,16,19,24,30,36),'hpf',sep=''),'2d','3d','4d','5d')
heatmap3( plot_mx,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=c(colorRampPalette(c('blue',"white","red"))(499)),color_key_label='cor',color_key_label_cex=1.3,margins=c(0,0),color_key_axis_cex=1.3,key_mar=c(3, 0, 1, 1),labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=1.2,cexCol=1, labRow=st_anno,labCol=NA)
# Hierarchical clustering and dendrogram
tmp<-unique(read.table('~/baolab/xuy/ciona/cluster2/product/code/type_all_col3.txt',sep='\t',as.is=T,comment.char='')[,2])
col6<-c('lightgrey','red',tmp[3],'purple',tmp[2],tmp[5],'blue','green4','orchid','turquoise','sienna',tmp[9],'yellow2','hotpink','navy','steelblue','skyblue','pink','black',tmp[4],rainbow(7))
plot_mx<- self_tf_cor
rownames(plot_mx)<-letters[1:nrow(plot_mx)]->colnames(plot_mx)
hc<- hclust(dist(plot_mx),method='ward.D2')
dd <- as.dendrogram(hc)
dd.reorder <- sort(dd)
dd.reorder<-set(dd.reorder,"branches_k_color",  k = 4, value=col6[c(2:4,6)])
# Put this dendrogram on top of heatmap
par(las=1,mar=c(0,0,0,0),lwd=4)
plot(dd.reorder, yaxt='n')

# Frog
load('xt_all_exp.rdata')
go_tf<- read.table('../ensembl_GO0003700.txt',sep='\t',as.is=T,header=T) # frog TFs
go_tf[go_tf[,2]=='',2]<-go_tf[go_tf[,2]=='',1]
rownames(go_tf)<-go_tf[,1]
tf<- intersect( rownames(xt_all_exp), go_tf[,1])
tf_exp<- xt_all_exp[tf,]
self_tf_cor<- cor(tf_exp[apply(tf_exp,1,max)>10,])
plot_mx<- self_tf_cor
#pdf('../result/xt_self_tf_cor.pdf',width=7,height=7)
heatmap3( plot_mx,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=c(colorRampPalette(c('blue',"white","red"))(499)),color_key_label='cor',color_key_label_cex=1.3,margins=c(0,0),color_key_axis_cex=1.3,key_mar=c(3, 0, 1, 1),labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=1.2,cexCol=1, labRow=colnames(xt_all_exp),labCol=NA)
#dev.off()
# Hierarchical clustering and dendrogram
plot_mx<- self_tf_cor
rownames(plot_mx)<-letters[1:nrow(plot_mx)]->colnames(plot_mx)
hc<- hclust(dist(plot_mx),method='ward.D2')
dd <- as.dendrogram(hc)
dd.reorder <- sort(dd)
dd.reorder<-set(dd.reorder,"branches_k_color",  k = 4, value=col6[c(2:4,6)])
# Put this dendrogram on top of heatmap
par(las=1,mar=c(0,0,0,0),lwd=4)
plot(dd.reorder, yaxt='n')

#' Correlation between zebrafish and frog by homologous TFs
hs_orth<-read.table('../../orth_hs_xt_dr_mm.txt',sep='\t',as.is=T,header=T) # ortholog from Ensembl
length(xt_exp<- rownames(xt_all_exp)[apply(xt_all_exp,1,max)>10])
length(dr_exp<- rownames(dr_all_exp)[apply(dr_all_exp,1,max)>10])
dr_st_anno<- scan('../../White/code/dr_st_anno.txt',what=character(0))
xt_dr_orth<- hs_orth[hs_orth[,3]!=''&hs_orth[,5]!='',3:6]
xt_dr_orth<- xt_dr_orth[ !duplicated(apply(xt_dr_orth[,c(1,3)],1,paste,collapse=' ')), ]
xd_uni_xt<- names(table(xt_dr_orth[,1]))[table(xt_dr_orth[,1])==1] # XT genes with unique DR orthlog
xd_uni_dr<- names(table(xt_dr_orth[,3]))[table(xt_dr_orth[,3])==1]
xd_uni_exp<- xt_dr_orth[xt_dr_orth[,1] %in% xd_uni_xt & xt_dr_orth[,3] %in% xd_uni_dr,]
dim(xd_uni_exp<- xd_uni_exp[ (xd_uni_exp[,1]%in% xt_exp | xd_uni_exp[,3] %in% dr_exp) & xd_uni_exp[,1] %in% rownames(xt_all_exp) & xd_uni_exp[,3] %in% rownames(dr_all_exp), ])
xd_uni_tf_exp_cor<- apply( xt_all_exp[xd_uni_exp[ xd_uni_exp[,1]%in% tf,1],], 2, function(x){ 
  apply( dr_all_exp[xd_uni_exp[xd_uni_exp[,1]%in% tf,3],], 2, function(y){ cor(x,y,method='spearman') } )
})
plot_mx<-xd_uni_tf_exp_cor
xt_st<- c(rep(1,7), rep(2,3), rep(3,7), rep(4,6))
dr_st<- c(rep(1,5), rep(2,3), rep(3,7), rep(4,3))
heatmap3( plot_mx,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=colorRampPalette(c("blue","white","red"))(499),color_key_label='Cor',color_key_label_cex=1,margins=c(0,0),color_key_axis_cex=1,key_mar=c(3, 0, 1, 1),labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=1,cexCol=1 ,labRow=NA,labCol=NA,ColSideColors=col6[c(2:4,6)][xt_st], RowSideColors=col6[c(2:4,6)][dr_st] )
# Dynamic time warping to align between zebrafish vs. frog
library(dtw)
xd_dtw<-dtw(x=1-xd_uni_tf_exp_cor)
par(cex=2,las=1,mar=c(1,1,1,1),lwd=5,pch=16)
plot(xd_dtw$index2,18-xd_dtw$index1,frame=F,type='l',xaxt='n',yaxt='n',xlab='XT',ylab='DR', lty=1) # Put the alignment on top of heatmap


# Mouse: Boroviak et al. 2015, Pijuan-Sala et al. 2019, Cao et al. 2019
# Because mouse data are from three different sources, we do not calculate correlation between timepoints from different datasets and only perform clustering withing each dataset
load('../../Boroviak_mouse/code/borov.rdata')
length(tf<-unique(read.table('../../../Imai_table/mouseTF_GO0003700.txt',sep='\t',header=T,as.is=T)[,1])) # mouse TFs
length(comg<-intersect(intersect( rownames(m2_ganno), rownames(m1_all_exp)), m3g )) # genes avaiable in all three datasets
rownames(m2_all_exp)<- rownames(m2_ganno)
dim(com_exp<- cbind(m1_all_exp[comg,], m2_all_exp[comg,-10], m3_all_exp[comg,])) # merged dataset of 3 mouse datasets
colnames(com_exp)[14:18]<- paste('E', colnames(com_exp)[14:18], sep='')
# expressed genes in at least 1 dataset of mouse
com_expg<- unique(c(comg[apply(m1_all_exp[comg,],1,max)>10], comg[apply(m2_all_exp[comg,],1,max)>0.1], comg[apply(m3_all_exp[ comg,],1,max)>0.01]))
self_tf_cor<- cor(com_exp[intersect(tf,com_expg),],method='spearman')
# Remove inter-dataset correlation and plot heatmap
plot_mx<-self_tf_cor
plot_min<- .4
plot_mx[c(1:4),-(1:4)]<- plot_min
plot_mx[c(5:13),-(5:13)]<-  plot_min
plot_mx[c(14:18),-(14:18)]<- plot_min
heatmap3( plot_mx,scale='none',dendrogram='none',trace='none',Rowv=F,Colv=F,symkey=F,density.info="none",keysize=1,col=c(col6[1],colorRampPalette(c('blue',"white","red"))(499)[1:420],'red' ),color_key_label='cor',color_key_label_cex=1.3,margins=c(0,0),color_key_axis_cex=1.3,key_mar=c(3,0,1,1),labRow_pos=c(2),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=1.2,cexCol=1,labCol=NA)


# rmarkdown::render('lin28_fig5e.r')