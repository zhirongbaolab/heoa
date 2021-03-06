####### A single cell transcriptome atlas of human early embryogenesis

####################
## Load dataset
####################

# Read raw and normalized data matrix
data_path<-'data/'
#data_path<-'../cross_embryo/result/'
load(file=paste(data_path,'all_emb/raw_mx_c333.rdata',sep=''))
load(file=paste(data_path,'all_emb/raw_mx_c333_total.rdata',sep='')) # total UMIs of each cells
load(file=paste(data_path,'all_emb/norm_mx_c333.rdata',sep=''))
raw_mx<- raw_mx_c333
norm_mx<- norm_mx_c333
rm( raw_mx_c333, norm_mx_c333)

# gene annotation
load(file=paste(data_path,'all_emb/gene_anno.rdata',sep=''))
ge2an<-as.character(gene_anno[,2])
names(ge2an)<- gene_anno[,1]

####################
## Load functions
####################

source('function.r') # all functions used in iterative clustering


####################
# Preprocessing
####################

# Red blood cells have much fewer expressed genes than other cells. To avoid the potential distortion on clustering by the large population of red blood cells, we identified red blood cells by the expression of HBA1, excluded them from clustering.

# Cells with extremely high expression of hemoglobin were considered as erythroids
hbg<-get_id(c("HBA1","HBA2","HBE1","HBG1","HBG2","HBZ")) # Signatures of erythroid
hbg_sum<- apply( raw_mx[hbg,],2,sum)
pdf('plot/hbg_umi.pdf', width=6,height=8)
par(cex=2,las=1,mar=c(6,4,1,1),mfrow=c(4,3),pch=16)
for( i in 1:12){
emb<-  c('ht7', 'tv7', 'h0', 'l0', 't0', 'v0', 'h5', 'lv5', 't5', 'vl5', 'lv6', 'v6') [i]
plot(raw_total[ind], hbg_sum[ind], xlab='Total UMIs',ylab='Total UMIs of Hemo genes',main=emb,frame=F, cex=.3)
}
dev.off()


####################
## HVG detection
####################

# prepare genes to be removed: hemoglobin genes, MT genes, sex-specific genes, cell cycle genes, and batch-effect genes.
cc_gene<- scan('list/total_cc.txt',what=character(0)) # cell cycle genes, merged from: pubmed ID 25378319; pubmed ID 30452682
hb_gene<- scan('list/hb_gene.txt',what=character(0)) # hemoglobin genes
mt_gene<-names(ge2an)[grepl(pattern='^MT-', x= ge2an)] # MT genes
bt_gene<- unique(get_id(scan('list/batch.txt',what=character(0))))  # batch-effect genes (including sex-specific genes), merged from: pubmed ID 30096314; pubmed ID 31835037
rm_gene<- unique(c(cc_gene, hb_gene, mt_gene, bt_gene))

# Calculate HVGs by sample 
emb<- c('ht7', 'tv7', 'h0', 'l0', 't0', 'v0', 'h5', 'lv5', 't5', 'vl5','h9a','h9b','t9','l9', 'lv6', 'v6')
all_hvg<- lapply( emb , function(x){
  mx<- norm_mx[, get_sam(colnames(norm_mx)) %in% x]
  hvg<-setdiff(get_hvg(mx_norm=mx,path1='plot/',part=x,conf_level=0.99), rm_gene)  # based on Poisson distribution
  return(hvg)
})


####################
## Umap visualization
####################

# load developmental systems
load('data/ds_list.rdata') # list of cells in each developmental system
# load TFs for visualization
hvg_tf<-scan('data/tf291.txt',what=character(0))
# downsample cell types with too many cells
ds<-scan('data/downsample.txt',what=character(0))
# load color scheme
load('data/part_col.rdata')
emb_id<-names(part_col)

# Total Umap
library(Seurat)
seu <- CreateSeuratObject(norm_mx[,setdiff(colnames(norm_mx),ds)] , min.cells = 0, min.features = 0)
seu <- NormalizeData(object = seu)
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs=50,features = hvg_tf, verbose = F)
seu <- RunUMAP(object = seu, dims = c(1:50), reduction='pca', min.dist=.1, n.neighbors=20)
umap<-seu@reductions$umap@cell.embeddings
rand_ind<-sample.int(dim(umap)[1])
pdf('plot/total_umap.pdf',width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
plot(umap[rand_ind,1],umap[rand_ind,2],pch=16,xlab='umap-1',ylab='umap-2',frame=F,col=part_col[ get_sam(rownames(umap)) ][rand_ind],cex=0.1)
dev.off()

# plot gene expression on Umap
plot_genes_on_tsne(tsne=umap, mx=norm_mx, genes=c(get_id(c('PAX6','MITF'))), file_name='test',path='plot/',plot_cex=.2,is_order=F) # test
plot_genes_on_tsne(tsne=umap, mx=norm_mx, genes=hvg_tf, file_name='test',path='plot/',plot_cex=.2,is_order=T)


####################
## Clustering of developmental systems
####################
# remove hemoglobin genes, MT genes, sex-specific genes, cell cycle genes, and batch-effect genes during clustering
# also remove HOX genes during clustering
hoxg<- names(grep('^HOX',ge2an,value=T))
# clustering
r2_res<-mapply( function(x,y){
  res<-r2_wrap( cell=x, mx=norm_mx[,x], remove_gene=c(rm_gene,hoxg), path='plot/', return_r1=y ) 
  return(res)
}, x=ds_list[setdiff(1:14,9)], y=c(F,F,T,T,F,F,T,T,T,T,F,F,T), SIMPLIFY=F)
# Limb was first separated by stage and clustering was performed on each stage. See Figure 5.
# y: the number of rounds of clustering. F: 2 rounds; T: 1 round.

# format clusters in a linear list
all_clu<- do.call( 'c', mapply( function(x,y){
  linear_r2_clu( x, name=y)[[1]]
}, x=r2_res, y=names(r2_res), SIMPLIFY=F) )


####################
## Calculate signature genes of clusters
####################
# To annotate each cluster, signature genes of clusters were calculated by Seurat.
clu_sig<- findMarkerbySeurat(norm_mx[,unlist(all_clu)],cluster.info=all_clu , clusterInfo_format='list', top_n=20)


####################
## Differentially expressed genes (DEGs)
####################
# Calculate mean expression and detected fraction of each gene in each cluster
tp_mn<- sapply( all_clu, function(x){
  apply( norm_mx[, x], 1, mean )
})
tp_fr<- sapply( all_clu, function(x){
  apply( norm_mx[, x], 1, function(y){ sum(y>0)/length(y) } )
})
tp_zs<- t(apply( tp_mn,1, function(x){
  if(sum(x>0)==0) return( rep(0,length(x)) )
  else return((x-mean(x))/sd(x))
}))
colnames(tp_zs)<-colnames(tp_mn)
# Only consider genes that are expressed in >=1 cluster as the candidate of DEGs
in_gene<- rownames(norm_mx)[ apply( (tp_mn[,]>=0.5&tp_fr[,]>=.4),1,sum )>0 ]
in_gene<-setdiff( in_gene, rm_gene)

# Method 1: general DEGs that have large variance on expression across cell types (not assigned to any cell type)
deg1<- in_gene[ apply(abs(tp_zs[in_gene,]),1, function(x){ x[x<4]<-0;sum(x) })>=15 | apply(abs(tp_zs[in_gene,]),1, function(x){ sum(x>7) })>=1 ]

# Method 2: confident DEGs that are not expressed in enough number of cell types to be convincingly called as a DEG (assigned to cell type)
low_cut<-0.3
mb_cut<- 70
length(deg2<- in_gene[ apply(tp_mn[in_gene,]<low_cut, 1, sum)>= mb_cut ]) # enough number of cell types without expression
deg2_type<- lapply( deg2, function(x){
  d<- tp_mn[x,]
  zs<- (d - mean(d[d<.5]))/sd(d[d<.5])
  return(colnames(tp_mn)[tp_mn[x,]>=0.5 & tp_fr[x,]>=.4 & zs>=7])
}) # For each DEG, call cell types that expressed ths DEG.
