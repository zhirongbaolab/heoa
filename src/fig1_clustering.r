####### A single cell transcriptome atlas of human early embryogenesis

####################
## Load dataset
## On computer node LILAC2
####################

# Read raw and normalized data matrix
data_path<-'data/'
#data_path<-'../cross_embryo/result/'
load(file=paste(data_path,'raw_mx.rdata',sep='')) # raw matrix
load(file=paste(data_path,'allc185140_original_tot.rdata',sep='')) # total UMIs of each cells
load(file=paste(data_path,'norm_mx.rdata',sep=''))
raw_mx<- raw_mx
norm_mx<- norm_mx
rm( raw_mx, norm_mx)

# gene annotation
load(file=paste(data_path,'gene_anno.rdata',sep=''))
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
all_sam<- c('h0','h22','h5','h9a','h9b','ht7','l0','l21','l22','l9','lv5','lv6','t0','t21','t22','t5','t9','tv7','v0','v21','v22','v6','vl5') # all samples
hbg<-get_id(c("HBA1","HBA2","HBE1","HBG1","HBG2","HBZ")) # Signatures of erythroid
hbg_sum<- apply( raw_mx[hbg,],2,sum)
pdf('plot/hbg_umi.pdf', width=6,height=8)
par(cex=2,las=1,mar=c(6,4,1,1),mfrow=c(4,3),pch=16)
for( i in 1:23){
emb<-  all_sam[i]
plot(raw_total[ind], hbg_sum[ind], xlab='Total UMIs',ylab='Total UMIs of Hemo genes',main=emb,frame=F, cex=.3)
}
dev.off()

####################
## Control cells by total UMIs, genes, MT
####################
# Quality filtering:
# High bound: UMIs < .9 quantile; Low bound: 1000 genes
# AND MT% < 10%
all_raw<- by(t(raw_mx), sapply(colnames(raw_mx),get_sam),t)  # raw matrix by sample
cid<- lapply(raw_mx, colnames) # cell id of each sample
mtg<- names(ge2an)[grep("^MT-",ge2an)]
det_mtg<-sapply(all_raw,function(x){ colSums(x[names(ge2an) %in% mtg,]) }) # total MT
det_mtp<-mapply(function(x,y){ x/y}, x=det_mtg, y=det_read)  # percentage of MT
det_gene<-sapply(all_raw,function(x){ apply(x>0,2,sum) })
det_read<-sapply(all_raw,colSums)
flc<- mapply(function(lv,lc,hv,hc,cl, mtv, mtc){ # filter cells
  l_bad<- cl[ lv<lc ]
  h_bad<- cl[ hv>quantile(hv, hc) | mtv>mtc ]
  res<- setdiff( cl, c(l_bad, h_bad) )
  return(res)
}, lv=det_gene, lc=1000, hv=det_read, hc=.9, cl=cid, mtv=det_mtp, mtc=.1 )

####################
## Control cells by doublets ratio
####################
# run in scrublet.txt in command line
# read score back and set cutoff to 0.4 accroding to histgram
scr_res<-lapply( all_sam, function(x){
  cl<- unlist(cid[which(body==x)])
  sc<- scan( paste('../result/scrublet2022/',x,'_score.txt',sep=''), what=numeric(0))
  pr<- scan( paste('../result/scrublet2022/',x,'_predict.txt',sep=''), what=logical(0))
  print( (length(cl)==length(sc)) & (length(cl)==length(pr)) )
  res<- cbind(sc,pr)
  rownames(res)<-cl
  return(res)
})
scr_cut<- 0.4
scr_db<- unlist(sapply(scr_res, function(x){ rownames(x)[x[,1]> scr_cut ] }))
flc<- setdiff( flc, scr_db) # remove doublets before clustering

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
all_hvg<- lapply( all_sam , function(x){
  mx<- norm_mx[, get_sam(colnames(norm_mx)) %in% x]
  hvg<-setdiff(get_hvg(mx_norm=mx,path1='plot/',part=x,conf_level=0.99), rm_gene)  # based on Poisson distribution
  return(hvg)
})

####################
## Identification of developmental system by TFs in HVGs (level-1)
####################
tf<-scan('data/humanTF_GO0003700_clear.txt',what=character(0)) # from Emsembl
hvg_tf<- intersect(tf, unlist(all_hvg)) # 288
write(hvg_tf, file='data/tf288.txt') # output TFs in HVGs
library(Seurat)
features<- hvg_tf
print(paste('features:',length(features)))
seu<- lapply( flc_mx, function(x){ CreateSeuratObject( x, min.cells = 0, min.features = 0 )}) # 'flc_mx': filtered matrix by batch
seu <- lapply(X = seu, FUN = function(x) {
    x <- NormalizeData(x)
})
seu <- lapply(X = seu, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchor <- FindIntegrationAnchors(object.list = seu, anchor.features = features, reduction = "rpca", reference=1, dims=1:50 ) # batch correct
seu_comb <- IntegrateData(anchorset = anchor)
DefaultAssay(seu_comb) <- "integrated"
# Run the standard workflow for visualization and clustering
seu_comb <- ScaleData(seu_comb, verbose = FALSE)
seu_comb <- RunPCA(seu_comb, npcs = 50, verbose = FALSE, features=features )
seu_comb <- RunUMAP(seu_comb, reduction = "pca", dims = 1:30, n.neighbors = 20, min.dist = 0.1)
umap<-seu_comb@reductions$umap@cell.embeddings
# get clusters
resos<-seq(.1,1,.1)
clu_res<- numeric(0)
seu_comb<-FindNeighbors(seu_comb, dims = 1:30)
for(reso in resos){
seu_comb<-FindClusters(seu_comb, dims.use =1:30, k.param = 30, resolution=reso)
seuclu<-as.numeric(seu_comb@active.ident)
clu_cen<- cbind( tapply(umap[,1], seuclu, mean), tapply(umap[,2], seuclu, mean) )
pdf(paste('../result/filter2022/',batch,'_umap_reso',reso,'.pdf',sep=''),width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
plot(umap[rand_ind,1],umap[rand_ind,2],pch=16,xlab='umap-1',ylab='umap-2',frame=F,col=col6[-1][as.factor(seuclu)][rand_ind],cex=0.1,main=length(unique(seuclu)))
text(clu_cen[,1], clu_cen[,2], label= 1:length(unique(seuclu)), cex=.8)
dev.off()
clu_res<- cbind( clu_res, seuclu )
print(paste('done with reso',reso))
}
colnames(clu_res)<- resos
rownames(clu_res)<- rownames(umap)
ds_list<-tapply( clu_res[,5], clu_res[,5], function(x){x}) # list of developmental system

####################
## Umap visualization
####################

# load color scheme
load('data/part_col.rdata')
emb_id<-names(part_col)

# Total Umap
rand_ind<-sample.int(dim(umap)[1])
pdf('plot/total_umap.pdf',width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
plot(umap[rand_ind,1],umap[rand_ind,2],pch=16,xlab='umap-1',ylab='umap-2',frame=F,col=part_col[ get_sam(rownames(umap)) ][rand_ind],cex=0.1)
dev.off()

# plot gene expression on Umap
plot_genes_on_tsne(tsne=umap, mx=norm_mx, genes=c(get_id(c('PAX6','MITF'))), file_name='test',path='plot/',plot_cex=.2,is_order=F) # test
plot_genes_on_tsne(tsne=umap, mx=norm_mx, genes=hvg_tf, file_name='test',path='plot/',plot_cex=.2,is_order=T)


####################
## Clustering in each developmental system (level-2)
####################
# remove hemoglobin genes, MT genes, sex-specific genes, cell cycle genes, and batch-effect genes during clustering
hoxg<- names(grep('^HOX',ge2an,value=T))
# clustering
r2_res<-mapply( function(x,y){
  res<-r2_wrap( cell=x, mx=norm_mx[,x], remove_gene=rm_gene, path='plot/', return_r1=y ) 
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
## Preparation for the calculation of differentially expressed genes (DEGs)
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
#################### continue in ficlu.Rmd