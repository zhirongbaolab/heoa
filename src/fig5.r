library(dplyr)
load('fl.rdata')
load('col6.rdata')
load('clu_ann.rdata')
load('ge2an.rdata')
load('basic.rdata')
source('functions.r')


#' [Fig.5B & Fig.S5A]  
#' UMAP of forelimb cells from each stage  

par(mai=c(0.5,0.5,0.5,0.5))
plot_cluster_on_scatter(coord = fl5_umap[[3]],point_size = 0.5,cluster_info = fl5_do[16:27],
                        col = unname(col6)[c(16,3,10,7,6,15,8,11,14,4,9,17)]) # stage3 500x500
plot_cluster_on_scatter(coord = fl5_umap[[2]],point_size = 0.5,cluster_info = fl5_do[7:15],
                        col = unname(col6)[c(3,6,4,7,11,10,9,8,15)])# stage2 470x470
plot_cluster_on_scatter(coord = fl5_umap[[1]],point_size = 1,cluster_info = fl5_do[1:6],
                        col=unname(col6)[c(7,6,3,4,10,9)]) # stage1 400x400


#' [Fig. 5D]  
#' To check the cell states of p/q (annotated as UZ) and i/j (UZ-adjoining domain) domain across stages, we project cells into a two-dimensional space.  
#' Color hexbins by domain & stages, and specific gene expressions.

load('pz_spatemp_cor.rdata') # pz_spatemp_cor; pz_id 

#' Using different color to show the domain identity of cells: cells from p/q (red), cells from i/j (blue).  
#' Using different brightness to show the stage of cells: cells from earlier stage with lighter color.  
#' Use the expression of known markers/ FGF-dependent genes as controls. 
#' Plot the pseudo-PD value of each stage to see the continuous cell states.  

# 1) domain identity: assign domain(color) value 
pz_id2=pz_id
pz_id2[pz_id2%in%c('1i','2i','3i')]=1
pz_id2[pz_id2%in%c('1c','1d','2c','2d','3c','3d')]=2


# 2) stage info: assign stage value
stage=get_stage(rownames(pz_spatemp_cor))
stage[stage=='st12']=1;stage[stage=='st13']=2;stage[stage=='st15']=3


# 3) known markers
gene=c('HOXA11','HOXA13','SPRY1','TFAP2A') %>% get_id
mx=as.matrix(fl_mx[gene,rownames(pz_spatemp_cor)]); rownames(mx)=ge2an[rownames(mx)]


# Hexbin plot
library(hexbin)
library(ggplot2)
h=hexbin(x=pz_spatemp_cor[,1], y=pz_spatemp_cor[,2], xbins=17, IDs = TRUE)
gg_input=cbind.data.frame( x=pz_spatemp_cor[,1], y=pz_spatemp_cor[,2], cID=h@cID,
                           stage=as.numeric(stage), ident=as.numeric(pz_id2),
                           t(mx)) %>% tidyr::gather(var, val, -x, -y, -cID) 
cID_ident = gg_input %>% dplyr::filter(var=='ident') %>% group_by(cID) %>% summarise(mean=mean(val))
cID_stage = gg_input %>% dplyr::filter(var=='stage') %>% group_by(cID) %>% summarise(mean=mean(val))
cID_HOXA11 = gg_input %>% dplyr::filter(var=='HOXA11') %>% group_by(cID) %>% summarise(mean=mean(val))
cID_HOXA13 = gg_input %>% dplyr::filter(var=='HOXA13') %>% group_by(cID) %>% summarise(mean=mean(val))
cID_TFAP2A = gg_input %>% dplyr::filter(var=='TFAP2A') %>% group_by(cID) %>% summarise(mean=mean(val))
cID_SPRY1 = gg_input %>% dplyr::filter(var=='SPRY1') %>% group_by(cID) %>% summarise(mean=mean(val))

gg_input2 = cbind.data.frame(cID=h@cell, x=hcell2xy(h)$x, y=hcell2xy(h)$y,
                             ident=cID_ident$mean, stage=cID_stage$mean,
                             HOXA11=log2(1+cID_HOXA11$mean),
                             HOXA13=log2(1+cID_HOXA13$mean),
                             TFAP2A=log2(1+cID_TFAP2A$mean),
                             SPRY1=log2(1+cID_SPRY1$mean))
ggplot(data=gg_input2) +
  geom_hex(aes(x=x, y=y, fill=ident,alpha=stage), col="gray70", stat="identity", size=0) + 
  scale_fill_gradient(low="red2", high="blue3") + scale_alpha_continuous(range = c(0.1, 1))+
  ggtitle("")  + xlab("") + ylab("") + 
  theme_void() + theme(legend.position="none")

ggplot(data=gg_input2) +
  geom_hex(aes(x=x, y=y, fill=HOXA11), col="grey30", stat="identity", size=0) + 
  scale_fill_gradient(low="grey99", high="brown") +
  ggtitle("")  + xlab("") + ylab("") + 
  theme_void() + theme(legend.position="none")

ggplot(data=gg_input2) +
  geom_hex(aes(x=x, y=y, fill=HOXA13), col="grey30", stat="identity", size=0) + 
  scale_fill_gradient(low="grey99", high="brown") +
  ggtitle("")  + xlab("") + ylab("") + 
  theme_void() + theme(legend.position="none")

ggplot(data=gg_input2) +
  geom_hex(aes(x=x, y=y, fill=TFAP2A), col="grey30", stat="identity", size=0) + 
  scale_fill_gradient(low="grey99", high="brown") +
  ggtitle("")  + xlab("") + ylab("") + 
  theme_void() + theme(legend.position="none")

ggplot(data=gg_input2) +
  geom_hex(aes(x=x, y=y, fill=SPRY1), col="grey30", stat="identity", size=0) + 
  scale_fill_gradient(low="grey99", high="brown") +
  ggtitle("")  + xlab("") + ylab("") + 
  theme_void() + theme(legend.position="none")


#' [Fig 5E]  
#' Only cells (domain p/q & i/j) from CS15-16 were used and the trajectory was reconstructed with diffusion map (R package destiny, k=400). The co-changing genes along the trajectory were defined in the same way as the AP-related genes in neural tube section.   
do_pheatmap <- function(genes,cells,mx=fl_mx,clusterRow=T) {
  x=mx[genes,cells]; rownames(x)=ge2an[rownames(x)]
  if(length(grep("TSPAN13",rownames(x)))==1){x['TSPAN13',]= sapply(x['TSPAN13',],function(y) ifelse(y<2,0,y))}
  x=apply(x,2,function(i) ifelse(i>10,10,i))
  info=data.frame(stage=get_stage(cells),embryo=get_embryo(cells),body=get_body_part(cells),
                  ident=get_clu_ind(fl5_do[c('32','33','21','22','26','12','13','15')])[cells])
  rownames(info)=cells
  ann_col=list(ident=c('32'='blue3', '33'='red2','21'='blue3','22'='blue3',
                       '26'='red2','12'='blue3','13'='blue3','15'='red2'),
               stage=c(st12='grey',st13='pink',st15='red4'),embryo=embryo_col,body=body_col)
  
  pheatmap(log2(1+x),cluster_cols = F,cluster_rows = clusterRow,
           annotation_col = info[,4,drop=F],annotation_colors = ann_col,annotation_legend = F,
           show_colnames = F,border_color = F,fontsize_row = 9,
           color = colorRampPalette(c('blue','white','red'))(100))
}




### [ part 1: prepare cell and gene ]
cell=unlist(fl5_do[c('32','33')])
hvg.1206=do_hvg(mx=fl_mx,cell = cell[get_embryo(cell)=='E1206'],rm.hox_gene = F)[[1]]
hvg.0809=do_hvg(mx=fl_mx,cell = cell[get_embryo(cell)=='E0809'],rm.hox_gene = F)[[1]]
hvg.1205=do_hvg(mx=fl_mx,cell = cell[get_embryo(cell)=='E1205'],rm.hox_gene = F)[[1]]
gene=names(table(c(hvg.0809,hvg.1205,hvg.1206)))[table(c(hvg.0809,hvg.1205,hvg.1206))==3]


### [ part 2: dm compute ]
library(destiny)
library(gridExtra)
library(pheatmap)
library(matrixStats)

dm=DiffusionMap(log1p(t(as.matrix(fl_mx[gene,cell]))),k = 400)
dpt=DPT(dm)
#grid.arrange(plot(dpt), plot(dpt,col_by='branch'),plot(dpt,col=get_embryo(cell)),plot(dpt,col=get_clu_ind(fl5_do)[cell]), ncol = 2)
#do_pheatmap(genes = gene,cells=cell[order(dpt[tips(dpt)[[2]], ])]) #check pseudo-PD


### [ part 3: co-PD genes ]
cell.dpt= dpt[tips(dpt)[[2]],]
input=unique(unlist(c(hvg.0809,hvg.1205,hvg.1206,get_id(c('HOXA13','HOTTIP',"KIAA1715",'TUBB6','C3orf58','SCX')))))
gene.c=apply(fl_mx[input,cell],1,function(i){cor(i,cell.dpt)});  gene.c[is.na(gene.c)]=0
g=names(gene.c)[gene.c<=(mean(gene.c)-1*sd(gene.c)) | gene.c>= (mean(gene.c)+1.2*sd(gene.c))]
g=setdiff(c(g,get_id('SCX')),get_id(c('PITX1','XXbac-BPG32J3.19','TPM1','ID2','PRRX1','DUSP6','PGF'))) #remove unwanted


# manually order genes
binary2=c('HOXA13','HOTTIP','KIAA1715','TSPAN13','HOXD13','HOXA11','LIX1','LIMCH1')
gradient2=c('DLX5','MSX1','HEY1','TFAP2A','BAMBI','LHX2','TFAP2B','MSX2','WNT5A','EPHA4','TUBB6','C3orf58',
            'FGF7','SNAI2','IFT57','SAT1','PMAIP1','LHX9','MAP1B','TGFBI','HGF','SCX','SHOX',
            'EDN3','SHOX2','IGFBP5','MEIS2','IGDCC3','CXCL12','SFRP1','FRZB','CRABP1')
g2=c(binary2,gradient2) # st3
do_pheatmap(genes = get_id(g2),cells = cell[order(dpt[tips(dpt)[[2]], ])],clusterRow = F)



#' [Fig S6]  
#' Digital in situ  

# Select genes for plotting
gene=c('EMX2','PBX1','SHH','MSX1','HAND2') %>% get_id
fl_mn=avgExp(mx = fl_mx,clu.list = fl5_do,gene = gene ) #%>% log1p
rownames(fl_mn)=ge2an[rownames(fl_mn)]
colnames(fl_mn)=fl_clu_no2do[colnames(fl_mn)]
fl_mn=apply(fl_mn,2,function(x) ifelse(x>2,2,x))


# taking "EMX2" as an example.
g='EMX2' 
max=2
par(mai=c(0.4,0.4,0.4,0.4),mfrow=c(1,3)) # 400x400


# stage1
domain=c('1a','1d','1i','1j','1p','1z')
col=get_col(fl_mn[g,domain],max = max )
plot(NA,xlim=c(0,400),ylim=c(-400,0),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
for(i in 1:6) {
  cs=read.table(paste('st1/',domain[i],'.csv',sep=''),sep=',',header = T,row.names = 1 )
  polygon(cs$X+120, -cs$Y,col=col[i],lwd=2,border = 'grey30')
}

# stage2
domain=c('2a','2d','2e','2f','2h','2i','2j','2pq','2z')
col=get_col(fl_mn[g,domain],max = max )
plot(NA,xlim=c(0,400),ylim=c(-400,0),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
for(i in 1:9) {
  cs=read.table(paste('st2/',domain[i],'.csv',sep=''),sep=',',header = T,row.names = 1 )
  polygon(cs$X, -cs$Y,col=col[i],lwd=2,border = 'grey30')
}

# stage3
domain=c('3a','3b','3c','3d','3e','3f','3h','3ij','3m1','3m2','3n1','3n2','3pq','3z')
domain.1=c('3a','3b','3c','3d','3e','3f','3h','3ij','3m','3m','3n','3n','3pq','3z')
col=get_col(fl_mn[g,domain.1],max = max )
plot(NA,xlim=c(0,450),ylim=c(-450,0),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
for(i in 1:14) {
  cs=read.table(paste('st3/',domain[i],'.csv',sep=''),sep=',',header = T,row.names = 1 )
  polygon(cs$X, -cs$Y,col=col[i],lwd=2,border = 'grey30')
}


