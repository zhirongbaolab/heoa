# functions used to analyse the single-cell RNA-seq dataset in "human atlas paper"

plot_cluster_on_scatter <- function(coord,cluster_info,point_size=0.5,col=NULL,main='',legend=F,border=F,label_clusters=F,labels=NULL,label_every_cell=F,randomPt=T){
  # coord is (x,y) position of each cell
  # cluster_info should be 'list' or 'vector' format
  
  library(matrixStats);library(crayon)
  cell=rownames(coord)
  
  
  if(legend){par(mai=c(0.5,0.5,0.5,1.5),xpd=TRUE)} else {par(mai=c(0.5,0.5,0.5,0.5))}
  if(border) {
    plot(coord,type = 'n',main=main,bty='n')
  } else {
    plot(coord,type='n',frame=F,xaxt='n',yaxt='n',xlab='',ylab='',main=main)
  }
  
  if(class(cluster_info)=="list") {
    cluster_info=cluster_info
  } else {
    cluster_info=tapply(names(cluster_info),cluster_info, function(x) x ) 
  }
  
  if(length(unlist(cluster_info))!=nrow(coord)) {
    cat(red('Warning: \n[coord] has'), nrow(coord), red("cells \n[cluster_info] has"), length(unlist(cluster_info)), red('cells'), "\n")
    cluster_info=lapply(cluster_info,function(x) intersect(x,cell))
    cluster_info=cluster_info[sapply(cluster_info,length)>0]
  }
  
  
  if(is.null(col)){
    col=rainbow(length(cluster_info))
    } else if (!is.null(col) & is.null(names(col))){
      col=col
    } else {col=col[names(cluster_info)]}
  
  if(randomPt){
    rdm=cell[sample(length(cell))]
    color=rep(NA,length(rdm)); for(i in 1:length(cluster_info)){ color[rdm%in%cluster_info[[i]]]=col[i] }
    points(coord[rdm,],col=color,pch=16,cex=point_size) 
  } else {
    for(i in 1:length(cluster_info)){
      mx.sub=coord[cell%in%cluster_info[[i]],,drop=F]
      points(mx.sub,col=col[i],pch=16,cex=point_size)
    } 
  }
  
  
  if(label_every_cell) { text(coord,col='black',labels = rownames(coord),cex=0.2) }
  if(legend) { legend('right',inset=c(-0.3,0),col=col[1:length(cluster_info)],legend=names(cluster_info),pch=16,bty='n') }
  if(label_clusters) { 
    for(i in 1:length(cluster_info)) {
      mx.sub=coord[cell%in%cluster_info[[i]],,drop=F]
      mx.sub.mean=colMedians(mx.sub)
      if(is.null(labels)) {labels=names(cluster_info)}
      text(x=mx.sub.mean[1],y =mx.sub.mean[2] ,col='black',labels = labels[i])
    }
  }
}

get_embryo <- function(cellname=cells){
  embryo=sapply(cellname, function(x){
    if(substr(x,start=1,stop=3) %in% c('ht7','tv7')) return('E1107')
    else if(substr(x,start=1,stop=3) %in% c('h0_','t0_','l0_','v0_')) return('E0710')
    else if(substr(x,start=1,stop=3) %in% c('h9a','h9b','t9_','l9_','v9_'))return('E0809')
    else if(substr(x,start=1,stop=3) %in% c('h5_','t5_','lv5','vl5')) return('E1205')
    else if(substr(x,start=1,stop=3) %in% c('h6_','t6_','v6_','lv6')) return('E1206')
  })
  return(unlist(embryo))
}

get_id <- function(x) {unlist(sapply(x, function(y){names(ge2an)[ge2an==y]}))}

get_stage <- function(cellname=cells){
  stage=sapply(cellname, function(x){
    if(substr(x,start = 1,stop = 3) %in% c('h0_','t0_','v0_','l0_')) return("st13")
    else if(substr(x,start = 1,stop = 3) %in% c('ht7','tv7')) return("st12")
    else if(substr(x,start = 1,stop = 3) %in% c('h9a','h9b','t9_','l9_','v9_','h5_','t5_','lv5','vl5','h6_','t6_','v6_','lv6')) return("st15")
  })
  return(unlist(stage))
}

do_hvg <- function(mx,cell,cut_num=50,rm.gene=NULL,rm.hox_gene=T,rm.broad_sam=T){
  
  hvg=find_hvg_by_emb(cell, out_path='./',cut_num=cut_num,xmin=-2.5, adj_int=0,all_mx = mx[,cell])
  
  hvg2=hvg[apply(mx[hvg,cell],1,function(x){sum(x>0)/length(x)})<0.7];length(hvg2) #remove HVG >= 70%
  hvg2=setdiff(hvg2,lysis_gene); length(hvg2);  #remove lysis_gene
  hvg2=setdiff(hvg2,sex_gene); length(hvg2) #remove sex_gene
  if(rm.hox_gene) {hvg2=setdiff(hvg2,hox_gene)} #remove hox_gene
  if(length(rm.gene)>0) {hvg2=setdiff(hvg2,rm.gene)} # remove particular gene
  if(rm.broad_sam) {broad_sam=get_broad_by_sam(x = cell,g = hvg,mx = mx); hvg2=setdiff(hvg2,broad_sam)} # remove sam broad
  sum(hvg2 %in% tf) #no of tf
  
  return(list(hvg=hvg,hvg2=unique(hvg2)))
}

####### A single cell transcriptome atlas of human early embryogenesis

####### Custom functions used in analysis

# get gene symbol by ID for a group of genes
get_id<-function(x){ sapply(x, function(y){names(ge2an)[ge2an==y][1]})}
get_mm_id <- function(x) {unlist(sapply(x, function(y){names(mm.ge2an)[mm.ge2an==y]}))}

# get embryo or sample by cell name
get_emb<-function(x){
  emb<-gsub("[A-z]", "", sapply(x,function(y){ strsplit(y,split='_')[[1]][1] }) )
  return(emb)
}
get_sam<-function(x){
  return(sapply(x,function(y){ strsplit(y,split='_')[[1]][1] }))
}

# Calculate HVGs by Poisson distribution (ref: https://www.biorxiv.org/content/10.1101/2020.03.02.966440v1.full)
library(DDoutlier)
library(dbscan)
library(Seurat)
# Parameters that may need adjustment
# bou_int: interval for fit a boundary for threshold in Method 1. Increase or decrease it to make the threshold fit the curve.
# path1: output directory
# part: output name
# xtem: minimal average expression for genes to be considered
get_hvg<-function(mx_norm=mx_norm,is_plot=T,cut_agg=0,bou_int=0,path1,part,xtem=100,conf_level=0.99,add_pt='',title='',is_adj=F,adj_info, add_pt2='',tran_plot=F, tran_xlm=25){
# gene measurement
if(!is_adj) ginfo<-cbind(apply(mx_norm,1,mean),apply(mx_norm==0,1,sum),apply(mx_norm,1,mean),apply(mx_norm,1,var)/(apply(mx_norm,1,mean)^2)) # mean, count-0, mean on norm, CV2 on norm
else ginfo<- cbind( adj_info[,1], adj_info[,2], adj_info[,1], apply(mx_norm[rownames(adj_info),],1,var)/(apply(mx_norm[rownames(adj_info),],1,mean)^2) )

# HVG of method 1
mean_cut<-10000
x4fit<- log2(ginfo[ ginfo[,1]< mean_cut & ginfo[,1]>0 &ginfo[,2]>0 ,1] ) # log2 on mean
y4fit<- ginfo[ ginfo[,1]< mean_cut & ginfo[,1]>0 &ginfo[,2]>0 ,2]

# define boundary points
if(bou_int==0) bou_int<- dim(mx_norm)[2]/200
if(bou_int<1) bou_int<-1
num0_agg<-KNN_AGG(cbind(x4fit,y4fit),k_min=3,k_max=6)
if(cut_agg==0) cut_agg<- sort(num0_agg)[round(quantile(1:length(x4fit),probs=0.99))]
if(cut_agg<10) cut_agg<-10
print(paste('cut_agg:',cut_agg,'bou_int:',bou_int))
#plot(1:length(num0_agg),sort(num0_agg)) # find the elbow
inlier<-names(x4fit)[num0_agg<=cut_agg]
#bou<-tapply(x4fit[inlier],y4fit[inlier],function(x){names(x)[which.max(x)] })
# instead of take a max(x) on every y, take a x on every 'bou_int' y
bou<-unlist(mapply(function(y1,y2){
  x<-x4fit[inlier]
  y<-y4fit[inlier]
  xx<-x[ y1<= y & y < y2 ]
  return(names(xx)[which.max(xx)])
},y1=seq(1,dim(mx_norm)[2],bou_int) ,y2=c(seq(1,dim(mx_norm)[2],bou_int)[-1],dim(mx_norm)[2]) ))

# smoothed y
smo_curve<-loess.smooth(x4fit[bou],y4fit[bou],span=0.25)
smooth_y<-sapply(x4fit,function(x,smo=smo_curve){
  sx<-smo$x
  sy<-smo$y
  diff_x<-abs(x-sx)
  xx<- sx[order(diff_x)[1:2]]
  yy<- sy[order(diff_x)[1:2]]
  a<-(yy[1]-yy[2])/(xx[1]-xx[2])
  b<-yy[1]-a*xx[1]
  y<-a*x+b
  return(y)
})

# terminal condition
ytem<-dim(mx_norm)[2]*0.1 # 0-count > 10%
if(xtem==100) xtem<-log2( 10*4/ (dim(mx_norm)[2]) ) # minimal cell population 10 at the level of 4 counts (similar level with ciona)
hvg1<-names(y4fit)[y4fit>smooth_y] # right side of smoothed curve
print(paste('hvg1:',length(hvg1 <- hvg1[x4fit[hvg1]> xtem & y4fit[hvg1]> ytem ])))

if(is_plot){
print(paste('print in:',paste(path1,part,'_method1_hvg.pdf',sep='')))
pdf(paste(path1,part,'_method1_hvg.pdf',sep=''),height=12,width=12)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3)
plot(x4fit,y4fit,ylab='number of 0-count cells',xlab='log2 mean of count',frame=F,cex=0.2,pch=20,col='grey',xaxt='n',yaxt='n',main=title)
#points(x4fit[bou],y4fit[bou],col='green',cex=0.2)
lines(loess.smooth(x4fit[bou],y4fit[bou],span=0.25)$x,loess.smooth(x4fit[bou],y4fit[bou],span=0.25)$y,col='red',lwd=3)
points(x4fit[hvg1],y4fit[hvg1],col='black',cex=0.2)
if(add_pt[1]!='') points(x4fit[add_pt],y4fit[add_pt],col='blue',cex=0.5)
if(add_pt2[1]!='') points(x4fit[add_pt2],y4fit[add_pt2],col='red',cex=0.5)
axis(2,tck=-0.02,lwd=3)
axis(1,tck=-0.02,lwd=3)
lines(c(xtem,xtem),c(0,dim(mx_norm)[2]),col='green3',lty=2,lwd=3)
lines(c(-10,5),c(ytem,ytem),col='green3',lty=2,lwd=3)
dev.off()
}

if(tran_plot){
pdf(paste(path1,part,'_method1_b.pdf',sep=''),height=12,width=12)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3)
plot( 2^(x4fit),log2(y4fit),ylab='log2 number of 0-count cells',xlab='mean of count',frame=F,cex=0.2,pch=20,col='grey',xaxt='n',yaxt='n',main=title, xlim=c(0,tran_xlm) )
points( (2^x4fit)[hvg1],log2(y4fit)[hvg1],col='black',cex=0.2)
if(add_pt[1]!='') points( (2^x4fit)[add_pt], log2(y4fit)[add_pt],col='blue',cex=0.5)
if(add_pt2[1]!='') points( (2^x4fit)[add_pt2], log2(y4fit)[add_pt2],col='red',cex=0.5)
axis(2,tck=-0.02,lwd=3)
axis(1,tck=-0.02,lwd=3)
dev.off()
}


# HVG of method 2
x4fit_o2<-log2(ginfo[ ginfo[,3]>0 & ginfo[,3]<mean_cut,3])
y4fit_o2<-log2(ginfo[ ginfo[,3]>0 & ginfo[,3]<mean_cut,4])
fit<-lm(y4fit_o2 ~ x4fit_o2)
x1<- seq(min(x4fit_o2),max(x4fit_o2),length.out=100)
y1<-predict(fit,newdata=data.frame(x4fit_o2=x1),interval="prediction",level=conf_level) # for plot confidence line
y2<-predict(fit,newdata=data.frame(x4fit_o2=x4fit_o2),interval="prediction",level=conf_level) # all genes for choosing genes
hvg2<-rownames(y2)[y4fit_o2[rownames(y2)]>y2[,3]]
print(paste('hvg2:',length(hvg2<- hvg2[x4fit_o2[hvg2]>xtem & y4fit_o2[hvg2]<(max(y4fit_o2)-0.01)])))

if(is_plot){
pdf(paste(path1,part,'_method2_hvg.pdf',sep=''),height=12,width=12)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3)
plot(x4fit_o2,y4fit_o2,frame=F,cex=0.2,xlab='log2 mean of count',ylab='log2 CV2',col='grey',xaxt='n',yaxt='n',main=title)
#abline(fit,col='red')
#lines(log2(c(0.1,0.1)),c(-6,8),col='green3',lty=2)
#lines(log2(c(100,100)),c(-6,8),col='green3',lty=2)
lines(x1,y1[,3],lty=1,col="red",lwd=5)
points(x4fit_o2[hvg2],y4fit_o2[hvg2],col='black',cex=0.3)
if(add_pt[1]!='') points(x4fit_o2[add_pt],y4fit_o2[add_pt],col='blue',cex=0.3)
if(add_pt2[1]!='') points(x4fit_o2[add_pt2],y4fit_o2[add_pt2],col='red',cex=0.3)
lines(c(xtem,xtem),c(-5,15),col='green3',lty=2,lwd=3)
axis(2,tck=-0.02,lwd=3)
axis(1,tck=-0.02,lwd=3)
dev.off()
}
return(union(hvg1,hvg2))
}

# Find hvg by sample: used in iterative clustering
find_hvg_by_emb<-function( all_cell, cut_num=1000, out_path='',xmin=log2(0.1),adj_int=0,ex_emb='',all_mx=norm_mx,rmg=rm_gene){
all_hvg<-list('')
for( emb in emb_id){
  if(emb %in% ex_emb) next
  cells<- all_cell[ sapply(all_cell, function(x){ strsplit(x, split='_')[[1]][1] })==emb ]
  if(length(cells)<cut_num) next
  emb_mx<- all_mx[, cells]
  if( length(cells)> 7000 ) int<-50 else if(length(cells)> 4000) int<-20 else if(length(cells)>2000) int<-10 else int<-0
  if( sum( emb %in% names(adj_int)) >0) int<- adj_int[emb]
  print( paste('int:',int))
  hvg<-setdiff(get_hvg(mx_norm=emb_mx,path1=out_path,part=emb,cut_agg=0,bou_int=int,conf_level=0.999,xtem=xmin), rmg) 
  print( paste(emb,':', length(hvg)))
  all_hvg<- c( all_hvg, list(hvg))
}
all_hvg<-all_hvg[-1]
all_hvg<-unique(unlist(all_hvg))
}

# Find hvg by embryo: 1) by sample; 2) by early and late embryos; Used in iterative clustering.
find_hvg_by_emb2<-function( all_cell, cut_num=100, out_path='',xmin=log2(0.1),adj_int=0,ex_emb='',all_mx=norm_mx,rmg=rm_gene, early_num=1, late_num=1){
early<- emb_id[ gsub("[A-z]", "", sapply(emb_id,function(y){ strsplit(y,split='_')[[1]][1] }) ) %in% c('7','0')]
late<- emb_id[ ! gsub("[A-z]", "", sapply(emb_id,function(y){ strsplit(y,split='_')[[1]][1] }) ) %in% c('7','0')]
sam_id<-sapply(all_cell, function(x){ strsplit(x, split='_')[[1]][1] })
sam_id[sam_id=='h9b']<-'h9a'
pass<- table(sam_id)
pass<- names(pass)[pass>=cut_num]
all_early<- all_cell[sam_id %in% early]
if(sum(pass %in% early)>=early_num & sum(pass %in% late)>=late_num ){
all_hvg<-list('')
for( emb in pass){
  if(emb %in% ex_emb) next
  cells<- all_cell[ sam_id==emb ]
  if(length(cells)<cut_num) next
  emb_mx<- all_mx[, cells]
  if( length(cells)> 7000 ) int<-50 else if(length(cells)> 4000) int<-20 else if(length(cells)>2000) int<-10 else int<-0
  if( sum( emb %in% names(adj_int)) >0) int<- adj_int[emb]
  print( paste('int:',int))
  hvg<-setdiff(get_hvg(mx_norm=emb_mx,path1=out_path,part=emb,cut_agg=0,bou_int=int,conf_level=0.999,xtem=xmin), rmg) 
  print( paste(emb,':', length(hvg)))
  all_hvg<- c( all_hvg, list(hvg))
}
all_hvg<-all_hvg[-1]
all_hvg<-unique(unlist(all_hvg))
way<-paste(1,paste(pass,collapse=','),sep=',')
}else if( length(all_early)>=80 & sum(pass %in% late)>=1 ){
  all_hvg<-list('')
  emb_mx<-all_mx[, all_early]
  if( length(all_early)> 7000 ) int<-50 else if(length(all_early)> 4000) int<-20 else if(length(all_early)>2000) int<-10 else int<-0
  hvg<-setdiff(get_hvg(mx_norm=emb_mx,path1=out_path,part='early',cut_agg=0,bou_int=int,conf_level=0.999,xtem=xmin), rmg) 
  all_hvg<- c( all_hvg, list(hvg))
  for( emb in pass[pass %in% late]){
  if(emb %in% ex_emb) next
  cells<- all_cell[sam_id==emb ] 
  if(length(cells)<cut_num) next
  emb_mx<- all_mx[, cells]
  if( length(cells)> 7000 ) int<-50 else if(length(cells)> 4000) int<-20 else if(length(cells)>2000) int<-10 else int<-0
  if( sum( emb %in% names(adj_int)) >0) int<- adj_int[emb]
  print( paste('int:',int))
  hvg<-setdiff(get_hvg(mx_norm=emb_mx,path1=out_path,part=emb,cut_agg=0,bou_int=int,conf_level=0.999,xtem=xmin), rmg) 
  print( paste(emb,':', length(hvg)))
  all_hvg<- c( all_hvg, list(hvg))
  }
  all_hvg<-all_hvg[-1]
  all_hvg<-unique(unlist(all_hvg))
  way<-paste(2,'early',paste(pass,collapse=','),sep=',')
}else{
  if( length(all_cell)> 7000 ) int<-50 else if(length(all_cell)> 4000) int<-20 else if(length(all_cell)>2000) int<-10 else int<-0
  all_hvg<-setdiff(get_hvg(mx_norm=all_mx[,all_cell],path1=out_path,part='all',cut_agg=0,bou_int=int,conf_level=0.999,xtem=xmin), rmg) 
  way<-paste(3,'all',sep=',')
}
return(list(all_hvg,way))
}


# 2-round of clustering
r2_wrap<- function( cell, mx, remove_gene, path='', hvg_adj=0, phe_r=c(0.5,0.8),which_phe=1, give_r1_hvg=0, cut_num=50, return_r1=F, is_r2_hvg_by_emb=T,r1_hvg_cut_num=300, r2_hvg_early=1, r2_hvg_late=1, rm_hox=F, is_merge=F, r1_pc_dim=1:20){
  if( rm_hox ) remove_gene<- c(remove_gene,names(grep('^HOX',ge2an,value=T)) )
  report<-paste(path,'report.txt',sep='')
  r1_brg<-get_broad_by_sam( g=rownames(mx) , x=cell, mx=mx)
  if(give_r1_hvg[1]==0){
    r1_hvg<- find_hvg_by_emb( cell, out_path=paste(path,'r1_',sep=''), cut_num=r1_hvg_cut_num,xmin=-2.5, adj_int=hvg_adj,rmg=remove_gene,all_mx=mx ) 
    r1_hvg<-setdiff( r1_hvg, r1_brg)
    write(r1_hvg, file=paste(path,'r1_hvg.txt',sep=''))
  }else r1_hvg<-give_r1_hvg
  if( rm_hox) r1_hvg<-setdiff( r1_hvg, names(grep('^HOX',ge2an,value=T)))
  r1_umap<-cal_umap_and_clu_seu3(x=cell,hvg=r1_hvg,npc=50,dims=r1_pc_dim,mx=mx,nnei=20,md=.1, is_umap=T, reso=phe_r)  
  plot_emb_source(umap=r1_umap[[1]],path=path,file_name='r1',plot_cex=.5)
  r1_clu<-r1_umap[[3]]
  for( i in 1:length(phe_r)) map_phe_clu(tsne=r1_umap[[1]], col=col6[-1], path1=paste(path,'r1_r',phe_r[i], sep=''),file_name='',ind=r1_clu[,i],plot_cex=.5)
  r1_clu_ind<- r1_clu[,which_phe]
  r1_mer_clu<- tapply( rownames(r1_umap[[1]]), r1_clu_ind, function(x){x})
  r1_res<-list(r1_umap,r1_mer_clu)
  if(return_r1) return(r1_res)
  write( paste('r1', length(cell),length(table(r1_clu_ind)),length(r1_mer_clu),sep='\t'), file=report)
  # round2
  r2_res<-list('')
  for( i in 1:length(r1_mer_clu)){
    if( length(r1_mer_clu[[i]])<=cut_num ){
      not_div<-list(r1_mer_clu[[i]])
      names(not_div[[1]])<-'1'
      r2_res[i]<-list(list(0,not_div))
      write( paste( paste('clu',i,sep=''), length(r1_mer_clu[[i]]), 1, 1,sep='\t'), file=report,append=T)
      next
    }
    r2_brg<-get_broad_by_sam( g=rownames(mx) , x=r1_mer_clu[[i]])
    if(is_r2_hvg_by_emb) r2_hvg<- find_hvg_by_emb2( r1_mer_clu[[i]], out_path=paste(path,'r2/clu',i,sep=''), cut_num=100,xmin=-2.5, adj_int=0,rmg=remove_gene,all_mx=mx, early_num=r2_hvg_early, late_num=r2_hvg_late ) 
    else r2_hvg<-setdiff(get_hvg(mx_norm=mx[,r1_mer_clu[[i]]], path1=paste(path,'r2/clu',i,sep=''), part='' ,cut_agg=0,bou_int=0,conf_level=0.999,xtem=-2.5), remove_gene) 
    r2_hvg_way<-r2_hvg[[2]]
    r2_hvg<-setdiff(r2_hvg[[1]], c(r1_brg,r2_brg))
    print( paste('clu',i,'HVG:',length(r2_hvg)))
    r2_umap<-cal_umap_and_clu_seu3(x=r1_mer_clu[[i]],hvg=r2_hvg,npc=50,dims=1:20,mx=mx,nnei=20,md=.1, is_umap=T, reso=phe_r[which_phe])  
    plot_emb_source(umap=r2_umap[[1]],path=path,file_name=paste('r2/clu',i,sep=''),plot_cex=1)
    r2_clu<-r2_umap[[3]]
    map_phe_clu(tsne=r2_umap[[1]], col=col6[-1], path1=paste(path,'r2/clu',i,'_r',phe_r, sep=''),file_name='',ind=r2_clu,plot_cex=1)
    r2_clu_ind<- r2_clu
    r2_mer_clu<- tapply( rownames(r2_umap[[1]]), r2_clu_ind, function(x){x})
    r2_res[i]<-list(list(r2_umap,r2_mer_clu))
    print(paste('FINISH clu',i,':',length(r2_mer_clu)))
    write( paste( paste('clu',i,sep=''), length(r1_mer_clu[[i]]), length(table(r2_clu_ind)),length(r2_mer_clu),r2_hvg_way,sep='\t'), file=report,append=T)
  }
  return(c(r1_res,r2_res))
}

# sub-function of 'r2_wrap': remove broad expressed genes during clustering 
get_broad_by_sam<-function(x,g,mx=norm_mx,num=100, perc=0.7){
  x_sam<-sapply(x,function(y){ strsplit(y,split='_')[[1]][1] })
  check_sam<- names(table(x_sam))[ table(x_sam)>=num ]
  if(length(check_sam)==0) return(NA)
  res<- unique(unlist(lapply( check_sam, function(y){
    broad<-g[ apply(mx[g, x[x_sam==y] ]>0, 1, sum)/sum(x_sam==y)>= perc ] 
    return(broad)
  })))
  return(res)
}

# sub-function of 'r2_wrap': calulate Umap and do clustering
cal_umap_and_clu_seu3<-function(x,hvg,npc=50,dims=1:30,mx=norm_mx,nnei=20,md=.1, is_umap=T, kr=30, reso=c(0.7,.8,.9), cycle=2){
seu <- CreateSeuratObject(mx[,x] , min.cells = 0, min.feature = 0)
seu <- NormalizeData(object = seu)
seu <- ScaleData(seu)
seu <- RunPCA(object = seu, features = hvg, do.print = F)
if(length(dims)>length(hvg)){
   dims<-1:(length(hvg)-1)
   print(paste('MY WARNING: HVGs size:',length(hvg)))
   print(paste('change dims:'))
   print(dims)
}
if(is_umap){
  seu <- RunUMAP(object = seu, dims = dims)
  umap<-seu@reductions$umap@cell.embeddings
}else{
  seu <- RunTSNE(object = seu, dims = dims)
  umap<-seu@reductions$tsne@cell.embeddings
}
pca<-seu@reductions$pca
seu=FindNeighbors(seu, dims = dims)
res <- rep(0, length(x))
if(cycle==1) ncyc<-length(kr) else ncyc<-length(reso)
for(k in 1:ncyc){ 
  if(cycle==1){
    seu=FindClusters(seu, resolution =reso, dims.use =dims, k.param = kr[k] )
  }else{
    seu=FindClusters(seu, resolution =reso[k], dims.use =dims, k.param = kr )
  }
  one_res<- as.numeric(seu@active.ident)
  names(one_res)<- names(seu@active.ident)
  res<- cbind( res, one_res)
}
res<-res[,-1]
all_res<-list( umap, pca, res)
return(all_res)
}

# sub-function of 'r2_wrap': plot Umap colored by embryo
plot_emb_source<-function(umap,path,file_name,plot_cex=.5,legend_left=T, is_umi=F, is_body=T, legend_pos='', first=''){
rand_ind<-sample.int(dim(umap)[1])
part_col<-col6[c(2,8,8,6,4,10,2,8,6,4,10,2,8,6,4,2,8,4,6)]
names(part_col)<-emb_id
pdf(paste(path,file_name,'_umap_emb.pdf',sep=''),width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
plot(umap[rand_ind,1],umap[rand_ind,2],pch=16,xlab='umap-1',ylab='umap-2',frame=F,col=part_col[ sapply(rownames(umap),function(x){strsplit(x,split='_')[[1]][1]}) ][rand_ind],cex=plot_cex)
if( first[1]!='') points(umap[first,1],umap[first,2],pch=16,col=part_col[ sapply(rownames(umap[first,]),function(x){strsplit(x,split='_')[[1]][1]}) ],cex=plot_cex)
dev.off()
if(is_body){
  part_col<-col6[c(3,3,3,3,3,7,11,11,5,5,10,2,2,2,2,6,6,6,5)]
  names(part_col)<-emb_id
  pdf(paste(path,file_name,'_umap_body.pdf',sep=''),width=12,height=12)
  par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
  plot(umap[rand_ind,1],umap[rand_ind,2],pch=16,xlab='umap-1',ylab='umap-2',frame=F,col=part_col[ sapply(rownames(umap),function(x){strsplit(x,split='_')[[1]][1]}) ][rand_ind],cex=plot_cex)
  if(legend_left) legend('topleft',legend=c('head','trunk','limb','vis','HeTr','TaVi','LiVi'), col=part_col[c(1,12,7,16,6,11,10)], bty='n',pch=16)
  dev.off()
}
#total UMI
if(is_umi){ 
data<- all_raw_total[rownames(umap)]
data[data<5000]<-5000
data[data>20000]<-20000
pdf(paste(path,file_name,'_umap_umi.pdf',sep=''),width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
tmp<-scatter_fill2(umap[,1],umap[,2],data[rownames(umap)],nlevels=999,title='total raw UMIs',xlab='UMAP-1',ylab='UMAP-2',pch=16,plot_cex1=plot_cex,plot_cex2=plot_cex,frame=F,is_0f=T, col=-1 )
dev.off()
}
}


# sub-function of 'r2_wrap': map all cluster on Umap
map_phe_clu<-function(tsne,col,is_plot=T,path1,mx,add_pt='',file_name='ite_clu',is_text=T,ind,plot_cex=.3, is_outlier=F, add_col='red'){
if(!is_outlier) center<-cbind(tapply(tsne[,1], ind[rownames(tsne)] ,median),tapply(tsne[,2],ind[rownames(tsne)],median))[as.character(1:length(unique(ind))),] else center<-cbind(tapply(tsne[,1], ind[rownames(tsne)] ,median),tapply(tsne[,2],ind[rownames(tsne)],median))[as.character(0:(length(unique(ind))-1)),]
clu_num<-length(unique(ind))
if(sum(ind==0)>0){
  center<-center[-1,]
  clu_num<-clu_num-1
}
if( length(unique(ind))==1 ) center<-rbind(center,center)
if(is_plot){
pdf(paste(path1,file_name,'.pdf',sep=''),width=12,height=12)
par(cex=2,las=1,mar=c(4,4,1,2),lwd=3)
if(!is_outlier) plot(tsne[,1],tsne[,2],pch=20,xlab='Umap-1',ylab='Umap-2',frame=F,col=col[as.numeric(ind[rownames(tsne)])],cex=plot_cex,main=paste('clu: ',clu_num,sep='')) else plot(tsne[,1],tsne[,2],pch=20,xlab='Umap-1',ylab='Umap-2',frame=F,col=col[as.numeric(ind[rownames(tsne)])+1],cex=plot_cex,main=paste('clu: ',clu_num,sep='')) 
if(add_pt[1]!='') points(tsne[add_pt,1],tsne[add_pt,2],cex=plot_cex, col=add_col,pch=20)
else if(is_text){
  if(sum(ind==0)>0) text(center[,1],center[,2],label=1:(length(unique(ind))-1),xpd=NA )
  else text(center[,1],center[,2],label=1:(length(unique(ind))),xpd=NA )
}
dev.off()
}
}

# sub-function of 'r2_wrap': reform result of r2_wrap into linear
linear_r2_clu<-function( clu, name ){
  res<-list('')
  for(j in 3:length(clu)){
    r2_cell<-do.call('c',clu[[j]][2])
    if( r2_cell[[1]][1]==0 ) r2_cell<- clu[[2]][j-2]
    res<-c( res, r2_cell )
    clu_name<-paste(name,j-2,1:length(clu[[j]][[2]]),sep='-')
    names(res)[ (length(res)-length(clu_name)+1):length(res) ]<-clu_name
  }
  res<-res[-1]
  ind<- rep( 1:length(res), time=sapply(res,length))
  names(ind)<-unlist(res)
  return(list(res,ind))
}

# Plot gene expression on Umap by batch
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
  pdf(paste(path,file_name,'.pdf',sep=''),height=12,width=9)
  par(las=1,mar=c(5,1,1,4),lwd=2, mfrow=c(4,3))
  if(length(genes)==1) res<-plot_mk_in_cluster4(x=plot_gene,plot_max=data_max[1],gene_name=title_lab,title=title_lab,tsne_res=tsne,count_mx=mx[,rownames(tsne)],nl=999,is_0f=zerof,plot_cex1=plot_cex, plot_cex2=plot_cex, xlm=xx, ylm=yy)
  else res<-mapply(plot_mk_in_cluster4,x=plot_gene,plot_max=data_max,gene_name=title_lab,title=title_lab,MoreArgs=list(tsne_res=tsne,count_mx=mx[,rownames(tsne)],nl=999,is_0f=zerof,plot_cex1=plot_cex, plot_cex2=plot_cex, xlm=xx, ylm=yy))   
  dev.off()
  if(is_ret) return(plot_gene)
}

# sub-function of 'plot_genes_on_tsne'
#source('scatter_fill2.r')
plot_mk_in_cluster4<-function(x,is_plot=T,noise=-1,gene_name='',folder='',tsne_res,count_mx,plot_cex1=0.5,plot_cex2=0.5,title='',plot_max=0,is_0f=F,nl,is_axis='n',main_cex=1,ret_data=F,xlm=c(0,0),ylm=c(0,0)){
    if(x==''){
       plot(0,0,col='white',xaxt='n',yaxt='n',frame=F,main= '',xlab='',ylab='')
       return(paste(x,gene_name,'blank'))
    }
    if( sum(x %in% rownames(count_mx))==0){
       plot(0,0,col='white',xlab='',ylab='',frame=F,xaxt='n',yaxt='n',main=paste(ge2an[x], 'no detection'))
       return(paste(x,gene_name,'no detection'))
    }else if( sum(count_mx[x,]>0)==0 ){
       plot(0,0,col='white',xlab='',ylab='',frame=F,xaxt='n',yaxt='n',main=paste(ge2an[x], 'no detection'))
       return(paste(x,gene_name,'no detection'))
    }
    if( sum( x %in% rownames(count_mx)) == 0){
       return(paste(x,gene_name,'no detection'))
    }
    if(sum(as.numeric(count_mx[x,])>0)==0) return(paste(x,gene_name,'all 0'))
    plot_pch<-16
    data<-as.numeric(count_mx[x,])
    if(plot_max==0){ 
      plot_max<-max(data)
      if(plot_max>10){
        m1<- quantile(data[data>0], probs=0.90)
  m2<- max(data)*.45
  plot_max<- max(m1,m2)
      }            
    }
    data[data> plot_max]<-plot_max
    if(nl==0) nl=plot_max+1
    if(nl<99) nl=99
    plot_min<-0
    data[data<plot_min]<-plot_min 
    if(xlm[1]==0&xlm[2]==0){
      xlm<-range(tsne_res[,1])
      ylm<-range(tsne_res[,2])
    }
    res<-scatter_fill2(tsne_res[,1],tsne_res[,2],data,zlim=c(plot_min,plot_max),nlevels=999,title=title,xlab='',ylab='',pch=plot_pch,plot_cex1=plot_cex1,plot_cex2=plot_cex2,frame=F,is_0f=is_0f, col=-1,is_axis=is_axis,main_cex=main_cex,xlim=xlm, ylim=ylm )
    if(!ret_data) return(c(plot_min, plot_max)) else return(res)
}


# a wrap for Seurat function 'FindAllMarkers'
findMarkerbySeurat <- function(mx,cluster.info,clusterInfo_format='list',cells.1=NULL,top_n=10) {
  library(dplyr);library(Seurat)
  if(clusterInfo_format=='list') {ident=get_clu_ind(cluster.info, cell=unlist(cluster.info))}
  if(clusterInfo_format=='vector') {ident=cluster.info}
  cell=names(ident)
  seu=CreateSeuratObject(counts = mx[,cell],min.cells = 0,min.features = 0,project = "seu")
  seu=NormalizeData(seu,normalization.method = "LogNormalize",scale.factor = 1e4)
  seu=ScaleData(seu)
  seu@active.ident=as.factor(ident)
  if(is.null(cells.1)) {seu.marker=FindAllMarkers(seu,only.pos = T); seu.marker$gene=ge2an[seu.marker$gene]; y=seu.marker%>%group_by(cluster)%>%top_n(top_n,wt = avg_logFC) }
  if(!is.null(cells.1)) {seu.marker=FindMarkers(seu,cells.1 = ident.1,only.pos = T); rownames(seu.marker)=ge2an[rownames(seu.marker)]; y=seu.marker[1:min(top_n,nrows(seu.marker)),]}
  return(list(seu.marker,y[,c('cluster','gene')]))
}

# get cluster index based on list
get_clu_ind <- function( clu, cell){
  clu_ind=rep( names(clu), time=sapply(clu,length))
  names(clu_ind)=unlist(clu)
  clu_ind=clu_ind[cell]
  return(clu_ind)
}

do_pheatmap <- function(genes,cells,mx=fl_mx,clusterRow=T) {
  x=mx[genes,cells]; rownames(x)=ge2an[rownames(x)]
  x['TSPAN13',]= sapply(x['TSPAN13',],function(y) ifelse(y<2,0,y))
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

get_body_part <- function(cellname=cells){
  body=sapply(cellname,function(x){
    if(substr(x,start = 1,stop = 3) %in% c("h0_","h5_","h6_","h9a","h9b")) return("head")
    else if(substr(x,start = 1,stop = 3) %in% c("t0_","t5_","t9_","t6_")) return("trunk")
    else if(substr(x,start = 1,stop = 3) %in% c("l0_","l9_")) return("limb")
    else if(substr(x,start = 1,stop = 3) %in% c("v0_","v9_","v6_")) return("vi")
    else if(substr(x,start = 1,stop = 3) %in% c("ht7")) return("ht")
    else if(substr(x,start = 1,stop = 3) %in% c("tv7")) return("tv")
    else if(substr(x,start = 1,stop = 3) %in% c("lv5","lv6","vl5")) return("lv")
  })
  return(unlist(body))
}

mono_pseudo <- function(mx, cell, gene, is_order=F, ncen=0){
  mono<- newCellDataSet( mx[,cell],expressionFamily=negbinomial.size())
  mono <- estimateSizeFactors(mono)
  mono <- estimateDispersions(mono,remove_outliers =F)
  mono <- setOrderingFilter(mono, gene)
  set.seed(4837)
  if(ncen==0) mono<-reduceDimension(mono, max_components = 2, reduction_method = 'DDRTree')
  else mono<-reduceDimension(mono, max_components = 2, reduction_method = 'DDRTree', ncenter=ncen)
  if(is_order) mono <- orderCells(mono)
  return(list(t(mono@reducedDimS),mono@phenoData@data))
}


get_smooth <- function(mx,bin=ncol(mx)/20,onlyLeft=F){
  if(!onlyLeft) {mx.s=sapply(1:ncol(mx),function(i){rowMeans(as.matrix(mx[,max(1,(i-bin)):min(ncol(mx),(i+bin))]))})}
  if(onlyLeft) {mx.s=sapply(1:ncol(mx),function(i){rowMeans(as.matrix(mx[,max(1,(i-bin)):i]))})}
  colnames(mx.s)=colnames(mx)
  return(mx.s)
}

get_scale<- function(mx,by_row=T) {
  if(by_row) {out=t( apply(mx,1,function(x){ (x-min(x)) / (max(x)-min(x)) }) )}
  if(!by_row) {out=apply(mx,2,function(x){ (x-min(x)) / (max(x)-min(x)) })}
  return(out)
}

avgExp <- function(mx,clu.list,gene=rownames(mx)) {
  if(class(clu.list) %in% c("list","array")) {
    clu.list=clu.list
  } else {
    clu.list=tapply(names(clu.list),clu.list, function(x) x ) 
  }
  tp_mn=sapply( clu.list, function(x) { apply( mx[gene, x,drop=F], 1, mean ) } )
  return(tp_mn)
}
pctExp <- function(mx,clu.list,gene=rownames(mx)) {
  if(class(clu.list) %in% c("list","array")) {
    clu.list=clu.list
  } else {
    clu.list=tapply(names(clu.list),clu.list, function(x) x ) 
  }
  tp_fr=sapply( clu.list, function(x) { apply( mx[gene, x,drop=F], 1, function(y) { sum(y>0)/length(y) } ) } )
  return(tp_fr)
}

plot_pair_bubble_sideByside <- function( pct.mx1, pct.mx2, avg.mx1, avg.mx2, scale1=3, scale2=3, color1="red", color2="blue", srt= 90, gap=0.25,sep.line=F,name1="human",name2="mouse",ylab_xadj=0,labelDiffSites=F){
  if(ncol(pct.mx1)!=ncol(avg.mx1)) { stop("Different colnums of pct.mx1 and avg.mx1 !")}
  if(ncol(pct.mx2)!=ncol(avg.mx2)) { stop("Different colnums of pct.mx2 and avg.mx2 !")}
  if(ncol(pct.mx1)!=ncol(pct.mx2)) { stop("Different colnums of pct.mx1 and pct.mx2 !")}
  if(labelDiffSites){
    diff.mx=get_diff_sites(pct.m=pct.mx2,pct.h=pct.mx1,avg.m=avg.mx2,avg.h=avg.mx1,gene=rownames(pct.mx1))
    diff.mx=as.data.frame(diff.mx)
  }
  avg.mx1=log2(1+avg.mx1); avg.mx2=log2(1+avg.mx2)
  x_num=nrow(pct.mx1); y_num=ncol(pct.mx1)
  x_lab=rownames(pct.mx1); y_lab=colnames(pct.mx1)
  
  par(las=1,mar=c(6,6,1,1),lwd=2,pch=16)
  plot(NA,frame=F,xlim=c(1.8,2*x_num+4.6),ylim=c(1+gap,(1+gap)*y_num),xaxt='n',yaxt='n',xlab='',ylab='',main='') 
  box()
  points(x=rep(2*x_num+2,7),y=seq(from=((1+gap)*y_num-1),to=((1+gap)*y_num-8.5),by = -1.25),cex=scale1*c(0.75,0.75,0,0.25,0.5,0.75,1),col=c(color1,color2,rep("black",5)))
  text( x=rep(2*x_num+4.4,7),y=seq(from = ((1+gap)*y_num-1),to=((1+gap)*y_num-8.5),by = -1.25),labels = c(name1,name2,'','25%','50%','75%','100%'))
  
  for(i in 1:y_num) {
    avg.mx1.grade=as.numeric(cut(x=c(min(avg.mx1),avg.mx1[,i,drop=T],max(avg.mx1)), breaks=20)[2:(x_num+1)])
    cols1=mapply(function(color1, value1) {return(colorRampPalette(c('grey',color1))(20)[value1])},color1=color1, value1=avg.mx1.grade)
    avg.mx2.grade=as.numeric(cut(x=c(min(avg.mx2),avg.mx2[,i,drop=T],max(avg.mx2)), breaks=20)[2:(x_num+1)])
    cols2=mapply(function(color2, value2) {return(colorRampPalette(c('grey',color2))(20)[value2])},color2=color2, value2=avg.mx2.grade)
    points(x=2*(1:x_num)-1,y=(1+gap)*rep(y_num+1-i,x_num), cex=scale1*pct.mx1[,i,drop=T],xpd=NA,col=cols1)
    points(x=2*(1:x_num),y=(1+gap)*rep(y_num+1-i,x_num),cex=scale2*pct.mx2[,i,drop=T],xpd=NA,col=cols2)
  }
  
  text( 2*(1:x_num), (par('usr')[3]-par('usr')[4])/500 , label=x_lab,xpd=NA, pos=2,srt=srt) # xlab
  text( (par('usr')[2]-par('usr')[1])/1000+ylab_xadj,  (1+gap)*(y_num:1),label=y_lab,xpd=NA, pos=2) # ylab
  
  abline(v=2*x_num+1)
  
  if(sep.line){abline(v=2*(1:(x_num-1))+0.5, lwd=0.5,lty=2)}
  
  if(labelDiffSites) {
    for(i in 1:y_num) {
      y_low=(2+gap)*(y_num+1-i)-1.5; y_high=(2+gap)*(y_num+1-i)+0.5
      site=as.numeric(which(diff.mx[,i]== -1))
      for(j in 1:length(site)){
        polygon(x=c(site[j]-0.5,site[j]+0.5,site[j]+0.5,site[j]-0.5),
                y=c(y_low,y_low,y_high,y_high),border ='blue')
      }
      
      site=as.numeric(which(diff.mx[,i]== 1))
      for(j in 1:length(site)){
        polygon(x=c(site[j]-0.5,site[j]+0.5,site[j]+0.5,site[j]-0.5),
                y=c(y_low,y_low,y_high,y_high),border ='red')
      }
      
    }
  }
}

get_col <- function(val,low='grey99',high='red3',max=NULL) {
  grade=as.numeric(cut(c(val,max), breaks=100))
  col=mapply(function(x){return(colorRampPalette(c(low,high))(100)[x])},x=grade)
  return(col[1:length(val)])
}

plot_bubble <- function( pct.mx, avg.mx, scale=5, color="red",xlab.srt=90,fontsize_row=1,feaEQmx=T,xlab_xadj=-0.05,addLine=F) {
  #' determine xy range, with genes in x axis, types in y axis '
  if(ncol(pct.mx)!=ncol(avg.mx)) { stop("Different colnums of pct.mx and avg.mx ==|") }
  if(nrow(pct.mx)!=nrow(avg.mx)) { stop("Different rownums of pct.mx and avg.mx ==|") }
  avg.mx=log2(1+avg.mx)
  x_num=nrow(pct.mx); y_num=ncol(pct.mx)
  if(feaEQmx) {x_lab=rownames(pct.mx); y_lab=colnames(pct.mx)}
  else{x_lab=ge2an[rownames(pct.mx)];y_lab=colnames(pct.mx)}
  
  
  par(las=1, mar=c(6,10,0,1), lwd=2, pch=16)
  plot(NA, frame=F, xlim=c(1,x_num), ylim=c(1,y_num), xaxt='n', yaxt='n', xlab='', ylab='', main='') 
  #' plot dots by rows/genes '
  if(addLine) {abline(v=1:x_num,h=1:y_num,lty=2,col='grey',lwd=0.5)}
  for(i in 1:y_num) {
    avg.exp.grade=as.numeric(cut(x=c(min(avg.mx),avg.mx[,i,drop=T],max(avg.mx)), breaks=20)[2:(x_num+1)])
    cols=mapply(function(color, value) {return(colorRampPalette(c('grey',color))(20)[value])},color=color, value=avg.exp.grade)
    points(x=1:x_num, y=rep(y_num+1-i,x_num), cex=scale*pct.mx[,i,drop=T], xpd=NA, col=cols)
  }
  text( (1:x_num)+xlab_xadj, (par('usr')[3]-par('usr')[4])/100 , label=x_lab, xpd=NA, pos=2, srt=xlab.srt,cex=fontsize_row)
  text( (par('usr')[2]-par('usr')[1])/200, y_num:1, label=y_lab, xpd=NA, pos=2)
 
}
