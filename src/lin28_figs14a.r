#' Punctuated embryonic stages in vertebrates.
#' TF waves of punctuated embryonic stages in zebrafish and frog

# Zebrafish
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
# Load functions
source('../../lin28_function.r')
# Define stage-specific TFs by 1) this stage > all other stages; 2) high correlation with a stardard vector of stage-specific gene expression
tmp<-unique(read.table('~/baolab/xuy/ciona/cluster2/product/code/type_all_col3.txt',sep='\t',as.is=T,comment.char='')[,2])
col6<-c('lightgrey','red',tmp[3],'purple',tmp[2],tmp[5],'blue','green4','orchid','turquoise','sienna',tmp[9],'yellow2','hotpink','navy','steelblue','skyblue','pink','black',tmp[4],rainbow(7))
load(file='../../White/code/dr_all_exp.rdata')
go_tf<- read.table('../../White/ensembl_GO0003700.txt',sep='\t',as.is=T,header=T) # zebrafish TFs
rownames(go_tf)<-go_tf[,1]
tf<- intersect( rownames(dr_all_exp), go_tf[,1])
tf_exp<- dr_all_exp[tf,]
st_anno<- c(paste( c(0,0.75,2.25,3,4.3,5.25,6,8,10.3,16,19,24,30,36),'hpf',sep=''),'2d','3d','4d','5d')
# 1) this stage > all other stages
tf_cmp<- cbind( cal_wilp2( mx1=tf_exp[,6:8], mx2=tf_exp[,9:15],cut_mean=10,small_val=1), # s2-s3
cal_wilp2( mx1=tf_exp[,9:15], mx2=tf_exp[,16:18],cut_mean=10,small_val=1), # s3-s4
cal_wilp2( mx1=tf_exp[,6:8], mx2=tf_exp[,16:18],cut_mean=10,small_val=1), # s2-s4
cal_wilp2( mx1=tf_exp[,1:5], mx2=tf_exp[,6:8],cut_mean=10,small_val=1), # s1-s2
cal_wilp2( mx1=tf_exp[,1:5], mx2=tf_exp[,9:15],cut_mean=10,small_val=1), # s1-s3
cal_wilp2( mx1=tf_exp[,1:5], mx2=tf_exp[,16:18],cut_mean=10,small_val=1) # s1-s4
) 
st_tf<- list( rownames(tf_cmp)[ tf_cmp[,4]> 1 & tf_cmp[,5]> 1 & tf_cmp[,6]> 1 ], # stage 1 TFs
rownames(tf_cmp)[ tf_cmp[,1]> 1 & tf_cmp[,3]> 1 & tf_cmp[,4]< -1 ],  # stage 2 TFs
rownames(tf_cmp)[ tf_cmp[,1]< -1 & tf_cmp[,2]> 1 & tf_cmp[,5]< -1 ],  # stage 3 TFs
rownames(tf_cmp)[ tf_cmp[,2]< -1 & tf_cmp[,3]< -1 & tf_cmp[,6]< -1 ] )  # stage 4 TFs
# 2) high correlation with a vector of stage-specific gene
(st_vec<-lapply(1:4, function(x){ vec<-rep(0,4); vec[x]<-1 ; res<-rep(vec,c(5,3,7,3)); return(res) }))
st_tf_cor<- mapply(function(x,y){
  res<-apply(dr_all_exp[x,], 1, function(z){ cor(z,y, method='pearson') })
  return(res)
}, x=st_tf[1:4], y=st_vec, SIMPLIFY=F )
tf_cor_cut<- 0.85
st_tf_corft <-sapply(st_tf_cor, function(x){  names(x)[x> tf_cor_cut ] })
for( i in 1:4){
  plot_dynamic( mx=dr_all_exp, gene=st_tf_corft[[i]], st_lab=st_anno, st_col=rep(col6[c(2:4,6)],c(5,3,7,3)), path='../result/s2_s3/', file_name=paste('corft_stage',i,'_tf',sep=''), name=paste('stage ',i,' TFs n=',length(st_tf_corft[[i]]),sep=''),lab_pos=.2)
}
dr_st_tf_corft<- st_tf_corft
# Plot each group of TFs as a wave
dr_st_tf_val <-lapply( dr_st_tf_corft, function(x) {
  val<-plot_dynamic( mx=dr_all_exp, gene=x, ret_val=T, file_name='tmp', path='../result')
  res<- cbind( apply(val,2,mean), apply(val,2,function(x){ sd(x) } ))
  return(res)
})
plot_col<- col6[c(2:4,6)]
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3,pch=16)
plot(0,0,xlim=c(1,ncol(dr_all_exp)),ylim=c(0,1),col='white',frame=F,xlab='',ylab='relative expression',main='',xaxt='n')
for(i in 1:4){
  polygon(x = c( 1:(dim(dr_st_tf_val[[i]])[1]) , rev(1:(dim(dr_st_tf_val[[i]])[1])) ), y = c( dr_st_tf_val[[i]][,1]+dr_st_tf_val[[i]][,2], rev(dr_st_tf_val[[i]][,1]-dr_st_tf_val[[i]][,2]) ),col = adjustcolor(plot_col[i], alpha.f = 0.5), border = NA)
}
text( 1:ncol(dr_all_exp), (-par('usr')[4]+par('usr')[3])*.2, label=st_anno, xpd=NA, cex=.8,srt=90, col=rep(col6[c(2:4,6)],c(5,3,7,3)) )
# add LIN28A expression
lin_exp<-dr_all_exp[c('ENSDARG00000004328','ENSDARG00000016999','ENSDARG00000052511'),]
lines(1:ncol(dr_all_exp), lin_exp[1,]/max(lin_exp[1,]), cex=1,xpd=NA,col='black',lwd=5)

# Frog
load('../ensembl100/code/map100.rdata')
# 1) this stage > all other stages 
tf_cmp<- cbind( cal_wilp2( mx1=tf_exp[,8:10], mx2=tf_exp[,11:17],cut_mean=10,small_val=1), # s2-s3
cal_wilp2( mx1=tf_exp[,11:17], mx2=tf_exp[,18:23],cut_mean=10,small_val=1), # s3-s4
cal_wilp2( mx1=tf_exp[,8:10], mx2=tf_exp[,18:23],cut_mean=10,small_val=1), # s2-s4
cal_wilp2( mx1=tf_exp[,1:7], mx2=tf_exp[,8:10],cut_mean=10,small_val=1), # s1-s2
cal_wilp2( mx1=tf_exp[,1:7], mx2=tf_exp[,11:17],cut_mean=10,small_val=1), # s1-s3
cal_wilp2( mx1=tf_exp[,1:7], mx2=tf_exp[,18:23],cut_mean=10,small_val=1) # s1-s4
) 
st_tf<- list( rownames(tf_cmp)[ tf_cmp[,4]> 1 & tf_cmp[,5]> 1 & tf_cmp[,6]> 1 ], # stage 1 TFs
rownames(tf_cmp)[ tf_cmp[,1]> 1 & tf_cmp[,3]> 1 & tf_cmp[,4]< -1 ],  # stage 2 TFs
rownames(tf_cmp)[ tf_cmp[,1]< -1 & tf_cmp[,2]> 1 & tf_cmp[,5]< -1 ],  # stage 3 TFs
rownames(tf_cmp)[ tf_cmp[,2]< -1 & tf_cmp[,3]< -1 & tf_cmp[,6]< -1 ] )  # stage 4 TFs
# common TFs in 1 stages
tf_st_mn<- t(apply(tf_exp, 1, function(x){ tapply(x, rep(1:4,c(7,3,7,6)),mean) }))
st_tf<- c( st_tf, list( setdiff(rownames(tf_st_mn)[ apply(tf_st_mn[,1:4]>=10,1,sum)==4],unlist(st_tf)) ))
# 2) high correlation with a stardard vector of stage-specific gene expression
(st_vec<-lapply(1:4, function(x){ vec<-rep(0,4); vec[x]<-1 ; res<-rep(vec,c(7,3,7,6)); return(res) }))
st_tf_cor<- mapply(function(x,y){
  res<-apply(xt_all_exp[x,], 1, function(z){ cor(z,y, method='pearson') })
  return(res)
}, x=st_tf[1:4], y=st_vec, SIMPLIFY=F )
st_tf_corft <-sapply(st_tf_cor, function(x){  names(x)[x>tf_cor_cut ] })
xt_st_tf_corft<- st_tf_corft
# Plot each group of TFs as a wave
xt_st_tf_val <-lapply( xt_st_tf_corft, function(x) {
  val<-plot_dynamic( mx=xt_all_exp, gene=x, ret_val=T, file_name='tmp', path='../result')
  res<- cbind( apply(val,2,mean), apply(val,2,function(x){ sd(x) } )) 
  return(res)
})
plot_col<- col6[c(2:4,6)]
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3,pch=16)
plot(0,0,xlim=c(1,ncol(xt_all_exp)),ylim=c(0,1),col='white',frame=F,xlab='',ylab='relative expression',main='',xaxt='n')
for(i in 1:4){
  polygon(x = c( 1:(dim(xt_st_tf_val[[i]])[1]) , rev(1:(dim(xt_st_tf_val[[i]])[1])) ), y = c( xt_st_tf_val[[i]][,1]+xt_st_tf_val[[i]][,2], rev(xt_st_tf_val[[i]][,1]-xt_st_tf_val[[i]][,2]) ),col = adjustcolor(plot_col[i], alpha.f = 0.5), border = NA)
}
text( 1:23, (-par('usr')[4]+par('usr')[3])*.25, label=colnames(xt_all_exp), xpd=NA, cex=.8,srt=90, col=rep(col6[c(2:4,6)],c(7,3,7,6)) )
la_exp<- xt_all_exp[ 'ENSXETG00000012324', ]
lines(1:length(la_exp), la_exp/max(la_exp), cex=1,xpd=NA,col='black',lwd=5)


# rmarkdown::render('lin28_figs14a.r')