#' LIN28A mRNA dynamics in vertebrate embryos

#' For human and mouse scRNA-seq, in each cell type, LIN28A expression is showed in normalized UMIs +/- SE
# Load human data
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
knitr::opts_chunk$set(warning = FALSE)
load('hs_data.rdata') # cell type, gene expression, annotation
tmp<-unique(read.table('~/baolab/xuy/ciona/cluster2/product/code/type_all_col3.txt',sep='\t',as.is=T,comment.char='')[,2])
col6<-c('lightgrey','red',tmp[3],'purple',tmp[2],tmp[5],'blue','green4','orchid','turquoise','sienna',tmp[9],'yellow2','hotpink','navy','steelblue','skyblue','pink','black',tmp[4],rainbow(7))
emb_col<-col6[c(2,3,3,6,4,10,2,3,6,4,10,2,3,6,4,2,3,4,6)]
emb_id<-c('h0','h9b','h9a','h5','h6','ht7','l0','l9','lv5','lv6','tv7','t0','t9','t5','t6','v0','v9','v6','vl5')
names(emb_col)<-emb_id
emb_col2<-emb_col
names(emb_col2)<- gsub("[A-z]", "", names(emb_col) )

# Calculate mean and SE in each cell type
get_id<-function(x){ sapply(x, function(y){names(ge2an)[ge2an==y][1]})}
get_emb<-function(x){
  emb<-gsub("[A-z]", "", sapply(x,function(y){ strsplit(y,split='_')[[1]][1] }) )
  return(emb)
}
la_ct_mn<-t(sapply( names(mg_clu), function(x){
  cl<-mg_clu[[x]]
  emb<-get_emb(cl)
  emb[emb %in% c(5,6)]<-'9'
  mn<-tapply( la_exp[get_id('LIN28A'),cl], emb, mean)
  res1<- rep(0,3)
  names(res1)<-c('7','0','9')
  res1[names(mn)]<-mn
  mn<-tapply( la_exp[get_id('LIN28B'),cl], emb, mean)
  res2<- rep(0,3)
  names(res2)<-c('7','0','9')
  res2[names(mn)]<-mn
  res<-c(res1,res2)
  return(res)
}))
la_ct_ci<-t(sapply( names(mg_clu), function(x,is_ci=F){ # 95%CI in each cell type by embryo or SE
  cl<-mg_clu[[x]]
  emb<-get_emb(cl)
  emb[emb %in% c(5,6)]<-'9'
  if(is_ci) sd<-tapply( la_exp[get_id('LIN28A'),cl], emb, function(y){ sd(y)*1.96/sqrt(length(y)) }) # 95% CI = 1.96*se, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1255808/, estimate the mean
  sd<-tapply( la_exp[get_id('LIN28A'),cl], emb, function(y){ sd(y)/sqrt(length(y)) })  # SE
  sd[is.na(sd)]<-0
  res1<- rep(0,3)
  names(res1)<-c('7','0','9')
  res1[names(sd)]<-sd
  sd<-tapply( la_exp[get_id('LIN28B'),cl], emb, function(y){ sd(y)/sqrt(length(y)) })  # SE
  sd[is.na(sd)]<-0
  res2<- rep(0,3)
  names(res2)<-c('7','0','9')
  res2[names(sd)]<-sd
  res<-c(res1,res2)
  return(res)
}))
#pdf('../result/lin28a_type_by_emb_ci.pdf',width=25) # use SE
#```{r, fig.width = 15}
par(cex=2,las=1,mar=c(2,4,0,0),lwd=3,pch=16)
j<-0
for(i in 1:dim(la_ct_mn)[1]){
  plot_mn<- la_ct_mn[i,1:3]
  plot_ci<- la_ct_ci[i,1:3]
  plot_ci<-plot_ci[plot_mn!= -.1]
  plot_mn<-plot_mn[plot_mn!= -.1]
  plot_mn[plot_mn>3]<-3
  if( T ) j<-j+1 else{
    j<-j+3 # separator
    lines( c(j-1.5,j-1.5), c(0,3), lty=2)
  }
  if(i==1) plot( rep(j,length(plot_mn)), plot_mn, col=emb_col2[ names(plot_mn) ] , xlab='', ylab='UMI mean', main='', frame=F, xlim=c(1,length(mg_clu)),ylim=c(0,3),xaxt='n',xpd=NA)
  else points( rep(j,length(plot_mn)), plot_mn, col=emb_col2[ names(plot_mn) ] ,xpd=NA)
  ind<- names(plot_mn) %in% c('9')
  if(sum(ind)>0) arrows( rep(j,length(plot_mn[ind])), plot_mn[ind], rep(j,length(plot_mn[ind])), plot_mn[ind]+plot_ci[ind], length = 0.1, angle = 90, col=emb_col2[names(plot_mn)[ind]],xpd=NA  )
  ind<- names(plot_mn) %in% c('0','7')
  if(sum(ind)>0) arrows( rep(j,length(plot_mn[ind])), plot_mn[ind], rep(j,length(plot_mn[ind])), plot_mn[ind]-plot_ci[ind], length = 0.1, angle = 90, col=emb_col2[names(plot_mn)[ind]],xpd=NA  )
}
for(j in col_sep) lines( c(j-1.5,j-1.5), c(0,3), lty=2)
#```
#dev.off()
# How strong and how many types are changing in human
la_jump<- t(apply(la_ct_mn, 1, function(x){
  res1<-ifelse( x[1]==0, (x[3]+.1)/(x[1]+.1), x[3]/x[1] )
  mn1<-x[1];mn2<-x[3]  
  res2<- ifelse( min(mn1,mn2)<0.1, log2((mn1+0.1)/(mn2+0.1)), log2(mn1/mn2))
  res<-c(res1,res2)
  return(res)
}))
# histogram on fold difference
plot_data<-1/la_jump[,1]
plot_data[plot_data>20]<-20
par(cex=2,las=1,mar=c(4,4,1,1),lwd=5,pch=16)
hist( plot_data, main='Human', ylab='Number of cell types', xlab='Fold difference', breaks=0:20)
text( 18,-8, label='>', xpd=NA)

# Load mouse data
load('../../../Shendure/code/shendure_lin28_data.rdata')
mm_col<-col6[c(10,2,3,11,7)]
names(mm_col)<-c( '9.5',   '10.5',   '11.5',  '12.5', '13.5') 
par(cex=2,las=1,mar=c(6,4,0,0),lwd=3,pch=16)
j<-0
for(i in 1:dim(la_mt_mn)[1]){
  plot_mn<- la_mt_mn[i,1:3] # only 3 stages
  plot_ci<- la_mt_ci[i,1:3]
  plot_ci<-plot_ci[plot_mn!= -.1]
  plot_mn<-plot_mn[plot_mn!= -.1]
  j<-i
  if(i==1) plot( rep(j,length(plot_mn)), plot_mn, col=mm_col[ names(plot_mn) ] , xlab='', ylab='UMI mean', main='', frame=F, xlim=c(1,38),ylim=c(0,max(la_mt_mn)),xaxt='n',xpd=NA)
  else points( rep(j,length(plot_mn)), plot_mn, col=mm_col[ names(plot_mn) ] ,xpd=NA)
  ind<- names(plot_mn) %in% c('11.5','12.5','13.5')
  if(sum(ind)>0) arrows( rep(j,length(plot_mn[ind])), plot_mn[ind], rep(j,length(plot_mn[ind])), plot_mn[ind]+plot_ci[ind], length = 0.1, angle = 90, col=mm_col[names(plot_mn)[ind]],xpd=NA  )
  ind<- names(plot_mn) %in% c('9.5','10.5')
  if(sum(ind)>0) arrows( rep(j,length(plot_mn[ind])), plot_mn[ind], rep(j,length(plot_mn[ind])), plot_mn[ind]-plot_ci[ind], length = 0.1, angle = 90, col=mm_col[names(plot_mn)[ind]],xpd=NA  )
}
text( 1:(dim(la_mt_mn)[1])+0.5, -0.005, label=rownames(la_mt_mn), cex=.5,xpd=NA,srt=90,pos=2)
for(j in c(16,20,36,37))   lines( c(j+.5,j+.5), c(0,1), lty=2)
# How strong and how many types are changing in mouse
la_jump<- t(apply(la_mt_mn, 1, function(x, adj=0.1/50){
  res1<-ifelse( x[1]==0, (x[3]+adj)/(x[1]+adj), x[3]/x[1] )
  mn1<-x[1];mn2<-x[3]  
  res2<- ifelse( min(mn1,mn2)<adj, log2((mn1+adj)/(mn2+adj)), log2(mn1/mn2))
  res<-c(res1,res2)
  return(res)
}))
# histogram on FD
plot_data<-1/la_jump[,1]
plot_data[plot_data>20]<-20
par(cex=2,las=1,mar=c(4,4,1,1),lwd=5,pch=16)
hist( plot_data, main='Mouse', ylab='Number of cell types', xlab='Fold difference', breaks=0:20)
text( 18,-3.7, label='>', xpd=NA)

#' LIN28A expression in a broad time windown in zebrafish, frog, and mouse.
# Zebrafish: White et al. 2017
raw<- read.table('../../White/elife-30860-supp3-v1.tsv',sep='\t',as.is=T,header=T,comment.char='',quote='')
lin<-c('ENSDARG00000004328','ENSDARG00000016999','ENSDARG00000052511') # Zebrafish lin28A,lin28a, lin28B
st<- sapply( colnames(raw)[-(1:8)], function(x){ substring(x,1,nchar(x)-2) })
st_anno<- c(paste( c(0,0.75,2.25,3,4.3,5.25,6,8,10.3,16,19,24,30,36),'hpf',sep=''),'2d','3d','4d','5d')
lin_mx<-raw[raw[,1] %in% lin, -(2:8)]
rownames(lin_mx)<-lin_mx[,1]
lin_mx<-lin_mx[,-1]
lin_exp<- t(apply( lin_mx, 1, function(x){ tapply(x, st, mean) }))[,unique(st)]
#pdf('../result/lin28_exp_f2_b.pdf',width=10, height=5)
par(cex=2,las=1,mar=c(2,4,1,2),lwd=5,pch=16)
  plot( 1:ncol(lin_exp), lin_exp[1,], type='l', frame=F, main='', xlab='',ylab='TPM', ylim=c(0, max( lin_exp[!is.na(lin_exp)] ) ) ,xaxt='n',xpd=NA, col='black') 
  points(1:ncol(lin_exp), lin_exp[1,], cex=1,xpd=NA,col='black')
  text( 1:ncol(lin_exp), (-par('usr')[4]+par('usr')[3])*.13, label=st_anno, xpd=NA, cex=.8,srt=90)
  legend( 'topright', legend=c('DR-LIN28A'), col=c('black'),lty=c(1),xpd=NA,bty='n')
#dev.off()

# Frog: Tan et al. 2013
# Load frog gene matrix after remapping it to Ensembl100 (see 'remap_Tan_frog' in this repository)
load('../../Tan_frog/ensembl100/code/xt_data.rdata')
lina_exp<- xt_all_exp[ c('ENSXETG00000012324','ENSXETG00000035037'), ] # frog lin28a and lin28b
#pdf('../result/newgenome_lin28a_exp_f2_b.pdf',width=10, height=5)
par(cex=2,las=1,mar=c(2,4,1,2),lwd=5,pch=16)
  plot( 1:ncol(lina_exp), lina_exp[1,], type='l', frame=F, main='', xlab='',ylab='TPM', ylim=c(0, max( lina_exp[!is.na(lina_exp)] ) ) ,xaxt='n',xpd=NA,xlim=c(1,23)) # 
  points(1:ncol(lina_exp), lina_exp[1,], cex=1,xpd=NA)
#  legend( 'topright', legend=c('XT-LIN28A'), col=c('black'),lty=c(1),xpd=NA,bty='n')
  text( 1:23, (-par('usr')[4]+par('usr')[2])*.12, label= colnames(xt_all_exp), xpd=NA, cex=.5,srt=90)
#dev.off()

# Mouse: Boroviak et al. 2015, Pijuan-Sala et al. 2019, Cao et al. 2019
# Because mouse data are from three different sources, we do normalizaton between them and do not draw line between datasets. The Boroviak and Pijuan-Sala datasets were normalized by ribosome genes. The Pijuan-Sala and Cao datasets were normalized by Lin28a expression at E8.5 (Pijuan-Sala) and E9.5 (Cao), because we found ribosome genes have large variation between cytosol (Pijuan-Sala) and nucleis (Cao).
# Boroviak data
raw<-read.table('~/baolab/xuy/human/lin28/Boroviak_mouse/table_s1.txt',sep='\t',as.is=T,header=T)
st<-sapply(colnames(raw)[c(4:18)],function(x){ paste(strsplit(x,split='.',fixed=T)[[1]][1:2],collapse='.') }) 
lin_exp<-rbind( raw[ raw[,1]=='Lin28a', 4:18], raw[ raw[,1]=='Lin28b', 4:18] )
rownames(lin_exp)<- c('Lin28a', 'Lin28b')
m1_lin<- t(apply( lin_exp,1,function(x){
  tapply(x, st, mean)
}))
(hmrb<-read.table('~/baolab/xuy/human/lin28/hmrb.txt',sep='\t',as.is=T,row.name=1)) # mouse ribosome genes
rb_exp<-t(sapply(hmrb[,2],  function(x){
  res<-raw[ raw[,1]==x, 4:18]
  return(res)
}))
m1_rb<- t(apply( rb_exp,1,function(x){
  tapply(as.numeric(x), st, mean)
}))
# Pijuan-Sala data
load('~/baolab/xuy/human/lin28/Pijuan-Sala_mouse/code/m2_lin.rdata')
m2_rb<- t(sapply(1:37, function(x){
  res<-as.numeric(scan(paste( '~/baolab/xuy/human/lin28/Pijuan-Sala_mouse/result/rb/', hmrb[x,1],'.txt',sep=''),what=character(0))[-1])
}))
rownames(m2_rb)<- hmrb[1:37,1]
m2_rb_mn<- apply(m2_rb,1,mean)
m1_rb_mn<- apply(m1_rb,1,mean)
# Find the normalization factor between Boroviak and Pijuan-Sala by linear relationship between ribosome genes from two datasets
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3,pch=16)
plot(m1_rb_mn[1:37],m2_rb_mn, xlab='Boroviak RB', ylab='Pijuan-Sala RB',main='', frame=F)
m12_lm<-lm( m2_rb_mn~m1_rb_mn[1:37] )
abline(m12_lm,col='red')
m12_fd<-summary(m12_lm)$coefficients[2,1]
m2_lin_r1<- m2_lin/m12_fd
# Cao data
m3_lin<- read.table('~/baolab/xuy/human/Shendure/shendure_lin28b_stage_mean.txt',sep=' ',as.is=T,row.name=1)
colnames(m3_lin)<- paste('E',c('9.5','10.5','11.5','12.5','13.5'),sep='')
lin_mg<- cbind( m1_lin, m2_lin_r1, m3_lin/(m3_lin[1,1]/m2_lin_r1[1,9]) ) # normalize Pijuan-Sala and Cao data by lin28a expression at E8.5 and E9.5
# plot Lin28a expression
plot_data<-lin_mg[1,]
par(cex=2,las=1,mar=c(2,4,1,2),lwd=5,pch=16)
plot( 1:ncol(plot_data), plot_data[1,], type='l', frame=F, main='', xlab='',ylab='relative expression',xaxt='n',xpd=NA, col='white',xlim=c(1,length(plot_data))) 
  points(1:ncol(plot_data), plot_data[1,], cex=1,xpd=NA,col='black')
  for(ind in list(1:4,5:13,14:18)) lines( (1:ncol(plot_data))[ind], plot_data[1,ind], col='black')
  text( 1:ncol(plot_data), (-par('usr')[4]+par('usr')[3])*.1, label= colnames(plot_data), xpd=NA, cex=.8,srt=90)
# legend( 'topleft', legend=c('MM-LIN28A'), col=c('black'),lty=c(1,1),xpd=NA,bty='n')


# rmarkdown::render('lin28_fig5ad.r')