#' LIN28A Figure 4: LIN28A appears to promote the translation of systemically changing genes. 

#' Compare systemically changing genes (human) to LIN28A targets and let-7 targets.
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
knitr::opts_chunk$set(warning = FALSE)
# Load data
load('temporal2.rdata') # systemically changing genes from 'lin28.r'
sapply(hs_glog<- c( list(unlist(hs_glo[1:2])), list(unlist(hs_glo[3:4])), list(ntg) ),length) # systemically down-regulated (from stage 3 to 4), systemically up-regulated, systemically unchanged genes
names(hs_glog)<-c('sys DN','sys UP','sys NC')
# LIN28A binding in human ESCs: Wilbert et al. 2012
length(yeo_hes<-unique(unlist(get_id(read.table('~/baolab/xuy/human/lin28/targets/gene_yeo_2012/hg18_refgene_lin28ES.txt',sep='\t',as.is=T)[,4]))))

# A function to test enrichment
test_enrich_logic2<- function( diag, test, total){ # hypergeometric test
  y<- diag[ diag %in% total]
  x<- test[ test %in% total]
    test_mx<- matrix( c( sum( x %in% y), sum( x %in% setdiff(total,y) ), sum( y %in% setdiff(total,x)), sum( setdiff(total,x) %in% setdiff(total,y) ) ), nr=2)
    pv<- fisher.test(test_mx)$p.value
    return( c(pv, sum(x%in%y),length(y), sum(x%in%y)/length(x), length(y)/length(total) ) )
}
tempg_glogl_yh<- sapply(hs_glog, function(x){ test_enrich_logic2( diag=yeo_hes, test= x, total=intersect(in_gene3,hg18_allg) ) })
plot_mx<- rbind( sapply(hs_glog,length)-tempg_glogl_yh[2,] , tempg_glogl_yh[2,])
plot_col<- rep('black', ncol(plot_mx))
plot_col[tempg_glogl_yh[1,]<0.001 & tempg_glogl_yh[4,]>tempg_glogl_yh[5,] ]<-'red' # significant cutoff p < 0.001
plot_col[tempg_glogl_yh[1,]<0.001 & tempg_glogl_yh[4,]<tempg_glogl_yh[5,] ]<-'blue'

# A function to plot the enrichment result
plot_enrich_tar<-function(mx, back, ratio,name, path='../result/',title='',perc_col,wd=7,adj=25,lab_cex=.5, main_cex=1){
par(cex=2,las=1,mar=c(2,4,4,0),lwd=1,pch=16 )
tmp<-barplot(mx, col=c('grey','red') , ylab='Number of genes',main=title, names.arg=rep('',ncol(mx)), cex.main=main_cex)
text( tmp, par('usr')[3]-(par('usr')[4]-par('usr')[3])/8, colnames(mx), xpd=NA,cex=lab_cex)
text( tmp, par('usr')[4]*1.15, label=paste(round(ratio,2),sep=''), xpd=NA,cex=1, col=perc_col)
text( par('usr')[1]-(par('usr')[2]-par('usr')[1])/8, par('usr')[4]*1.15, label= round( back, 2), cex=1, xpd=NA, col='darkgrey')
}
plot_enrich_tar(mx=plot_mx, back=tempg_glogl_yh[5,1], ratio=tempg_glogl_yh[4,], name='',title='LIN28A pull down',perc_col=plot_col,wd=5,adj=15,lab_cex=.8)

# let-7 targets: TargetScan
ts_let7<- unlist(unique(get_id( read.table(file='~/BaoLab/xuy/human/lin28/targets/let-7/targetscan/TargetScan7.2__let-7-5p_98-5p.predicted_targets.txt', sep='\t', as.is=T,header=T,comment.char='',quote='',fill=T)[,1] )))
tempg_glogl_l7ts5<- sapply( hs_glog, function(x){ test_enrich_logic2( diag=ts_let7, test= x, total=hg18_allg ) }) # background: whole genome
plot_mx<- rbind( sapply(hs_glog,length)-tempg_glogl_l7ts5[2,] , tempg_glogl_l7ts5[2,])
plot_col<- rep('black', ncol(plot_mx))
plot_col[tempg_glogl_l7ts5[1,]<0.001 & tempg_glogl_l7ts5[4,]>tempg_glogl_l7ts5[5,] ]<-'red'
plot_col[tempg_glogl_l7ts5[1,]<0.001 & tempg_glogl_l7ts5[4,]<tempg_glogl_l7ts5[5,] ]<-'blue'
plot_enrich_tar(mx=plot_mx, back=tempg_glogl_l7ts5[5,1], ratio=tempg_glogl_l7ts5[4,], name='',title='let-7 targets',perc_col=plot_col,wd=5,adj=15,lab_cex=.8)

#' Compare down-regulated genes in systemically changing genes (human&mouse) to mRNA and ribosome occupancy upon Lin28a KD
load('../../Shendure_mouse/code/mm_dp_glo.rdata')
hs_glo_ft<- mapply( function(x,y){ intersect(ge2an[x],y) }, x=list( unlist(hs_glo[1:2]), unlist(hs_glo[3:4])), y= mm_dp_glo )
# mRNA and ribosome occupancy upon Lin28a KD: Cho et al. 2012
lakd_rb<- read.table('../../Cho_lin28/mmc5.txt',sep='\t',as.is=T,comment.char='',fill=T,quote='',header=T) # Cho's data
lakd_rb<- lakd_rb[!duplicated(lakd_rb[,2]),]
lakd_rb_mn<-cbind(log2(lakd_rb[,18]+1), lakd_rb[,c(13)], lakd_rb[,21], lakd_rb[,12] ) # mRNA expression, LFD of ribosome occupancy, LFD of mRNA
rownames(lakd_rb_mn)<- lakd_rb[,2]

# A function for MA plot
plot_ma_with_glo2<- function(mx, ind1, ind2, glo, glo2, file_name='', name='', ylm=3, tar, is_p1=F,path='cross_species',title='',cols=c('blue','red'),x_name='',y_name='', is_red=F){
if(x_name=='')  x_name<-paste('Mean (',name,'&WT)',sep='')
if(y_name=='') y_name<-paste('log2 FD (',name,'/WT)',sep='') else if(y_name=='none') y_name<-''
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16)
x<- mx[,ind1]; y<-mx[,ind2]
plot( x,y, xlab=x_name ,ylab=y_name,main=title,frame=F,col='grey', cex=.2, ylim=c(-ylm,ylm),cex.main=1, xaxt='n')
axis(1,at=seq(0, round(max(x[!is.na(x)])), length.out=4))
res<-c(0,0)
if(is_red) len<-2 else len<-1
for(i in len:len ){
  plot_col<- cols[i]
  g<- intersect( glo[[i]], rownames(mx) )
  res[i]<-mean(mx[g,ind2])
  if(i==1) p1<-mx[g,ind2]
  else if(i==2) p2<-mx[g,ind2]
  x<- mx[g,ind1]; y<-mx[g,ind2]
  points( x, y, col=plot_col, cex=.5)
}
}
plot_ma_with_glo2( mx=lakd_rb_mn, ind1=1, ind2=2, glo=sapply( hs_glo_ft,function(x){ setdiff(unique(hs_orth[hs_orth[,2]%in%x,8]),'') }), file_name='lin28a_kd',name='RB_hs&mm_glo', path='cross_species/glo_4sp',title='' ,ylm=3, x_name='log2 RPKM of mRNA', y_name='none',is_red=F ) 
plot_ma_with_glo2( mx=lakd_rb_mn, ind1=1, ind2=3, glo=sapply( hs_glo_ft,function(x){ setdiff(unique(hs_orth[hs_orth[,2]%in%x,8]),'') }), file_name='lin28a_kd',name='mRNA_hs&mm_glo', path='cross_species/glo_4sp',title='' ,ylm=3, x_name='log2 RPKM of mRNA', y_name='log2 FD of mRNA',is_red=F ) 

# rmarkdown::render('lin28_fig4.r')