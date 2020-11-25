#' Fig. S6: digit in situ in forelimb
#' For any gene, one can plot its average normalized UMIs in each limb domain. Also see https://organroot.shinyapps.io/base/

#' Load gene expression data and limb domain data
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
library(plotrix)
load(file='limb_do_mn.rdata')
load(file='limb_do2tj.rdata')
load('ge2an.rdata') # gene annotation
    fl7_dc<-lapply( names(do2tj)[1:6],function(x){
      res<-read.table(paste('figs6/st1/',x,'.csv',sep=''), sep=',',as.is=T,skip=1)[,6:7]
      res[,2]<- 380-res[,2]
      return(res)
    })
    fl0_dc<-lapply( names(do2tj)[7:15],function(x){
      res<-read.table(paste('figs6/st2/',x,'.csv',sep=''), sep=',',as.is=T,skip=1)[,6:7]
      res[,2]<- 380-res[,2]
      return(res)
    })
    fl3_dc<-lapply( names(do2tj)[16:29],function(x){
      res<-read.table(paste('figs6/st3/',x,'.csv',sep=''), sep=',',as.is=T,skip=1)[,6:7]
      res[,2]<- 430-res[,2]
      return(res)
    })

#' A function to plot gene expression in domains of 3 limb buds
source('../../scatter_fill2.r')
plot_digit<-function(x, path='digit/',col_grad=c(colorRampPalette(c("white","red"))(20)), cut_min=0.3,zlim_mx=2,is_cs12=F,is_cs13=F,is_cs15=F ){
  if(is_cs12){
  val<- do_mn[ x, do2tj[1:6] ]
  val[val<cut_min]<-0
  if(sum(val>cut_min)==0) plot_col<-rep('white',length(val))
  else{
    plot_col<- scatter_fill2(z=val,col=col_grad, only_col=T, zlim=c(0,zlim_mx))
    plot_col<-plot_col[1:length(val)]
  }
  st<-'cs12'
#  pdf(paste(path,st,'_',ge2an[x],'.pdf',sep=''),width=3.74,height=5.46)
  par(cex=2,las=1,mar=c(0,0,0,0),lwd=3,pch=16)
  plot(0,0,xlab='',ylab='',main='',frame=F, xaxt='n',yaxt='n', xlim=c(0,450), ylim=c(0,380),col='white')
  for(i in 1:length(fl7_dc)){
    polygon( x=fl7_dc[[i]][,1], y=fl7_dc[[i]][,2], col=NA, border='black',lwd=4)
    polygon( x=fl7_dc[[i]][,1], y=fl7_dc[[i]][,2], col=plot_col[i], border=plot_col[i])
  }
#  dev.off()  
  }
  if(is_cs13){
  val<- do_mn[ x, do2tj[7:15] ]
  val[val<cut_min]<-0
  if(sum(val>cut_min)==0) plot_col<-rep('white',length(val))
  else{
    plot_col<- scatter_fill2(z=val,col=col_grad, only_col=T, zlim=c(0,zlim_mx))
    plot_col<-plot_col[1:length(val)]
  }
  st<-'cs13'
#  pdf(paste(path,st,'_',ge2an[x],'.pdf',sep=''),width=4.01,height=5.19) # 4.01*5.46/5.19= 4.218613
  par(cex=2,las=1,mar=c(0,0,0,0),lwd=2,pch=16)
  plot(0,0,xlab='',ylab='',main='',frame=F, xaxt='n',yaxt='n', xlim=c(0,450), ylim=c(0,380),col='white')
  for(i in 1:length(fl0_dc)){
    polygon( x=fl0_dc[[i]][,1], y=fl0_dc[[i]][,2], col=NA, border='black',lwd=4)
    polygon( x=fl0_dc[[i]][,1], y=fl0_dc[[i]][,2], col=plot_col[i], border=plot_col[i])
  }
#  dev.off()
  }
  if(is_cs15){
  val<- do_mn[ x, do2tj[16:29] ]
  val[val<cut_min]<-0
  if(sum(val>cut_min)==0) plot_col<-rep('white',length(val))
  else{
    plot_col<- scatter_fill2(z=val,col=col_grad, only_col=T, zlim=c(0,zlim_mx))
    plot_col<-plot_col[1:length(val)]
  }
  st<-'cs15'
#  pdf(paste(path,st,'_',ge2an[x],'.pdf',sep=''),width=5.99,height=5.94) # 5.99*5.46/5.94=5.50596
  par(cex=2,las=1,mar=c(0,0,0,0),lwd=3,pch=16)
  plot(0,0,xlab='',ylab='',main='',frame=F, xaxt='n',yaxt='n', xlim=c(0,450), ylim=c(0,430),col='white')
  for(i in 1:length(fl3_dc)){
    polygon( x=fl3_dc[[i]][,1], y=fl3_dc[[i]][,2], col=NA, border='black',lwd=4)
    polygon( x=fl3_dc[[i]][,1], y=fl3_dc[[i]][,2], col=plot_col[i], border=plot_col[i])
  }
#  dev.off()
  }
}
#' plot MSX2 as an example
plot_digit( names(ge2an)[ge2an=='MSX2'], is_cs12=T, is_cs13=T, is_cs15=T)
# rmarkdown::render('figs6.r')
