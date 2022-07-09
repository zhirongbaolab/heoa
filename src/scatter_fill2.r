# 2019.1.22, modify for multiple plots
scatter_fill2 <- function (x, y, z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 20,title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, 
                          axes = TRUE, frame.plot = axes, is_0f=F, plot_cex1,plot_cex2,col=-1,is_axis='n',is_legend=T,xlab='',ylab='',main_cex=1,only_col=F,...) 
{
# choose colors to interpolate
#print(zlim)
levels <- seq(zlim[1],zlim[2],length.out = nlevels)
if(col[1]==-1) col <- c('grey',colorRampPalette(c("dark green","yellow","red"))(nlevels))
colz <- col[cut(c(z,zlim[2]),nlevels)]  
if(only_col) return(colz)
#   
   # points
   plot(x,y,type = "n",xaxt=is_axis,yaxt=is_axis,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,bty="n",main=title,cex.main=main_cex)
   if(!is_0f){
     ind<- sample.int(length(x))
     points(x[ind],y[ind],col = colz[ind],xaxt='n',yaxt='n',xlab="",ylab="",bty="n",cex=plot_cex1,...)
   }else{
     points(x[z==0],y[z==0],col = colz[z==0],xaxt='n',yaxt='n',xlab="",ylab="",bty="n",cex=plot_cex1,...)
     points(x[z!=0],y[z!=0],col = colz[z!=0],xaxt='n',yaxt='n',xlab="",ylab="",bty="n",cex=plot_cex2,...)
   }
#   legend( par('usr')[2]+(par('usr')[2]-par('usr')[1])/20,par('usr')[4], signif(c(levels[1],min(z[z>levels[1]]),seq(levels[2],zlim[2],length.out=3)[2:3]),digit=2), col=c('grey','dark green',"yellow","red"),bty='n',pch=16,xpd=NA)
  if(is_legend) legend( par('usr')[2]-(par('usr')[2]-par('usr')[1])/20,par('usr')[4], signif(c(levels[1],zlim[2]),digit=2), col=col[c(1,length(col))] ,pch=16,xpd=NA)
  return(colz)
}