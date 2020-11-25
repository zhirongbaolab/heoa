####### LIN28 expression delineates developmental stages in vertebrate embryogenesis

####### Custom functions used in analysis

# calculate expressed genes
get_exp<-function(mx, cut1=10, cut2=5, ind){
  res1<- apply(mx[,ind],1,max)>= cut1
  inds<- lapply( 1:(length(ind)-1), function(x){ ind[x:(x+1)] })
  loose_res<- lapply( inds, function(x){ apply( mx[,x]>=cut2, 1, sum)==2 })
  loose_res<- Reduce('|',loose_res)
  res<- res1 | loose_res
  return(res)
}

# calculate DEGs between adjecent time points
cal_adj_deg<-function(mx, cut_fd=1/2, cut_exp=10, isg=F){
  res<-sapply( 1:(ncol(mx)-1), function(x){
    fd<- apply( mx[,x:(x+1)], 1, function(y){ min(y)/max(y) })
    mexp<- apply( mx[,x:(x+1)], 1, function(y){ max(y) })
    whichm<-apply( mx[,x:(x+1)], 1, function(y){ which.max(y) })
    num<- c( sum(fd< cut_fd & mexp>= cut_exp & whichm==2), sum(fd< cut_fd & mexp>= cut_exp & whichm==1) )
    if(isg){
      gres<- rep(0,nrow(mx))
      gres[fd< cut_fd & mexp>= cut_exp & whichm==2]<-1
      gres[fd< cut_fd & mexp>= cut_exp & whichm==1]<- -1
      return(gres)
    }
    return(num)
  })
  return(res)
}
# calculate FDs between adjecent time points
cal_adj_trd<-function(mx, cut_fd=1/2, cut_exp=10){
  res<-sapply( 1:(ncol(mx)-1), function(x){
    fd<- apply( mx[,x:(x+1)], 1, function(y){ min(y)/max(y) })
    mexp<- apply( mx[,x:(x+1)], 1, function(y){ max(y) })
    whichm<-apply( mx[,x:(x+1)], 1, function(y){ which.max(y) })
    sign1<- ifelse( whichm==2, -1, 1)
    fd2<- ifelse( mexp>=cut_exp&fd<cut_fd, fd*sign1, 1)
    fd2[is.na(fd2)]<-1
    return(fd2)
  })
  return(res)
}
# calculate FDs between adjecent time points: log2+1
cal_adj_trd2<-function(mx, cut_fd=1/2, cut_exp=10){
  res<-sapply( 1:(ncol(mx)-1), function(x){
    fd<-log2(mx[,x]+1)-log2(mx[,x+1]+1)
    return(fd)
  })
  return(res)
}

# calculate DEGs between adjecent time points, considering Wilcox p value
cal_wilp<-function(mx1,mx2){
  res<-mapply(function(x,y){
    wilcox.test(x,y)$p.value
  }, x=split( mx1, 1:nrow(mx1)), y=split( mx2, 1:nrow(mx2)) )
}
cal_ttestp<-function(mx1,mx2){
  res<-mapply(function(x,y){
    t.test(x,y)$p.value
  }, x=split( mx1, 1:nrow(mx1)), y=split( mx2, 1:nrow(mx2)) )
}
# calculate difference without Wilcox: for scRNA-seq
cal_wilp2<-function(mx1, mx2, genes,is_sam=F, cut_mean=0.1, input_mn=F, small_val=-1, is_pv=F ){
  if(!input_mn){
    mn1<- apply(mx1,1,mean)
    mn2<- apply(mx2,1,mean)
  }else{
    mn1<-mx1
    mn2<-mx2
  }
  min_mn<- apply( cbind( mn1, mn2 ), 1, min)
  if(small_val==-1) small_val<- cut_mean/5
  fd<- ifelse( min_mn< small_val, log2((mn1+small_val)/(mn2+small_val)), log2(mn1/mn2))
  fd[ apply( cbind( mn1, mn2 ), 1, max) < cut_mean ]<-0
  if(is_pv){
   pv<- mapply( function(x,y){
     if( mean(x)<cut_mean & mean(y)<cut_mean ) return(1)
     else return(wilcox.test( x, y)$p.value)
   },x=split(mx1,1:nrow(mx1)), y=split(mx2,1:nrow(mx2)) )
   return( cbind( fd, pv ))
  } else return(fd)
}


# the start of expression
get_start<-function(x, cut1=10, cut2=5){
  st1<- min(which(x>=cut1))
  if(is.na(st1)) st1<-length(x)
  st2<-0
  for( st2 in 1:(length(x)-1)){
    if( x[st2]>=cut2&x[st2+1]>=cut2) break
  }
  if( st2== (length(x)-1) & (x[st2]<cut2|x[st2+1]<cut2) ) st2<-length(x)
  if( st2<length(x) ) st2<-st2+1
  return(min(st1,st2))
}

# order a gene list in heatmap
order_glist<-function(gs, mx){
  god<- list('')
  for( i in 1:length(gs)){
    tmp<-heatmap3(mx[ gs[[i]] , ] ,labRow=NA,scale='none',dendrogram='none',trace='none',Rowv=T,Colv=F,symkey=F,density.info="none",keysize=1,col=colorRampPalette(c("blue","white","red"))(499),color_key_label='log2 FD',color_key_label_cex=1,margins=c(3,3),color_key_axis_cex=1,key_mar=c(5, 1, 2, 1),labCol='',labRow_pos=c(2,4),sepwidth=c(0.1,0.1),sepcolor='black',cexRow=.4,lhei=c(1,5) )
    god<- c( god, list(gs[[i]][rev(tmp$rowInd)]) )
  }
  return(god[-1])
}

# compare two lists of genes with trend
cmp_trend<-function(trd1, trd2, ind1, ind2, orth, hs21, hs22){
apply( trd1,2, function(x){
  res<-apply( trd2, 2, function(y){
    g1<- rownames(trd1)[abs(x)<1/2]
    g2<- rownames(trd2)[abs(y)<1/2]
    total_hs<- setdiff(unique(orth[ orth[,ind1]%in%g1 | orth[,ind2]%in%g2,1]),'')
    if(length(total_hs)==0) return(0)
    total_hs_g2<- sapply( total_hs, function(z){
      if(length(z)==0) return(0)
      fd<- y[rownames(trd2)%in%hs21[z]]
      fd_sum<- c( sum(fd>0 & fd!=1), sum(fd<0) )
      if( fd_sum[1]>fd_sum[2] ) return(1) else if( fd_sum[1]<fd_sum[2] ) return(-1) else return(0)      
    })
    total_hs_g1<- sapply( total_hs, function(z){
      if(length(z)==0) return(0)
      fd<- x[rownames(trd1)%in%hs22[z]]
      fd_sum<- c( sum(fd>0 & fd!=1), sum(fd<0) )
      if( fd_sum[1]>fd_sum[2] ) return(1) else if( fd_sum[1]<fd_sum[2] ) return(-1) else return(0)      
    })
    return( sum(total_hs_g2*total_hs_g1>0)/length(total_hs) )
  })
  print('x')
  return(res)
})
}

# get global genes from scRNA-seq data
get_glo<-function(x, rm_tp='', broad_ratio=3/4 ){ # increase 'broad_ratio' from 1/2 to 3/4, compare to human scRNA-seq
  x<-x[,setdiff(colnames(x),rm_tp) ]
  cut_num<-  ceiling(broad_ratio*ncol(x) )
  conflict<- floor( 1/3*ncol(x) )
  dn<-rownames(x)[apply(x,1,function(y){
    is_glo<-sum(y>log2(1.2))>cut_num & (sum(y>1)>0) & (sum(y< -1)==0) 
    if( sum(y>0)>=conflict & sum(y<0)>=conflict ) is_glo<-F
    return(is_glo)
  })]
  up<-rownames(x)[apply(x,1,function(y){
    is_glo<-sum(y< -log2(1.2))> cut_num & (sum(y< -1)>0 ) & (sum(y>1)==0) 
    if( sum(y>0)>=conflict & sum(y<0)>=conflict ) is_glo<-F
    return(is_glo)
  })]
  return(list(dn,up))
}
# average LFD in scRNA-seq
cal_sc_lfd<- function(diff, cut_lfd=log2(1.2)){
apply( diff, 1, function(x){ xx<-x[abs(x)>=cut_lfd ]; if(length(xx)==0) return(0) else return( mean(xx)) })
}

# get specific genes on scRNA-seq data
get_spec<-function(x, rm_tp='', spec_ratio=1/5 ){ 
  x<-x[,setdiff(colnames(x),rm_tp) ]
  cut_num<- floor( spec_ratio*ncol(x) )
  res<-rownames(x)[apply(x,1,function(y){
    is_spec<-sum(abs(y)>0)<=cut_num & (sum( abs(y) >1)>0)  
    return(is_spec)
  })]
  return(res)
}

# read global genes from scRNA-seq
read_gene<- function(x, path='../../Briggs/result/up_phase/'){
  dn<- scan( paste( path,x, '_DN.txt',sep=''),what=character(0))
  up<- scan( paste( path,x, '_UP.txt',sep=''),what=character(0))
  return(list(dn,up))
}

# plot dyanmic for a group of genes
plot_dynamic<- function( mx, gene, st_lab, st_col, path, file_name, name, lab_pos=.13, ret_val=F){
plot_mx<-mx[gene,]
plot_mx<- t(apply(plot_mx, 1, function(x){ (x-min(x))/(max(x)-min(x)) }))
if(ret_val) return(plot_mx)
file_name<- paste(path, '/', file_name,'.pdf', sep='')
pdf(file_name, width=10,height=5)
par(cex=2,las=1,mar=c(4,4,2,1),lwd=3,pch=16)
plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=name, frame=F, xaxt='n')
for( i in 2:nrow(plot_mx)) lines( 1:ncol(plot_mx), plot_mx[i,] )
text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos, label=st_lab, xpd=NA, cex=.8,srt=90, col=st_col )
dev.off()
}

# plot dyanmic for a group of genes one by one in 1 pdf
plot_one_dynamic<- function( mx, gene, st_lab, st_col, path, file_name, name, lab_pos=.13, ret_val=F){
plot_mx<-mx[gene,]
plot_mx<- t(apply(plot_mx, 1, function(x){ (x-min(x))/(max(x)-min(x)) }))
if(ret_val) return(plot_mx)
file_name<- paste(path, '/', file_name,'.pdf', sep='')
pdf(file_name, width=5,height=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=3,pch=16, mfrow=c(5,1))
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=name[1],, frame=F, xaxt='n')
  text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos, label=st_lab, xpd=NA, cex=.8,srt=90, col=st_col )
for( i in 2:nrow(plot_mx)){
  plot( 1:ncol(plot_mx), plot_mx[i,], type='l', xlab='', ylab='relative expression',main=name[i],, frame=F, xaxt='n')
  text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos, label=st_lab, xpd=NA, cex=.8,srt=90, col=st_col )
}
dev.off()
}

# plot dyanmic for a group of genes batch by batch in 1 pdf
plot_batch_dynamic<- function( mx, gene, st_lab, st_col, path, file_name, name, lab_pos=.13, ret_val=F){
plot_mxs<-lapply( gene,function(y){
  plot_mx<-mx[y,]
  plot_mx<- t(apply(plot_mx, 1, function(x){ (x-min(x))/(max(x)-min(x)) }))
  return(plot_mx)
})
if(ret_val) return(plot_mxs)
file_name<- paste(path, '/', file_name,'.pdf', sep='')
pdf(file_name, width=5,height=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16, mfrow=c(5,1))
  plot_mx<-plot_mxs[[1]]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=1, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
  text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos, label=st_lab, xpd=NA, cex=.8,srt=90, col=st_col )
for( i in 2:length(plot_mxs)){
  plot_mx<-plot_mxs[[i]]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=i, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
}
dev.off()
}

# plot dyanmic for a group of genes batch by batch in 1 pdf (dr and xt)
plot_batch_dynamic_drxt<- function( mx, gene, st_lab, st_col, path, file_name, name, lab_pos=.13, ret_val=F, st_lab2, st_col2, lab_pos2=.2){
plot_mxs<-mx
if(ret_val) return(plot_mxs)
file_name<- paste(path, '/', file_name,'.pdf', sep='')
pdf(file_name, width=7,height=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16, mfrow=c(5,2))
  plot_mx<-plot_mxs[[1]][,1:18]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=1, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
  text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos, label=st_lab, xpd=NA, cex=.8,srt=90, col=st_col )
  plot_mx<-plot_mxs[[1]][,-(1:18)]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=1, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
  text( 1:ncol(plot_mx), (-par('usr')[4]+par('usr')[3])*lab_pos2, label=st_lab2, xpd=NA, cex=.8,srt=90, col=st_col2 )

for( i in 2:length(plot_mxs)){
  plot_mx<-plot_mxs[[i]][,1:18]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=i, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
  plot_mx<-plot_mxs[[i]][,-(1:18)]
  plot( 1:ncol(plot_mx), plot_mx[1,], type='l', xlab='', ylab='relative expression',main=i, frame=F, xaxt='n')
  for( j in 2:nrow(plot_mx)) lines(1:ncol(plot_mx), plot_mx[j,])
}
dev.off()
}


# MA plot
plot_ma_with_glo<- function(mx, ind1, ind2, glo, glo2, file_name='', name='', ylm=3, tar, is_p1=F,path='cross_species', y_lab=''){
if(y_lab=='') y_lab=paste('log2 FD (',name,'/WT)',sep='')
pdf( paste('../result/',path,'/',file_name,'_',name, '_maplot_with_glo.pdf',sep='') ,height=7,width=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16 )
x<- mx[,ind1]; y<-mx[,ind2]
plot( x,y, xlab= paste('Mean (',name,'&WT)',sep=''),ylab=y_lab,main='',frame=F,col='grey', cex=.2, ylim=c(-ylm,ylm))
for(i in 1:4){
  plot_col<- ifelse( i<=2, 'blue','red')
  g<- intersect( glo[[i]], rownames(mx) )
  x<-mx[g,ind1]; y<-mx[g,ind2]
  points( x, y, col=plot_col, cex=.4)
}
dev.off()
print('finish 1')
if(is_p1) return(NA)
# conserved genes in human & mouse
pdf( paste('../result/',path,'/',file_name,'_',name, '_maplot_with_glo_con.pdf',sep=''), height=7,width=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16 )
x<- mx[,ind1]; y<-mx[,ind2]
plot( x,y, xlab= paste('Mean (',name,'&WT)',sep=''),ylab=paste('log2 FD (',name,'/WT)',sep=''),main='',frame=F,col='grey', cex=.2, ylim=c(-ylm,ylm))
for(i in 1:2){
  plot_col<- ifelse( i<=1, 'blue','red')
  g<- intersect( glo2[[i]], rownames(mx))
  x<-mx[g,ind1]; y<-mx[g,ind2]
  points( x, y, col=plot_col, cex=.4)
}
dev.off()
print('finish 2')
# global genes with lin28a target
pdf( paste('../result/',path,'/',file_name,'_',name, '_maplot_with_glo_tar.pdf',sep='') ,height=7,width=7)
par(cex=2,las=1,mar=c(4,4,1,1),lwd=1,pch=16 )
x<- mx[,ind1]; y<-mx[,ind2]
plot( x,y, xlab= paste('Mean (',name,'&WT)',sep=''),ylab=paste('log2 FD (',name,'/WT)',sep=''),main='',frame=F,col='grey', cex=.2, ylim=c(-ylm,ylm))
for(i in 1:4){
  plot_col<- ifelse( i<=2, 'blue','red')
  g<- intersect( intersect( glo[[i]], tar), rownames(mx))
  x<-mx[g,ind1]; y<-mx[g,ind2]
  points( x, y, col=plot_col, cex=.4)
}
dev.off()
print('finish 3')
}

# for list of DN and UP
plot_ma_with_glo2<- function(mx, ind1, ind2, glo, glo2, file_name='', name='', ylm=3, tar, is_p1=F,path='cross_species',title='',cols=c('blue','red'),x_name='',y_name='', is_red=F){
if(x_name=='')  x_name<-paste('Mean (',name,'&WT)',sep='')
if(y_name=='') y_name<-paste('log2 FD (',name,'/WT)',sep='') else if(y_name=='none') y_name<-''
pdf( paste('../result/',path,'/',file_name,'_',name, '_maplot_with_glo.pdf',sep='') ,height=5,width=5)
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
dev.off()
return( res ) 
}

# function to plot enrichment result
plot_enrich_tar<-function(mx, back, ratio,name, path='../result/',title='',perc_col,wd=7,adj=25,lab_cex=.5, main_cex=1){
pdf(paste(path,name,'.pdf',sep='') ,height=5,width=wd)
par(cex=2,las=1,mar=c(2,4,4,0),lwd=1,pch=16 )
tmp<-barplot(mx, col=c('grey','red') , ylab='Number of genes',main=title, names.arg=rep('',ncol(mx)), cex.main=main_cex)
text( tmp, par('usr')[3]-(par('usr')[4]-par('usr')[3])/8, colnames(mx), xpd=NA,cex=lab_cex)
text( tmp, par('usr')[4]*1.15, label=paste(round(ratio,2),sep=''), xpd=NA,cex=1, col=perc_col)
text( tmp, apply(mx,2,sum)+adj, label=apply(mx,2,sum), cex=.6,xpd=NA)
text( par('usr')[1]-(par('usr')[2]-par('usr')[1])/8, par('usr')[4]*1.15, label= round( back, 2), cex=1, xpd=NA, col='darkgrey')
dev.off()
}