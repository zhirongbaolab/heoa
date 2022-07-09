library(dplyr,quietly = T)
library(monocle,quietly = T)
library(colorspace)
library(pheatmap)
library(matrixStats,quietly = T)
library(Matrix)


# matrix, cells, and genes
load('neu2_3.rdata') #normed matrix
load('clu_ann.rdata') #annotation
load('ge2an.rdata') #gene annotation
load('mm.ge2an.rdata')
load('c333_deg.rdata') #type_deg
load('colours.rdata')
load('basic.rdata')
source('functions.r')


#' [Fig 4A]   
#' Generate pseudo-AP axis using neural progenitor cells from hindbrain and spinal cord.  
#' All expressed HOX genes and 5 housekeeping genes (RPL5, ELAVL1, ATP2C1, ARPC2, ARPC1A, to prevent omitting cells with zero HOX gene expression) were used to reconstruct the pseudo-axis using Monocle 2.  

clu=clu_1012[paste('neu2-1-',c(1:13,16:29),sep='')]
cell=unlist(clu);length(cell)
hox_gene=grep("^HOX",ge2an,value = T)
gene=c(names(hox_gene),get_id(c("RPL5","ELAVL1","ATP2C1","ARPC2","ARPC1A",'CNPY1')))

mx=exp_mx[gene,cell];dim(mx)
mx=mx[, colSums(as.matrix(mx))>1];dim(mx)
mx=apply(mx,2,function(x)return(ifelse(x>0.3,x,0))) 
mono=mono_pseudo(mx = log2(1+mx),cell = colnames(mx),is_order = T,gene = rownames(mx))
rm(mx)

x=exp_mx[,cell] # for heatmap plot and for co-AP gene compute
x=x[,order(mono[[2]]$Pseudotime)];rownames(x)=ge2an[rownames(x)]

#prepare annotation info and color
x.body=get_body_part(colnames(x))
x.info=cbind(PT=mono[[2]]$Pseudotime[order(mono[[2]]$Pseudotime)], 
             ident=get_clu_ind(clu)[colnames(x)], 
             tv=sapply(x.body, function(x)ifelse(x=='tv','tv','other')), 
             trunk=sapply(x.body, function(x)ifelse(x=='trunk','trunk','other')),  
             ht=sapply(x.body, function(x)ifelse(x=='ht','ht','other')),
             head=sapply(x.body, function(x)ifelse(x=='head','head','other')) ) %>% as.data.frame
x.info$PT=as.numeric(x.info$PT)

pt.col=mapply(function(x){return(colorRampPalette(c('white','pink','red'))(100)[x])},
              x=as.numeric(cut(x.info$PT, breaks=100)))
ident.col=c(rainbow_hcl(15,150,55)[c(1:4,6,8:15)],rainbow_hcl(15,150,55)[c(1:4,6:15)]);names(ident.col)=names(clu)
ann_colors2=list(PT=pt.col,ident=ident.col,
                 head=c(head="#82CA3F",other='white'),trunk=c(trunk="red",other='white'),
                 ht=c(ht="blue",other='white'),tv=c(tv="turquoise",other='white'))

gene.p1=c('HOXA2','HOXB2','HOXA3','HOXB4','HOXB5','HOXB6',"HOXC9",'HOXA10','HOXC10')

otx_cell=cell[exp_mx['ENSG00000165588',cell]>1.5] #rm cells from mid-hindbrain 
x.s=x[,ncol(x):1];x.s=x.s[,setdiff(colnames(x.s),otx_cell)] #reverse cell order
x.s=apply(x.s[c(gene.p1),],2,function(x)return(ifelse(x>2,x,0))) #set low-value to 0
x.s=get_smooth(x.s,bin = ceiling(length(cell)/10),onlyLeft = T) #smooth expression
x.s=apply(log2(1+x.s),2,function(x) ifelse(x>0.05,x,0))

pheatmap(x.s[c(gene.p1),],cluster_rows = F,cluster_cols = F,annotation_col = x.info[ colnames(x.s),],
         col=colorRampPalette(c("white","blue",'navy'))(100),border_color = F,show_colnames = F,
         fontsize_row = 6,gaps_row = c(1:9),annotation_colors = ann_colors2,annotation_legend = F) 
rm(x.s)

#' [Fig 4B]  
#' We then defined a group of AP-related genes by calculating the correlation between gene expression and cell ordering on the trajectory.    

gene=c(setdiff( ge2an[unlist(c333_type_deg)], grep('^HOX',ge2an,value = T) ),'RP11-357H14.17')
y=x[gene,setdiff(colnames(x),otx_cell)] #remove cells from hind- and mid-brain
y=y[rowMaxs(as.matrix(y))>2,];y=log1p(y) #choose expressing genes
y.pt=mono[[2]][colnames(y),'Pseudotime']

# For each embryo, calculate the correlation value. Then, select genes that show notable correlation in at least two embryos
t=invisible(lapply(c('E1107','E0710','E0809','E1205'),function(i) {
  yE=get_embryo(colnames(y))==i;
  cor.E=apply(y[,yE],1,function(x) {cor(x,y.pt[yE])});
  cor.E=cor.E[complete.cases(cor.E)]
  names(cor.E)[cor.E < mean(cor.E)-4.5*sd(cor.E)| cor.E > mean(cor.E)+3.5*sd(cor.E)] 
}))
g=names(table(unlist(t)))[table(unlist(t))>=2]

# manually order genes 
g=c('CNPY1','CNTNAP2','LIX1','CYP26C1','LGI1','RP11-834C11.4','HOTAIRM1','ZADH2','SKAP2',
    'RP11-834C11.6','FLJ12825','RAB38','SCUBE2','CALCB','SPARCL1','RP11-357H14.17')
mx= y[g,] %>% get_smooth(bin = ceiling(ncol(y)/30))
zs=t( apply( mx,1, function(x){ (x-mean(x))/sd(x) }) )
zs=apply(zs, 2,function(x) ifelse(x>2,2,x) )
zs=rbind(zs,rep(-2,ncol(zs)))
pheatmap(zs[,ncol(zs):1],cluster_cols = F,show_colnames = F,cluster_rows = F,gaps_row = 16,
         color = colorRampPalette(c('blue','white','red'))(100))
rm(zs,t)

#' 
#' [Fig S4A]  
#' Temporal expression pattern of the five lncRNAs and their upstream hox genes.   

cell.p=rev(setdiff(colnames(x),otx_cell))

# Average expression of hox 1-13
hox.mx=sapply(1:13,function(i){ x[grep(paste("^HOX[A-D]",i,"$",sep=""),ge2an,value = T),cell.p,drop=F] %>% as.matrix %>% colMeans   }) %>% t
rownames(hox.mx)=paste("HOX",1:13,sep="")

gene=c('RP11-834C11.4','HOTAIRM1','RP11-834C11.6','FLJ12825','RP11-357H14.17')
y= x[gene,cell.p] %>% rbind(hox.mx) %>% get_smooth %>% get_scale %>% get_smooth

col.p=c(col6[c(2:4,6:7)],colorRampPalette(c('grey90','black'))(13))

plot(NA,xlim=c(1,length(cell)),ylim=c(0,1),xlab='',ylab='',xaxt='n',cex.axis=0.8)
abline(v=1320,lty=2) 
r=c(7,9);invisible(lapply(r,function(i){lines(1:ncol(y),y[i,],col=col.p[i],lwd=2,lty=2)} ))
r=c(1:4);invisible(lapply(r,function(i){lines(1:ncol(y),y[i,],col=col.p[i],lwd=2)} ))
r=c(1:4,7,9);legend('topleft',legend = rownames(y)[r],lwd=3,col=col.p[r],bty='n',cex=0.7)

plot(NA,xlim=c(1,length(cell)),ylim=c(0,1),xlab='',ylab='',xaxt='n',cex.axis=0.8)
abline(v=1320,lty=2)
r=c(14,15);lapply(r,function(i){lines(1:ncol(y),y[i,],col=col.p[i],lwd=2,lty=2)} ) 
r=c(5);lapply(r,function(i){lines(1:ncol(y),y[i,],col=col.p[i],lwd=2)} )
r=c(5,14,15);legend('topleft',legend = rownames(y)[r],lwd=3,col=col.p[r],bty='n',cex=0.7) 


#' [Fig 4C]  
#' 
#' [Fig 4D]  
#' We compared our data with a recent dataset of E9.5-E13.5 mouse neural tube single-cell transcriptomes (Delile, J. et al. Development, 2019).  

# 1) use only first 3 stages of mouse data (E9.5, E10.5, E11.5)
# 2) Use only cells from spinal cord of our human data (dp1-dp6, p0-p3, pMN, rp,fp)

load('bubble_plot_genelist.rdata')

# 1) human progenitor-spinal (13col x bb_gene)
clu=clu_1012[paste('neu2-1-',1:13,sep='')]
names(clu)=c('RP','dp1','dp2','dp3','dp4','dp5','dp6','p0','p1','p2','pMN','p3','FP')
hg_prog_avg_mx=avgExp(mx = exp_mx,clu.list = clu,gene=get_id(bb_gene_prog))
hg_prog_pct_mx=pctExp(mx = exp_mx,clu.list = clu,gene=get_id(bb_gene_prog))
rownames(hg_prog_avg_mx)=ge2an[rownames(hg_prog_avg_mx)]
rownames(hg_prog_pct_mx)=ge2an[rownames(hg_prog_pct_mx)]

# 2) mouse progenitor-first3stages (13col x bb_gene)
load('mouse_prog_data.rdata')
cell=rownames(megadata.pro)[megadata.pro$timepoint%in%c('9.5','10.5','11.5')]
mx.prog.mm=apply(mx.progen[,cell],2,function(x){x/sum(x)*10000})
rm(mx.progen)

ident=megadata.pro[cell,'Type_step2'] 
names(ident)=colnames(mx.prog.mm)
clu.mm= tapply(names(ident),ident, function(x){x}) %>% as.list
type_name=c('RP','dp1','dp2','dp3','dp4','dp5','dp6','p0','p1','p2','pMN','p3','FP')
clu.mm=clu.mm[type_name]

mm_prog_avg_mx=avgExp(mx = mx.prog.mm,clu.list = clu.mm,gene = get_mm_id(bb_gene_prog))
mm_prog_pct_mx=pctExp(mx = mx.prog.mm,clu.list = clu.mm,gene = get_mm_id(bb_gene_prog))
rownames(mm_prog_avg_mx)=mm.ge2an[rownames(mm_prog_avg_mx)]
rownames(mm_prog_pct_mx)=mm.ge2an[rownames(mm_prog_pct_mx)]

#pdf('fig4d.pdf',width = 9.5,height = 4)
plot_pair_bubble_sideByside(pct.mx1 = hg_prog_pct_mx,pct.mx2 = mm_prog_pct_mx,
                            avg.mx1 = hg_prog_avg_mx,avg.mx2 = mm_prog_avg_mx,
                            scale1 = 2.8,scale2 = 2.8,sep.line = T,ylab_xadj = -0.5) 
#dev.off()
#rm(mx.prog.mm,megadata.pro)
rm(mm_prog_avg_mx,mm_prog_avg_mx,hg_prog_avg_mx,hg_prog_avg_mx)



#' [Fig S4B]  
#' Bubble plot of neurons

# 3) human neuron-spinal (12col x bb_gene)
clu=clu_1012[paste('neu3-1-',c(12:21,23),sep='')]
names(clu)=c('dl1','dl2','dl4','dl5','dl6','V0','V1','V2a','V2b','MN','V3')
hg_neu_avg_mx=avgExp(mx = exp_mx,clu.list = clu,gene=get_id(bb_gene_neu))
hg_neu_pct_mx=pctExp(mx = exp_mx,clu.list = clu,gene=get_id(bb_gene_neu))
rownames(hg_neu_avg_mx)=ge2an[rownames(hg_neu_avg_mx)]
rownames(hg_neu_pct_mx)=ge2an[rownames(hg_neu_pct_mx)]


# 4) mouse neuron-first3stages (12col x bb_gene)
load('mouse_neu_data.rdata')
mx.neuron.mm=mx.neuron[,megadata.neu$timepoint%in%c('9.5','10.5','11.5')] %>% apply(2,function(x){x/sum(x)*10000})
rm(mx.neuron)

ident=megadata.neu[megadata.neu$timepoint%in%c('9.5','10.5','11.5'),'Type_step2']
names(ident)=colnames(mx.neuron.mm)
clu.mm=tapply(names(ident),ident, function(x){x})
type_name=c('dl1','dl2','dl4','dl5','dl6','V0','V1','V2a','V2b','MN','V3')
clu.mm=clu.mm[type_name]

mm_neu_pct_mx=pctExp(mx = mx.neuron.mm,clu.list = clu.mm,gene = get_mm_id(bb_gene_neu))
mm_neu_avg_mx=avgExp(mx = mx.neuron.mm,clu.list = clu.mm,gene = get_mm_id(bb_gene_neu))
rownames(mm_neu_avg_mx)=mm.ge2an[rownames(mm_neu_avg_mx)]
rownames(mm_neu_pct_mx)=mm.ge2an[rownames(mm_neu_pct_mx)]

#pdf('figs4b.pdf',width = 9.5,height = 4)
plot_pair_bubble_sideByside(pct.mx1 = hg_neu_pct_mx,pct.mx2 = mm_neu_pct_mx,
                            avg.mx1 = hg_neu_avg_mx,avg.mx2 = mm_neu_avg_mx,
                            scale1 = 2.5,scale2 = 2.5,sep.line = T,ylab_xadj = -0.5) 

#dev.off()
rm(mx.neuron.mm,megadata.neu)


#' [Fig S4C]  
dp4=clu_1012[[rownames(clu_ann_1012)[grep('dp4',clu_ann_1012$type_anno)]]] 
dp1=clu_1012[[rownames(clu_ann_1012)[grep('dp1',clu_ann_1012$type_anno)]]]

mx=exp_mx[get_id(c('MSX1','ASCL1','MSX2')),colnames(exp_mx)%in%c(dp1,dp4)]
mx=t(as.matrix(log2(1+mx)))%>% as.data.frame()
colnames(mx)=c('MSX1','ASCL1','MSX2')

library(patchwork)
p1=ggplot(data = mx[dp1,],aes(x=MSX1,y=ASCL1,size=MSX2))+geom_point(colour=alpha('blue',0.5))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab("MSX1")+ylab("ASCL1")+xlim(0,4)+ylim(0,4)
p2=ggplot(data = mx[dp4,],aes(x=MSX1,y=ASCL1,size=MSX2))+geom_point(colour=alpha('red',0.5),fg='black')+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab("MSX1")+ylab("ASCL1")+xlim(0,4)+ylim(0,4)
p1+p2

# mouse


dp4.mm=rownames(megadata.pro)[megadata.pro$timepoint%in%c(9.5,10.5,11.5)&megadata.pro$Type_step2_unique=='dp4'] # 231
dp1.mm=rownames(megadata.pro)[megadata.pro$timepoint%in%c(9.5,10.5,11.5)&megadata.pro$Type_step2_unique=='dp1'] #174

mx.mm=mx.prog.mm[get_mm_id(c('MSX1','ASCL1','MSX2')),c(dp1.mm,dp4.mm)]
mx.mm=t(as.matrix(log2(1+mx.mm)))%>% as.data.frame()
colnames(mx.mm)=c('MSX1','ASCL1','MSX2')


p1=ggplot(data = mx.mm[dp1.mm,],aes(x=MSX1,y=ASCL1,size=MSX2))+
  geom_point(colour=alpha('blue',0.5))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("MSX1")+ylab("ASCL1")+xlim(0,4)+ylim(0,4)+
  scale_size_continuous(breaks=c(0,1,2,3,4),limits = c(0,3))

p2=ggplot(data = mx.mm[dp4.mm,],aes(x=MSX1,y=ASCL1,size=MSX2))+
  geom_point(colour=alpha('red',0.5))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("MSX1")+ylab("ASCL1")+xlim(0,4)+ylim(0,4)+
  scale_size_continuous(breaks = c(0,1,2,3,4),limits = c(0,3))
p1+p2



