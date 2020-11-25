#' [Fig 2A]  
#' Bubble plot of representative cell/tissue types.


library(pheatmap)
library(dplyr)
source('functions.r')
load('ge2an.rdata')
load('tp_mn.rdata')
load('tp_fr.rdata')
load('fig2a_bb.rdata')


node=type.gene$V1
gene2=tapply(type.gene$V2,INDEX = type.gene$V1,FUN = function(x){unlist(strsplit(x,split = '+',fixed = T))})
gene2=gene2[node]
gene2=unlist(get_id(unname(unlist(gene2))))


# Get average expression and avg fr by bubble_type 
type.avg.bb=avgExp(mx = tp_mn,gene = gene2,clu.list = na.omit(type.ann))
type.fr.bb=avgExp(mx = tp_fr,gene = gene2,clu.list = na.omit(type.ann))
type.avg.bb=type.avg.bb[gene2,node];rownames(type.avg.bb)=ge2an[rownames(type.avg.bb)]
type.fr.bb=type.fr.bb[gene2,node];rownames(type.fr.bb)=ge2an[rownames(type.fr.bb)]
type.avg.bb=apply(type.avg.bb,2,function(x)ifelse(x>5,5,x))
type.fr.bb=apply(type.fr.bb,2,function(x)ifelse(x>0.75,0.75,x))


# Bubble plot
plot_bubble(pct.mx = type.fr.bb,scale = 2.5,avg.mx = type.avg.bb,xlab.srt = 90,fontsize_row = 0.7,color = 'red',xlab_xadj = 0.9)



#' [Fig 2D]  
#' Spatial composition

load('clu_ann.rdata')
load('clu_orig.rdata')

# Calculate cell percentage with unexpected section origin
p=sapply(rownames(clu_orig), function(x) { 
  orig.r=get_body_part(clu_1012[[x]]) # real origin
  part=colnames(clu_orig)[3:6] # section names
  1-( sum(orig.r%in%part[clu_orig[x,3:6]==1])/sum(orig.r%in%part) ) } )

cluNew=clu_1012[rownames(clu_orig)]
body_distr=cbind(head=sapply(cluNew,function(x) unname(table(get_body_part(x))['head'])),
                 ht=sapply(cluNew,function(x) unname(table(get_body_part(x))['ht'])),
                 trunk=sapply(cluNew,function(x) unname(table(get_body_part(x))['trunk'])),
                 tv=sapply(cluNew,function(x) unname(table(get_body_part(x))['tv'])),
                 vi=sapply(cluNew,function(x) unname(table(get_body_part(x))['vi'])),
                 lv=sapply(cluNew,function(x) unname(table(get_body_part(x))['lv'])),
                 limb=sapply(cluNew,function(x) unname(table(get_body_part(x))['limb'])),
                 all=sapply(cluNew,function(x) sum(get_body_part(x)%in%c('head','trunk','vi','limb'))))
rownames(body_distr)=names(cluNew)
body_distr=apply(body_distr,2,function(x)ifelse(is.na(x)==T,0,x)) %>% as.data.frame()


# Only use samples with single origin
body_distr2=body_distr[,c('head','trunk','limb','vi','all')]
for(i in 1:4){body_distr2[,i]=body_distr2[,i]/body_distr2[,5]};body_distr2=body_distr2[,1:4]

# cumsum for plot
body_distr3=body_distr2;for(i in 1:3){body_distr3[,i]=rowSums(body_distr2[,i:4])}
body_distr3=body_distr3[complete.cases(body_distr3),]


chunk=clu_orig[rownames(clu_orig)%in%rownames(body_distr3),'Block']
p2=p[names(p)%in%rownames(body_distr3)]

t1=c('neural progenitor','neuron','epidermis')
t2=c('sensory neuron','schwann','craniofacial','head mesoderm','somite','IM')
t3=c('somatic LPM','limb','splanchnic LPM','PGC')
t4=c('endothelium','blood','endoderm','epithelium','fibroblast','miscellaneous')

ta=list(t1,t2,t3,t4)
na=c(0,sum(chunk%in%c(t1)),sum(chunk%in%c(t1,t2)),sum(chunk%in%c(t1,t2,t3)))
par(mfrow=c(4,1),mai=c(0.1,1,0.1,0.5))
invisible(lapply(1:4, function(i){
  t=ta[[i]];n=na[i]
  barplot(body_distr3$head[chunk%in%t],col='#8FAADC',border = '#8FAADC',width = 1,xlim = c(0,sum(chunk%in%t)),space = 0,ylim=c(0,1.1))
  barplot(body_distr3$trunk[chunk%in%t],col='#D0CECE',border = '#D0CECE',width = 1,space = 0,add=T)
  barplot(body_distr3$limb[chunk%in%t],col='#A9D18E',border = '#A9D18E',width = 1,space = 0,add=T)
  barplot(body_distr3$vi[chunk%in%t],col='#F4B183',border = '#F4B183',width = 1,space = 0,add=T)
  clip(0, sum(chunk%in%t), 0, 1)
  abline(v=c(0,cumsum(table(chunk)[t])),h=c(0,1))
  lines(x=((1:324)[chunk%in%t]-0.5-n),y=p2[chunk%in%t],col='red',lwd=2)
}))



