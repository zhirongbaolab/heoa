
#' [Fig.3A] Signaling center  
#' We examined the signaling capacity of the 9 known signaling centers identified in our data, namely, ZPA and AER in limb patterning, and 7 in neural tube and the brain. Signaling moleculars were download from Ramilowski 2015.

library(dplyr,quietly = T)
library(matrixStats)
source('functions.r')
load('ge2an.rdata')
load('tp_mn.rdata')
load('tp_fr.rdata')


#' The Type ID and Type anno of the 9 signaling centers
sc.tp=c('neu1-5-1'='ANR','neu2-1-30'='MHB','neu1-3-2'='ZLI','neu2-1-29'='FP.brain','neu2-1-13'='FP.spinal',
        'neu2-1-16'='RP.brain','neu2-1-1'='RP.spinal','epi-21'='AER','forelimb-2-7'='ZPA')  
sc.known.gene=c('FGF8','FGF17','FGF19','SHH','BMP4','FGF9','WNT3A','WNT4')


#' Filter ligands using 2 criterias:  
#' 1) robust expression in at least 1 signaling center (avg>=0.5 & fraction>=40%)  
#' 2) not broadly expressed (have at least 10% of all 333 cell types with <0.2 average expression)  

all_lig=read.table("all-ligands.txt",stringsAsFactors = F)$V1 %>% get_id
all_lig=intersect(all_lig,rownames(tp_mn))
mn=tp_mn[all_lig,names(sc.tp)]
fr=tp_fr[all_lig,names(sc.tp)]

g1=apply(as.matrix(mn>=0.5)+as.matrix(fr>=0.4),1,function(x){max(x)==2}) #1) 87genes
g2=apply(as.matrix(tp_mn[all_lig,]<0.2),1,function(x) {sum(x)>=0.1*ncol(tp_mn)}) #2) 663genes

mx=tp_mn[all_lig[g1&g2],names(sc.tp)] #1&2) 44 genes
mx=log2(1+mx);mx=apply(mx,2,function(x){ifelse(x>2,2,x)} )
rownames(mx)=ge2an[rownames(mx)];colnames(mx)=sc.tp[colnames(mx)]


#' Order genes according to 'known' or 'novel'
library(pheatmap)
p=pheatmap(mx,cluster_rows = T,cluster_cols = F,
           color =  colorRampPalette(c('blue','white','red2'))(100),
           fontsize_row = 7,show_colnames = T,border_color = F)
dev.off()
gene=rownames(mx)[p$tree_row$order]
gene=c(sc.known.gene,setdiff(gene,sc.known.gene))

#' Plot Fig3a
#pdf('fig3a.pdf',width = 3,height = 4)
pheatmap(mx[gene,],cluster_rows = F,cluster_cols = F,
         color =  colorRampPalette(c('blue','white','red2'))(100),fontsize=6,
         show_colnames = T,border_color = F,gaps_row = length(sc.known.gene))
#dev.off()
rm(g1,g2,mn,fr,p)



#' [Fig.3B]  
#' Expression of ligands, receptors, and antagonists of 5 major pathways (Fgf, Notch, Shh, TGF-beta, Wnt)  
#' Cell types are ordered by blocks (developmental systems)  
#' Only genes with robust expression were plotted   

load('tp_zs.rdata')
load('clu_ann.rdata')

#' Select genes with robust expression in at least 1 cell type (avg>=0.5 & fraction>=0.4) 
sig.all=read.table('signaling-gene.txt',header = T,row.names = 1,stringsAsFactors=F)
sig.gene=get_id(rownames(sig.all))

#' Generate expressed gene list for each type (avg>=0.5,fr>=40%)
be1=apply(tp_mn[sig.gene,],2,function(x){sig.gene[x>=0.5]}) # a list of genes with avg>=0.5 for each type
be2=apply(tp_fr[sig.gene,],2,function(x){sig.gene[x>=0.4]}) # a list of genes with fr>0.4 for each type
type_be=sapply(colnames(tp_mn),function(x){intersect(be1[[x]],be2[[x]])}) # each type have a list of genes which are considered as robust expression
length(unique(unlist(type_be))) # total 70 genes considered expressed in at least 1 cell type


#' Plot Fig3b heatmap using expressed genes. Genes/Rows are ordered manually. Cell types/Cols are ordered by developmental systems 
sig.plot=read.table('signaling-gene-plot.txt',row.names = 1,header = T,stringsAsFactors = F,sep='\t') # read ordered gene list for plot
gene=get_id(rownames(sig.plot))
type=rownames(clu_ann_1012)

mx=tp_zs[gene,type]; rownames(mx)=ge2an[rownames(mx)]
mx=apply(mx,2,function(x) ifelse(x>4,4,x))
mx=apply(mx,2,function(x) ifelse(x< -1,-1,x))


#' To better separate between blocks, add one blank column between blocks  
#' Gap rows to distinguish pathways  
gaps_row=cumsum(table(sig.plot$Pathway)[unique(sig.plot$Pathway)]) # gap between different pathway
gaps_col=cumsum(table(clu_ann_1012$Block)[unique(clu_ann_1012$Block)]) # gap between different block
gaps_col2=sapply(1:18, function(x) gaps_col[x]+x)
#c=paste(c(1,gaps_col[1:18]),gaps_col[1:19],sep=':') %>% paste(collapse = ',')
mx2=mx[,c(1:55,55:88,88:90,90:92,92:99,99:131,131:152,152:168,168:171,171:186,186:229,229:265,265:271,271:284,284:293,293:294,294:309,309:327,327:333)]
mx2[,gaps_col2]=NA; colnames(mx2)[gaps_col2]=NA

#' Plot heatmap
load('colours.rdata')
ann_col=list(Block=col.block)
pheatmap(mx2,cluster_cols = F,cluster_rows = F,show_colnames = F,
           gaps_row = gaps_row,color = colorRampPalette(c('blue','white','red'))(100),na_col = 'black',
           fontsize_row = 6,annotation_col =clu_ann_1012[,"Block",drop=F],annotation_colors = ann_col)


#'  
#' [Fig.S3A]  
#' We then calculated the percentage of cell types that show ligand, receptors, and antagonists expression.  
#' A cell type is considered with ligand(receptor, antagonist) expression of a particular pathway 
#' if it expresses any of the ligands(receptors, antagonists) of that pathway.
redun.type=tapply(sig.gene,paste(sig.all$pathway,sig.all$type,sep='-'),function(x)x) # group genes according to pathways and types
redun.mx=sapply(type_be,function(expr){ sapply(redun.type,function(t){sum(expr%in%t)} ) })
redun.mx=redun.mx[,colnames(tp_mn)] # reorder cols


#' Compute %cell-types express ligands (belonging to these 5 pathways)
lig <- redun.mx[grep('ligand',rownames(redun.mx)),] %>% colSums() 
rec <- redun.mx[grep('receptor',rownames(redun.mx)),] %>% colSums()
ant <- redun.mx[grep('antagonist',rownames(redun.mx)),] %>% colSums()

sum(lig>0)/length(lig) # 60.4% express at least 1 ligand
sum(rec>0)/length(lig) # 57.1% express at least 1 receptor
sum(ant>0)/length(ant) # 77.8% express at least 1 antagonist

#' Plot pie charts which show the expression fraction
par(mfrow=c(1,3),mai=c(.1,.1,.1,.1)) # 400x250
pie(c(sum(lig>0)/length(lig),1-sum(lig>0)/length(lig)),col=c(col_nrc[10],"white"),main='ligand',labels = "",radius = 1)
pie(c(sum(rec>0)/length(rec),1-sum(rec>0)/length(rec)),col=c(col_nrc[10],"white"),main='receptor',labels = "",radius = 1)
pie(c(sum(ant>0)/length(ant),1-sum(ant>0)/length(ant)),col=c(col_nrc[10],"white"),main='antagonist',labels = "",radius = 1)

#' Furthermore, we explored the expression fraction of lig/rec/ant in each developmental system
lig.b=tapply(lig,clu_ann_1012[names(lig),"Block"],function(x){c(sum(x>0)/length(x),sum(x==0)/length(x))}) %>% sapply(function(y){y}) 
lig.b=lig.b[,rev(unique(clu_ann_1012$Block))]
rec.b=tapply(rec,clu_ann_1012[names(rec),"Block"],function(x){c(sum(x>0)/length(x),sum(x==0)/length(x))}) %>% sapply(function(y){y}) 
rec.b=rec.b[,rev(unique(clu_ann_1012$Block))]
ant.b=tapply(ant,clu_ann_1012[names(ant),"Block"],function(x){c(sum(x>0)/length(x),sum(x==0)/length(x))}) %>% sapply(function(y){y}) 
ant.b=ant.b[,rev(unique(clu_ann_1012$Block))]

par(mfrow=c(1,4),mai=c(.3,.1,.1,.1)) # 300x300
barplot(lig.b,horiz = T,col = c(col_nrc[1],'white'),yaxt="n",xaxt="n"); axis(1, at=c(0,0.5,1),labels = c(0,0.5,1))
barplot(rec.b,horiz = T,col = c(col_nrc[2],'white'),yaxt="n",xaxt="n"); axis(1, at=c(0,0.5,1),labels = c(0,0.5,1))
barplot(ant.b,horiz = T,col = c(col_nrc[4],'white'),yaxt="n",xaxt="n"); axis(1, at=c(0,0.5,1),labels = c(0,0.5,1))
barplot(table(clu_ann_1012$Block)[rev(unique(clu_ann_1012$Block))],horiz = T,col=col_nrc[10],yaxt="n")   

#' Head mesoderm show unexpected complexity of signaling expression  
sum(clu_ann_1012$Block=='head mesoderm') # there are totally 21 cell types
hm.lig=redun.mx[grep('ligand',rownames(redun.mx)),rownames(clu_ann_1012)[clu_ann_1012$Block=='head mesoderm']]
sum(colSums(hm.lig)>0) #14 out of 21 show ligand expression
hm.rec=redun.mx[grep('receptor',rownames(redun.mx)),rownames(clu_ann_1012)[clu_ann_1012$Block=='head mesoderm']]
sum(colSums(hm.rec)>0) #11 out of 21 show receptor expression
hm.ant=redun.mx[grep('antagonist',rownames(redun.mx)),rownames(clu_ann_1012)[clu_ann_1012$Block=='head mesoderm']]
sum(colSums(hm.ant)>0) # 20 out of 21 show antagonist expression

#'  
#' [Fig.3C] 
#' We checked the expression breadth of ligands and receptors of different pathways  

t.pct=apply(redun.mx,1,function(x) {sum(x>0)/length(x)})
x=t.pct[grep('antagonist',names(t.pct),invert = T)]
par(mfrow=c(1,1),mai=c(.5,1,.5,.5))
barplot(cbind(x[grep('FGF',names(x))],x[grep('WNT',names(x))],x[grep('SHH',names(x))],
              x[grep('NOTCH',names(x))],x[grep('TGF',names(x))]),
        col = col_nrc[1:2],beside = T,ylim=c(0,0.49),ylab="% of cell-type expressed",
        names.arg = c('Fgf','Wnt','Shh','Notch','TGF-β'))
legend('topright',fill = col_nrc[1:2],bty = "n",legend = c("ligand",'receptor'))


#'  
#' [Fig.3D] 
#' We infered potential pathway crosstalks according to the co-expression of receptors.  
#' [Fig.3D letf panel] - cell types with crosstalks 

# Only select rows of receptors and assign to binary (0 or 1)
rec.mx=redun.mx[grep('receptor',rownames(redun.mx)),]
rec.mx=apply(rec.mx,2,function(x) ifelse(x>=1,1,0))

# types with 2 or more receptors per block ( number of types and percentage of types)
rec.sum= colSums(rec.mx)
pct=tapply( rec.sum, clu_ann_1012[colnames(redun.mx),"Block"], function(x) sum(x>=2)/length(x) )
num=tapply( rec.sum, clu_ann_1012[colnames(redun.mx),"Block"], function(x) sum(x>=2) )

# Barplot showing most prodominent cell types
par(mfrow=c(1,2),mai=c(.5,.15,.1,.1))
barplot(rev(num[unique(clu_ann_1012$Block)]),horiz = T,yaxt="n",col = col_nrc[3],xaxt='n',xlim=c(0,30))
axis(1,at=c(0,15,30),cex=0.4,labels = c(0,"",30))
barplot(rev(pct[unique(clu_ann_1012$Block)]),horiz = T,yaxt="n",col="grey40",xaxt='n',xlim=c(0,0.5))
axis(1,at=c(0,0.25,0.5),labels = c(0,"",0.5),cex=0.4)


#'  
#' [Fig.S3b]
#' Crosstalk distribution of all cell tyeps  

library(chorddiag)

# prepare linking matrix 
m=matrix(NA,5,5)
for(i in 1:5) {
  for(j in 1:5) {
    m[i,j]=sum(colSums(rec.mx[c(i,j),])==2) #rownames(clu_ann_1012)[clu_ann_1012$Block=='limb']
  }
}
for(i in 1:5){ m[i,i]=0} # remove self-linking

# Build the chord diagram
names=c('FGF','NOTCH','SHH','TGF','WNT')
dimnames(m) <- list(row = names,col = names)
chorddiag::chorddiag(m, groupColors = col_nrc2[5:9],groupnamePadding = 25,groupnameFontsize = 12,tickInterval = 1,groupThickness = 0.2,showGroupnames = T)


#' Whether the frequency of these crosstalks are expected? We calculated the probability.

num_tp=apply(rec.mx,1,function(x){sum(x>0)}) # for each pathway, the number of expressed cell types.
pct_tp=apply(rec.mx,1,function(x){sum(x>0)/length(x)}) # for each pathway, the frequency of receptorexpression

# all possible combinations of 5 pathways and the number/frequency of co-expressed cell types
C5.2=list(c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5))
C5.2.num=sapply(C5.2, function(i){ sum(colSums(rec.mx[i,])==2)})

out=sapply(1:10,function(i){
  a=pct_tp[ C5.2[[i]][1] ] # pathway 1 pct
  a.num=num_tp[ C5.2[[i]][1] ] # pathway 1 number
  b=pct_tp[ C5.2[[i]][2] ] # pathway 2 pct
  b.num=num_tp[ C5.2[[i]][2] ] # pathway 2 number
  axb=pct_tp[ C5.2[[i]][1] ] * pct_tp[ C5.2[[i]][2] ] # a-b co-expression expected percentage
  axb.real=C5.2.num[i] # a-b co-expression real number
  p=unname(choose(333,axb.real) * (axb)^axb.real * (1-axb)^(333-axb.real)) # probability of co-expression in current number of cell types
  unname(c(names[C5.2[[i]][1]],names[C5.2[[i]][2]],a,a.num,b,b.num,axb,axb.real,p))
})
rownames(out)=c('pathway1','pathway2','1.pct','1.num','2.pct','2.num','1x2.pct','1x2.real','p-v')
out=t(out)
#' The frequency of Fgf-Notch; Fgf-Shh; Fgf-Wnt; Notch-Tgf; Shh-Wnt co-expression are significantly higher than expected.  

#' [Fig.3D right panel] - co-expression detail of 3 specific developmental systems  

block='neural progenitor'
m=matrix(NA,5,5)
for(i in 1:5) {
  for(j in 1:5) {
    m[i,j]=sum(colSums(rec.mx[c(i,j),rownames(clu_ann_1012)[clu_ann_1012$Block==block]])==2) #
  }
}
for(i in 1:5){ m[i,i]=0} # remove self-linking

dimnames(m) <- list(row = names,col = names)
chorddiag::chorddiag(m, groupColors = col_nrc2[5:9],groupnamePadding = 20,groupnameFontsize = 12,tickInterval = 1,groupThickness = 0.2,showGroupnames = T)

block='head mesoderm'
m=matrix(NA,5,5)
for(i in 1:5) {
  for(j in 1:5) {
    m[i,j]=sum(colSums(rec.mx[c(i,j),rownames(clu_ann_1012)[clu_ann_1012$Block==block]])==2) #
  }
}
for(i in 1:5){ m[i,i]=0} # remove self-linking

dimnames(m) <- list(row = names,col = names)
chorddiag::chorddiag(m, groupColors = col_nrc2[5:9],groupnamePadding = 20,groupnameFontsize = 12,tickInterval = 1,groupThickness = 0.2,showGroupnames = T)

block='limb'
m=matrix(NA,5,5)
for(i in 1:5) {
  for(j in 1:5) {
    m[i,j]=sum(colSums(rec.mx[c(i,j),rownames(clu_ann_1012)[clu_ann_1012$Block==block]])==2) #
  }
}
for(i in 1:5){ m[i,i]=0} # remove self-linking

dimnames(m) <- list(row = names,col = names)
chorddiag::chorddiag(m, groupColors = col_nrc2[5:9],groupnamePadding = 20,groupnameFontsize = 12,tickInterval = 1,groupThickness = 0.2,showGroupnames = T)
