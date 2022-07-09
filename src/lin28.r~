####### LIN28 expression delineates developmental stages in vertebrate embryogenesis

####################
## Systemically changing genes in human scRNA-seq
####################

# load human scRNA-seq data
load(file=paste('data/norm_mx_c333.rdata',sep='')) # normalized matrix
load('data/tp_mn_fr.rdata') # From 'clustering.r', 'tp_mn': the mean expression of each gene in each cell type; 'tp_fr': the fraction of detected cell (>0) of each gene in each cell type.
all_deg<-scan('data/all_deg.txt',what=character(0)) # From 'clustering.r', all DEGs
cell_anno<-read.table('data/cell_annotate.csv',sep=',',as.is=T,header=T) 
all_clu<- tapply( cell_anno[,1], cell_anno[,2], function(x){x})  # cell list of each cell type
clu_anno<- do.call( rbind, lapply(names(all_clu), function(x){ as.character(cell_anno[ cell_anno[,2]==x,][1,4:3]) })) # Annotation of cell types
rownames(clu_anno)<- names(all_clu)


# Load functions
source('function.r') 
source('lin28_function.r') 

# To avoid distortion by cell types with low number of cells, cell types with low number of cells were merged with similar cell type by a greedy approach.
# Define similar cell types by correlation
tp_deg2_corpe<- apply( tp_mn[all_deg, ], 2, function(x){
  apply( tp_mn[all_deg, ], 2, function(y){ cor(x,y) })
})
# set cutoff on correlation for similar cell types and only consider cell types within the same developmental system
mg_chunk<- clu_anno[,2]
names(mg_chunk)<- rownames(clu_anno)
mg_chunk[mg_chunk=='epidermis']<-'epithelium'
tp_cls<-sapply( rownames(tp_deg2_corpe), function(x, cut=.7){
  res<- colnames(tp_deg2_corpe)[ mg_chunk==mg_chunk[x] & tp_deg2_corpe[x,]>=cut ]
  res<- res[ order( -tp_deg2_corpe[x,res] ) ]
  res<- setdiff( res, x)
  return(res)
}) # candidates for merging

# merge cell types
clu_emb<-t( sapply( all_clu, function(x){
  res<-rep(0,5)
  names(res)<-c('7','0','5','9','6')
  emb<-table(get_emb(x))
  res[ names(emb) ]<- emb
  return(res)
  }
))
cn1<-60;cn2<-40 # cell number cutoff, use this cutoff. 'cn1': CS12+CS13 > 60; 'cn2': at least two out of three > 40
tp_np<-rownames(clu_emb)[ ! apply( clu_emb, 1, function(x){ sum(x[1:2])>=cn1 & (sum(x[3:5]>=cn2)>=2) }) ] # cell types that do not pass the threshold of cell number 
length(tp_np<- setdiff( tp_np, rownames(clu_anno)[ clu_anno[,2] %in% 'limb'] )) # limb as special case that will be merged by stage
test_num<-function(x){ sum(x[1:2])>=cn1 & (sum(x[3:5]>=cn2)>=2) }
s2_num<- function(x){ c( sum(x[1:2]), sort(x[3:5],decreasing=T)[1:2] )}
# start from begin to end: greedy search
length(pool<- tp_np)
tp_gs<-t(sapply( rownames(tp_deg2_corpe), function(x){
  if(! x %in% tp_np ) return(c(1,'', s2_num(clu_emb[x,]) ))
  if(! x %in% pool) return( c(0,'',s2_num(clu_emb[x,]) )) # merged by other types already
  can<- intersect( tp_cls[[x]], pool )
  if( length( can )==0 ) return(c(0,'', s2_num(clu_emb[x,]) ))
  mg<-''
  pass<-0
  for( i in 1:length(can)){
    if( test_num(apply(clu_emb[c(x,can[1:i]),],2,sum)) ) break
  }
  mg<- paste( can[1:i], collapse=',')
  pass<- ifelse( test_num(apply(clu_emb[c(x,can[1:i]),],2,sum)),1,0)
  if(pass==1) pool<<- setdiff( pool, c(x,can[1:i])) # global variable
  return( c(pass,mg, s2_num(apply(clu_emb[c(x,can[1:i]),],2,sum)) ) )
}))
length(pool) # [1] 
tp_gs<- tp_gs[ rownames(tp_gs)[ clu_anno[,2]!='limb'], ] # exclude limb types for now
# merge limb cell types by stage
tp_gs<-rbind( tp_gs, 
c(1,paste(rownames(clu_anno)[sapply(rownames(clu_anno),function(x){ strsplit(x,split='-')[[1]][1] })=='forelimb'][-1],collapse=','),0,0,0),
c(1,paste(rownames(clu_anno)[sapply(rownames(clu_anno),function(x){ strsplit(x,split='-')[[1]][1] })=='hindlimb'][-1],collapse=','),0,0,0) )
rownames(tp_gs)[289:290]<- c('forelimb-1-1','hindlimb-1-1')

# result of merged cell types
mg_clu<- lapply( rownames(tp_gs)[tp_gs[,1]==1], function(x){
  if(tp_gs[x,2]=='') tp<-x else tp<- c(x, strsplit(tp_gs[x,2],split=',')[[1]])
  return(unlist(all_clu[tp]))
})
names(mg_clu)<- rownames(tp_gs)[tp_gs[,1]==1]

# calulate temporal genes: CS12/13 vs. CS15 in merged cell types
# remove batch-effect genes
bt_gene<- unique(get_id(scan('list/batch.txt',what=character(0))))  # batch-effect genes (including sex-specific genes), merged from: pubmed ID 30096314; pubmed ID 31835037
in_gene<- setdiff( rownames(norm_mx_c333), c(bt_gene))
length(in_gene3<- in_gene[ apply( (tp_mn[in_gene,]>=0.5&tp_fr[in_gene,]>=.4),1,sum )>0 ])
clu_temp<- lapply( names(mg_clu),function(x,g=in_gene3){
  c1<- mg_clu[[x]][ get_emb(mg_clu[[x]]) %in% c('0','7') ]
  c2<- mg_clu[[x]][ get_emb(mg_clu[[x]]) %in% c('5','9','6') ]
  res<-cal_wilp2( mx=norm_mx_c333, clu1=c1, clu2=c2, genes=g,is_pv=T )
})
names(clu_temp)<-names(mg_clu)
clu_temp_fd<- sapply( clu_temp,function(x){
  x[ p.adjust(x[,2],method='BH')>0.0001 | abs(x[,1])< log2(2), 1]<-0
  return(x[,1])
})
length(tempg<- setdiff(rownames(clu_temp_fd)[apply(clu_temp_fd!=0, 1, sum)>0],lysis_gene))
clu_temp_fdo<- sapply( clu_temp,function(x){  return(x[,1])  })

# classify temporal genes into 3 categories: specific, systemic genes, other genes
# specific genes: expressed sites < 20%
length(tempg_spec<- tempg[ apply(abs(clu_temp_fdo[tempg,])>0,1,sum)/ncol(clu_temp_fdo) < 0.20 ])
# other genes 1): not enough sites (50%) with valid change (1.2 fold)
length(tmpg<- setdiff(tempg,c(tempg_spec) )) # 
length(tempg_oth1<- tmpg[ apply(abs(clu_temp_fdo[tmpg,])>log2(1.2),1,sum)/ncol(clu_temp_fdo) < 0.5 ]) # 
# other genes 2): inconsistent within developmental system (up-regulated in > 1/3 cell types & down-regulated in > 1/3 cell types)
length(tmpg<- setdiff(tempg,c(tempg_spec) )) 
tmpg_info<- apply( clu_temp_fdo[tmpg,], 1, function(y, fd_cut=log2(1)){
  stat<-tapply(y,mg_chunk2[colnames(clu_temp_fdo)],function(x, perc_cut=1/3){
    if( length(x)<6 ) return(F)
    return( sum(x>fd_cut)/length(x)> perc_cut & sum(x< -fd_cut)/length(x)> perc_cut )
  })
  return(sum(stat))
})
length(tempg_oth2<- names(tmpg_info)[tmpg_info>0] )
# systemically changing genes: >2 fold in >1 cell type; > 1.2 fold in > 50% cell types; not inconsistent within developmental system
length(tempg_glo<- setdiff( tempg, c(tempg_oth1,tempg_oth2, tempg_spec)))
# systemically changing genes: down-regulated and up-regulated
non_neu<-! mg_chunk[colnames(clu_temp_fdo)] %in% c('neural progenitor','neuron')
hs_glo<- list( tempg_glo[apply( clu_temp_fdo[tempg_glo,non_neu],1, function(x){ sum(x>0)>sum(x<0) } )], 
tempg_glo[apply( clu_temp_fdo[tempg_glo,non_neu],1, function(x){ sum(x>0)<sum(x<0) } )] )

# a control group: systemically expressed genes that are not changing on temporal
clu_temp_fd_ns<- sapply( clu_temp,function(x){
  x[ p.adjust(x[,2],method='BH')>0.0001 & abs(x[,1])< log2(2), 1]<-0
  return(x[,1])
})
length(ntg1<- rownames(clu_temp_fd_ns)[ apply(clu_temp_fd_ns!=0,1,sum)==0 ])
length(ntg2<- rownames(clu_temp_fdo)[apply( clu_temp_fdo, 1, function(x){ sum(x!=0)/length(x)>.5 })])
length(ntg<- intersect(ntg1,ntg2)) 


####################
## Systemically changing genes in mouse (by bulk and scRNA-seq)
####################
# Data: Boroviak et al. 2015, Pijuan-Sala et al. 2019, and Cao et al. 2019
# Stage 1 vs. stage 2: DEGs between E4.5~5.5 and E2.5~3.5, and systemically expressed in scRNA-seq (at E6.5)
m1_raw<-read.table('data/lin28_data/Boroviak_mouse_table_s1.txt',sep='\t',as.is=T,header=T) # It can be downloaded from https://www.sciencedirect.com/science/article/pii/S1534580715006589?via%3Dihub
m1_s12_mn<- t(apply( m1_raw[,4:18], 1, function(x){ tapply(x, c(rep('s1',6),rep('s2',9)) , mean) }))
m1_s12_diff<- cbind(  log2(m1_s12_mn[,1]+1)-log2(m1_s12_mn[,2]+1), (m1_s12_mn[,1]+m1_s12_mn[,2])/2 )  # early - late, max mean
# read Pijuan-Sala et al. 2019
#!# Because the size is too large, raw data of Pijuan-Sala is not included in this repertoire. It can be downloaded from https://github.com/MarioniLab/EmbryoTimecourse2018
library(Matrix)
raw<-readMM('data/lin28_data/Pijuan-Sala_mouse/atlas/raw_counts.mtx')
print('read raw matrix')
sf<- scan('data/lin28_data/Pijuan-Sala_mouse/atlas/sizefactors.tab',what=numeric(0))
m2_ganno<-read.table('data/lin28_data/Pijuan-Sala_mouse/atlas/genes.tsv',sep='\t',as.is=T,row.name=1)
canno<-read.table('data/lin28_data/Pijuan-Sala_mouse/atlas/meta.csv',sep=',',as.is=T,header=T)
print('read cell annotation')
st<-canno[,4]
norm_mx<- t(t(as.matrix(raw))/sf)
print('normalize matrix')
st_norm_mx<- by( t(norm_mx), st, function(x){t(x)})
print('split by stage')
# mean by cell type  by stage
tp<-canno[,'celltype']
tp[ is.na(tp) ]<-'Outlier'
st_tp<- tapply( tp, canno[,'stage'], function(x){x})
st_tp_mn<-mapply( function(x,y){
  res<-by( t(x), y, function(z){ apply( t(z), 1, mean) })
  return(res)
},x=st_norm_mx, y=st_tp )
print('finish calculate mean')
e65_mn<- st_tp_mn[[1]]
rownames(e65_mn)<- rownames(m2_ganno)
e65_mn<- e65_mn[, setdiff(colnames(e65_mn),'Outlier')]
cute<- 0.85
print(c(cute,median(apply( e65_mn>=cute,2,sum))) ) # use number of expressed genes to choose a cutoff for expression in different dataset
length(e65_glo<- rownames(e65_mn)[apply( e65_mn>=cute, 1, sum)>= ceiling(ncol(e65_mn)*.5)]) # 
sapply(m1_s12_glo<- sapply( list( rownames(m1_s12_diff)[m1_s12_diff[,1]>1& m1_s12_diff[,2]>=10],
rownames(m1_s12_diff)[m1_s12_diff[,1]< -1& m1_s12_diff[,2]>=10] ), function(x){ intersect(x, e65_glo) }),length) # 178 254
mm_up_glo<- m1_s12_glo

# Stage 3 vs. stage 4: DEGs between E11.5 vs. E9.5~10.5 from scRNA-seq. 
#!# Since the mouse scRNA-seq data in this time window is based on nuclear sequencing and the low detection level makes the identification less reliable, we filtered the identified human systemically temporal genes for those whose mouse ortholog have consistent direction of changes and are systemically expressed in mouse scRNA-seq.
#!# Because the size is too large, raw data of Cao et al. 2019 is not included in this repertoire. It can be downloaded from https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing
raw<-readMM('data/lin28_data/cao_mouse/gene_count.txt')
print('read raw matrix')
m3g<-read.table('data/lin28_data/cao_mouse/gene_annotation.csv',sep='\t',as.is=T,row.name=1)
cell_info<- read.table( 'data/lin28_data/cao_mouse/cell_annotate.csv', header=T, sep=',', as.is=T)
cell_st<- cell_info[,'development_stage']
cell_mt<- cell_info[,'Main_cell_type']
cell_sf<- cell_info[,'Size_Factor']
m3g<- unname(sapply( genes, function(x){ strsplit(x, split='.', fixed=T)[[1]][1] }))
dim(hs_orth<-read.table('data/ensembl100_hs_xt_dr_mm_orth.txt',sep='\t',as.is=T,header=T)) # orthlog from ensembl100
glo_mm<- sapply( unlist(hs_glo), function(x){
  intersect( unique(hs_orth[hs_orth[,1]==x,7]),m3g)
})
val<- raw[which(m3g %in% unlist(glo_mm)),]
glo_st_tp<-tapply( cell_info[,1], cell_info[,'Main_cell_type'], function(x){ # only compare E9.5/10.5 (pool) vs. E11.5
    ind1<- x[cell_st[x] %in% c('9.5','10.5')]
    ind2<- x[cell_st[x] %in% c('11.5') ]
    if(length(ind1)==0) mn1<-rep(0,nrow(val)) else mn1<- apply( t(t(val[,cell_info[,1] %in% ind1])/cell_sf[ind1]), 1, mean)
    if(length(ind2)==0) mn2<-rep(0,nrow(val)) else mn2<- apply( t(t(val[,cell_info[,1] %in% ind2])/cell_sf[ind2]), 1, mean)
    res<-cbind(mn1,mn2)
    print(Sys.time())
    return(res)
})
glo_mn<-lapply(1:length(gind), function(x){
  res<-t(sapply( glo_st_tp, function(y){ y[x,] }))
  return(res)
})
glo_lfd<-t(sapply(glo_mn, function(x){ 
  xx<- (x[,2]+x[,1])/2
  res<- log2(xx/x[,3])
  res[ xx==0 ]<- -1
  res[ x[,3]==0 ]<- 1
  res[ x[,3]==0 & xx==0 ]<- 0
  return(res)
}))
# summary whether those genes are changing globally
hs_glo_mm_info<- mapply( function(x,y,cut_lfd=log2(1.2),lfd=glo_lfd ){ 
  res<-t(sapply(x, function(z){
    if( length(z)==0 ) return( c(-1,-1) )
    if( length(z)>1 ){
      if(y==1) z<-z[which.max(apply(lfd[z,],1,mean))]
      else z<-z[which.min(apply(lfd[z,],1,mean))]
    }
    mn_ind<- -c(12,15)
    xx<-lfd[z, mn_ind]
    return( c( mean(xx[abs(xx)>=cut_lfd]), sum( abs(lfd[z,])>= cut_lfd)*100/ncol(lfd) ))
  }))
  return(res)
}, x= list( glo_mm[1:sum(sapply(hs_glo,length)[1:2])], glo_mm[-(1:sum(sapply(hs_glo,length)[1:2]))] ) , y=c(1,2) )
mm_dp_glo<- list( rownames(hs_glo_mm_info[[1]])[hs_glo_mm_info[[1]][,1]>=log2(1.2)&hs_glo_mm_info[[1]][,2]>=50], rownames(hs_glo_mm_info[[2]])[hs_glo_mm_info[[2]][,1]<= -log2(1.2)&hs_glo_mm_info[[2]][,2]>=50] )


####################
## Systemically changing genes in frog (by bulk and scRNA-seq)
####################
# Data: Tan et al. 2013 and Briggs et al. 2018
#!# The gene mapping reported by paper was from old frog genome that we found there is error on gene annotatino, so we re-map to frog genome in Ensembl100. See "remap_Tan_frog/" folder.
get_htseq<-function(x){
  sig<-read.table(file=paste('../htseq_output/',x,'.txt',sep=''),as.is=T,header=F,row.name=1)
  sig_nm<-sig[-((dim(sig)[1]-4):dim(sig)[1]),1]
  names(sig_nm)<-rownames(sig)[-((dim(sig)[1]-4):dim(sig)[1])]
  return(sig_nm)
}
sam_info<-read.table('remap_Tan_frog/sam_info.txt',sep='\t',as.is=T,row.name=1)
mx<-sapply(rownames(sam_info)[-25],get_htseq)
lib<- apply(mx,2,sum)/1000000
norm_mx<- t(t(mx)/lib)
xt_all_exp<-t(apply( norm_mx[ ,],1,function(x){ tapply( x, sam_info[colnames(norm_mx),1], mean) })) 
xt_all_exp<-xt_all_exp[,c('2cell','4cell','8cell','16cell',colnames(xt_all_exp)[-(1:4)][order(as.numeric(sapply(substring(colnames(xt_all_exp)[-(1:4)],6), function(x){ strsplit(x,split='-')[[1]][1] })))] )]

# systemically expressed genes in scRNA-seq
#!# Because the size is too large, raw data of Briggs et al. 2018 is not included in this repertoire. It can be downloaded from GSE113074
raw<-read.table('data/GSE113074_Corrected_combined.annotated_counts.tsv',sep='\t',row.name=1,as.is=T) # [1]  26559 136966, take ~0.5 hr to read
st<- raw[7,]
tp<- raw[8,]
st14_mn<- sapply( names(st_tp[['Stage_14']])[-1], function(x){
  res<-apply( data.matrix(raw[-(1:9),tp==x]),1,mean)
  print(x)
  return(res)
})
st16_mn<- sapply( names(st_tp[['Stage_16']])[-1], function(x){
  res<-apply( data.matrix(raw[-(1:9),tp==x]),1,mean)
  print(x)
  return(res)
})
cute<- 0.15 # use number of expressed genes to choose a cutoff for expression in different dataset
sapply( cute, function(x){ median(apply(st16_mn>=x,2,sum)) })
length(s14_in_tp<- names(st_tp[['Stage_14']])[ st_tp[['Stage_14']]> 40 ])
s14_in_tp<- setdiff(s14_in_tp, 'Outlier')
length(xts_s14_glo<- rownames(st14_mn)[apply( st14_mn[,s14_in_tp]>=cute, 1, sum)>= floor(length(s14_in_tp)*.5) ])
length(s16_in_tp<- names(st_tp[['Stage_16']])[ st_tp[['Stage_16']]> 40 ])
s16_in_tp<- setdiff(s16_in_tp, 'Outlier')
length(xts_s16_glo<- rownames(st16_mn)[apply( st16_mn[,s16_in_tp]>=cute, 1, sum)>= floor(length(s16_in_tp)*.5) ])

# Stage 1 vs. stage 2 (as in 4 stages defined in paper): DEGs between stages 13~14 and stages 6~8 (as in frog embryonic stages), and systemically expressed in scRNA-seq (stage 14)
up_cmp<- cal_wilp2( mx1=xt_all_exp[,6:7], mx2=xt_all_exp[,9:10],cut_mean=10,small_val=1)
sapply( xtb_up_diff<- list( names(up_cmp)[up_cmp> 1], names(up_cmp)[up_cmp< -1] ),length) # Down and Up
xtb_up_glo<-xtb_up_diff
length(xtb_up_glo[[1]]<-unique(hs_orth2[ hs_orth2[,3] %in% xtb_up_glo[[1]] & hs_orth2[,4] %in% xts_s14_glo, 3]))
length(xtb_up_glo[[2]]<-unique(hs_orth2[ hs_orth2[,3] %in% xtb_up_glo[[2]] & hs_orth2[,4] %in% xts_s14_glo, 3]))

# Stage 3 vs. stage 4 (as in 4 stages defined in paper): DEGs between stages 31-34 vs. stages 15-18 (as in frog embryonic stages), and systemically expressed in scRNA-seq (stage 16)
dp_cmp<- cal_wilp2( mx1=xt_all_exp[,11:12], mx2=xt_all_exp[,18:19],cut_mean=10,small_val=1) # 
sapply( xtb_dp_diff<- list( names(dp_cmp)[dp_cmp> 1], names(dp_cmp)[dp_cmp< -1] ),length) # Down and Up
xtb_dp_glo<-xtb_dp_diff
length(xtb_dp_glo[[1]]<-unique(hs_orth2[ hs_orth2[,3] %in% xtb_dp_glo[[1]] & hs_orth2[,4] %in% xts_s16_glo, 3]))
length(xtb_dp_glo[[2]]<-unique(hs_orth2[ hs_orth2[,3] %in% xtb_dp_glo[[2]] & hs_orth2[,4] %in% xts_s16_glo, 3]))


####################
## Systemically changing genes in zebrafish (by bulk and scRNA-seq)
####################
# Data: White et al. 2017 and Wagner et al. 2018
# expression matrix in White et al. 2017
raw<- read.table('data/lin28_data/white_zebrafish_elife-30860-supp3-v1.tsv',sep='\t',as.is=T,header=T,comment.char='',quote='') # It can be download from https://elifesciences.org/articles/30860
st<- sapply( colnames(raw)[-(1:8)], function(x){ substring(x,1,nchar(x)-2) })
dr_all_exp<-t(apply( raw[,-(1:8)], 1, function(x){ tapply(x, st, mean) }))[,unique(st)]
rownames(dr_all_exp)<- raw[,1]

# systemically expressed genes in scRNA-seq
#!# Because the size is too large, raw data of Wagner et al. 2018 is not included in this repertoire. It can be downloaded from GSE112294
raw<- mapply( function(x,y){
  file_name<-paste(rp, 'GSM',y,'_', x,'hpf_nm.csv',sep='')
  mx<- read.csv(file_name)
  mx2<- mx[,-1]
  rownames(mx2)<-mx[,1]
  return(mx2)
},x=c('04','06','08','10','14','18','24'), y=3067189:3067195 )
tp48_mn<-mapply( function(x,y){ 
  res<-t(apply(x, 1, function(z){ tapply(z, y, mean) }))
  return(res)
}, x=raw[c(1,3)], y=clu_info[c(1,3)] )
drs_s8_glo<- rownames(tp48_mn[[2]])[apply( tp48_mn[[2]]>=0.5, 1, sum)>= ceiling(ncol(tp48_mn[[2]])*.5)]
s14_in_tp<- names(table(clu_info[[5]]))[ table(clu_info[[5]])>40 ]
s14_mn<-sapply( s14_in_tp , function(x){
  s2<- raw14[[1]][, clu_info[[5]]== x ]
  apply(s2,1,mean) } )
drs_s14_glo<- rownames(s14_mn)[apply( s14_mn>=0.5, 1, sum)>= ceiling(ncol(s14_mn)*.5)]

# Stage 1 vs. stage 2: DEGs between 6~8 hpf and 3~4.3 hpf, and systemically expressed in scRNA-seq (8 hpf)
up_cmp<- cal_wilp2( mx1=dr_all_exp[,4:5], mx2=dr_all_exp[,7:8],cut_mean=10,small_val=1) # 
sapply( drb_up_diff<- list( names(up_cmp)[up_cmp> 1], names(up_cmp)[up_cmp< -1] ),length) # Down and Up
drb_up_glo<-drb_up_diff
length(drb_up_glo[[1]]<-unique(hs_orth2[ hs_orth2[,5] %in% drb_up_glo[[1]] & hs_orth2[,6] %in% drs_s8_glo, 5])) # 128
length(drb_up_glo[[2]]<-unique(hs_orth2[ hs_orth2[,5] %in% drb_up_glo[[2]] & hs_orth2[,6] %in% drs_s8_glo, 5])) # 385

# Stage 3 vs. stage 4: DEGs between 3~4 d and 10.3~16 hpf, and systemically expressed in scRNA-seq ( hpf)
dp_cmp<- cal_wilp2( mx1=dr_all_exp[,9:10], mx2=dr_all_exp[,16:17],cut_mean=10,small_val=1) # 
sapply( drb_dp_diff<- list( names(dp_cmp)[dp_cmp> 1], names(dp_cmp)[dp_cmp< -1] ),length) # Down and Up
drb_dp_glo<-drb_dp_diff
length(drb_dp_glo[[1]]<-unique(hs_orth2[ hs_orth2[,5] %in% drb_dp_glo[[1]] & hs_orth2[,6] %in% drs_s14_glo, 5]))
length(drb_dp_glo[[2]]<-unique(hs_orth2[ hs_orth2[,5] %in% drb_dp_glo[[2]] & hs_orth2[,6] %in% drs_s14_glo, 5]))


####################
## Systemically changing genes in 4 species
####################
# unify namespaces
hs_dp_glo<-hs_glo
dr_dp_glo<- drb_dp_glo
dr_up_glo<- drb_up_glo
xt_up_glo<- xtb_up_glo
xt_dp_glo<- xtb_dp_glo
# change to human genes by orthlogs
to_hs<- function( gl, ind=0, name='',is_out=F, orth=hs_orth2 ){
  file1<- paste('result/gene/',name,'_DN.txt',sep='')
  file2<- paste('result/gene/',name,'_UP.txt',sep='')
  if( ind==0 ){ 
    g1<- sort(ge2an[gl[[1]]])
    g2<- sort(ge2an[gl[[2]]])
  }else{  
    # take the shortest human symbol for multiple orthologs
    g1<- unique(unlist(sapply( gl[[1]], function(x){ res<-unique(setdiff(orth[ orth[,ind]==x, 2],'')); if(length(res)<=1) return(res) else return( res[order(nchar(res))][1] )  })))
    g2<- unique(unlist(sapply( gl[[2]], function(x){ res<-unique(setdiff(orth[ orth[,ind]==x, 2],'')); if(length(res)<=1) return(res) else return( res[order(nchar(res))][1] )  })))
  }
  if(is_out){
    write( sort(g1), file=file1)
    write( sort(g2), file=file2)
  }
  return( list(g1,g2) )
}
allg_hs<-mapply( to_hs, gl=list(hs_dp_glo,mm_dp_glo,mm_up_glo,xt_dp_glo,xt_up_glo,dr_dp_glo,dr_up_glo), ind=c(0,7,7,3,3,5,5), name=c('hs_dp','mm_dp','mm_up','xt_dp','xt_up','dr_dp','dr_up'), MoreArgs=list(is_out=F, orth=hs_orth), SIMPLIFY=F)
names(allg_hs)<-c('hs_dp','mm_dp','mm_up','xt_dp','xt_up','dr_dp','dr_up')
