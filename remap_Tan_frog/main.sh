# Main script: map frog data (Tan et al. 2013) to frog genome in Ensembl100
# download Xenopus_tropicalis.Xenopus_tropicalis_v9.1.100.gtf and Xenopus_tropicalis.Xenopus_tropicalis_v9.1.dna.toplevel.fa from ensembl FTP. Modify date: 03/06/2020
# see download.sh

# build genome index
Bowtie2-build ../Xenopus_tropicalis.Xenopus_tropicalis_v9.1.dna.toplevel.fa xt91 > build.out 2>&1 &

# map samples in 2 batches
nohup ./map.sh > map.out 2>&1 & 
nohup ./map2.sh > map2.out 2>&1 & 

# read count by htseq
mkdir htseq_output
nohup ./count.sh > count.out 2>&1 & 
