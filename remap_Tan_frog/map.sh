for sample in $(cat sam_name.txt);do
    tophat2 --library-type fr-firststrand -p 8 -G Xenopus_tropicalis.Xenopus_tropicalis_v9.1.100.gtf -o tophat_output/${sample} xt91/xt91 ../raw/${sample}_1.fastq ../raw/${sample}_2.fastq
    echo finish_${sample}
done
