for sample in $(cat all_sam_name.txt);do
  python2.7 -m HTSeq.scripts.count -q -f bam -r pos -s reverse tophat_output/${sample}/accepted_hits.bam Xenopus_tropicalis.Xenopus_tropicalis_v9.1.100.gtf > htseq_output/${sample}.txt
  echo ${sample}
done
