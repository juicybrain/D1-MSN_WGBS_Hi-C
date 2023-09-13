### This file lists the key parameters used in the WGBS alignment and methylation call workflow
### trimming
trimgalore --paired  --clip_r1 10 --three_prime_clip_r1 15  $sample\_R1.fq.gz --clip_r2 15 --three_prime_clip_r2 3 $sample\_R2.fq.gz

### alignment to lambda DNA and conversion rate calculation
bismark   --bowtie2 --gzip --parallel 8  --bam  /home/yli/ref/genome/lambdaDNA -1 $wd/$sample\_R1_val_1.fq.gz -2 $wd/$sample\_R2_val_2.fq.gz -o ./${sample}_bismark_bowtie_lambda
# deduplication
deduplicate_bismark --bam ${sample}*.bam
#filter non-conversion reads
samtools view ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | awk 'BEGIN{OFS="\t"}{split($14,a,":")}{if(a[3] ~ /[a-z]/); else if(gsub(/[XH]/, "", a[3])>2)print $0}' | cut -f 1 | sort | uniq -c |awk '$1==2' |cut -f 2  > ${spl}\.list
samtools view -h ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | grep -vf  ${spl}\.list | samtools view -bS -o ${spl}\_filter.bam

# call methylation
bismark_methylation_extractor -p --bedGraph --gzip --CX --ignore_3prime 3 --ignore_3prime_r2 3 --cytosine_report --genome_folder ~/ref/genome/lambdaDNA/  ${spl}\_filter.bam
      zcat *.CX_report.txt.gz | grep $'\t'CG[ATCG]$ > CG.lambda.cov
      zcat *.CX_report.txt.gz | grep $'\t'CC[ATCG]$ > CC.lambda.cov
      zcat *.CX_report.txt.gz | grep $'\t'CT[ATCG]$ > CT.lambda.cov
      zcat *.CX_report.txt.gz | grep $'\t'CA[ATCG]$ > CA.lambda.cov
      
      cat CA.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CA", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
      cat CC.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CC", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
      cat CT.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CT", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
      cat CG.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CG", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt

### alignment to mm10 reference genome
bismark --bowtie2 --gzip -X 1000  --parallel 4 --bam  /home/yli/ref/genome/male.mm10  -1  $sample\_merge.R1.fq.gz  -2 $spl\_merge.R2.fq.gz  -o ${spl}\_bismark_bowtie_mm10
#deduplicaiton and filtering non-conversion reads
deduplicate_bismark --bam ${spl}*pe.bam
samtools view ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | awk 'BEGIN{OFS="\t"}{split($14,a,":")}{if(a[3] ~ /[a-z]/); else if(gsub(/[XH]/, "", a[3])>2)print $0}' | cut -f 1 | sort | uniq -c |awk '$1==2' |cut -f 2  > ${spl}\.list
samtools view -h ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | grep -vf  ${spl}\.list | samtools view -bS -o ${spl}\_filter.bam

#call methylation
bismark_methylation_extractor -p --bedGraph --gzip --parallel 8  --CX  --cytosine_report --genome_folder /home/yli/ref/genome/male.mm10  ${spl}\_filter.bam
