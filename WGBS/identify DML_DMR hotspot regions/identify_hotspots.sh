###  Identify broad regions enriched with DMLs
#
#  Author  Yuxiang Li
#

### Count DML numbers in each 10k-window bin of all genomic regions covered by WGBS (we excluded regions with less than 20% of average CG coverage) 
for spl in `ls dml*bed`; do   bedtools map -a <(sort -k1,1 -k2,2n  /home/yli/ref/bed/windows/D1_WGBS_covered_region.10kW5kS.auto.bed) -b <(sort -k1,1 -k2,2n ${spl}) -c 4 -o sum | awk -v OFS="\t" '{if($4>0)print $0; else print $1,$2,$3,0}' > ${spl%".bed"*}\.dml.counts; done
wait

### Calculate the p-value of each bin 
for spl in `ls *counts`; do  Rscript poisson.r ${spl}; done

wait


### use comb-p to combine regional p-values
func1(){
    cat ${1} |awk -v FS=" " -v OFS="\t" '{print $1,$2,$3,$6}' > ${1}\.comb-p
            wait
    python ~/git/combined-pvalues/cpv/peaks.py --seed ${2} --dist 200000  ${1}\.comb-p > ${1}\.res
}

for spl in `ls *txt`; do cat ${spl} | awk   '$4==0' | head -n 1 | awk -v name=${spl} '{printf name"\t"$6"\n"}' >> seed.txt; done

cat seed.txt | while read line; do func1 $line; done

for spl in `ls *txt`; do cat ${spl}|awk -v FS=" " -v OFS="\t" '{print $1,$2,$3,$6}' > ${spl}\.comb-p;  python ~/git/combined-pvalues/cpv/peaks.py --seed 0.01 --dist 100000  ${spl}\.comb-p > ${spl}\.res; done

wait
for spl in `ls *res`; do cat $spl | awk '$4<1e-10' > ${spl}\.1e-10.bed; done
