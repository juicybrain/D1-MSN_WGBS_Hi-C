#!/bin/bash


# binomial test
btest(){
/usr/bin/Rscript binomial_test_fast.r   ${1}\_CA.gz ${2}
}

# binomial test for simulated data
btest_simu(){
/usr/bin/Rscript binomial_test_fast.r   ${1}\_CA.gz.total.simu.cov  ${2}
}


# simulation data with same coverage as the real WGBS data
simu(){
/usr/bin/Rscript binomial_simulation.r  ${1}\_CA.gz ${2}
}

# Count positive mC with different p.values
count_p(){
cat ${1}\_CA.gzbinomial_test.cov |sed '1d'| cut -f 8|sort -g | uniq -c | awk -v FS=" " -v OFS="\t" '{print $1,$2}' > ${1}\.pCount.txt
cat ${1}\_CA.gz.total.simu.covbinomial_test.cov |sed '1d' | cut -f 8|sort -g | uniq -c | awk -v FS=" " -v OFS="\t" '{print $1,$2}' > ${1}\.simu.pCount.txt
}


# calculate fdr at different p.value threshold
cal_fdr(){
/usr/bin/Rscript fdr_calculation.r  ${1}\.pCount.txt  ${1}\.simu.pCount.txt
}

# simulation with binomial model
cat binomial_test.txt| sed '1d' | cut -f 1,3  | while read line
       do
       simu ${line}
       done

wait

# binomial_test for real data and simulated data
cat binomial_test.txt| sed '1d' | cut -f 1,3  | while read line
    do
    btest ${line}
    btest_simu ${line}
    done

wait

# calculate fdr and select p.value threshold which makes fdr<0.01
cat binomial_test.txt| sed '1d' | cut -f 1,3  | while read line
do
    count_p ${line}
    wait
    cal_fdr ${line}
done

# read file with selected p.value as column #1 and file name as column #2
check_fdr(){
zcat ${2}\_CA.gzbinomial_test.cov.gz | sed '1d' | awk -v a=${1} -v OFS="\t" '{if($8<a) print $0; else print $1,$2,$3,0,$5+$4,$6,$7,$8}' | gzip > ${2}\.fdr_fixed.CA.gz
}

# get mCA sites
cat fdr.res.txt | while read line
do
check_fdr ${line}
done

