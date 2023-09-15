#
#  Author Yuxiang Li
#
### prepare the methylation level in 51 scales from 1kb to 20Mb

Draw_windows(){
  j=1
  for i in  1       2       3       4       5       6       7       8       10      12      15      18      22      26      31      38      46      55      66      79      95      114     137     164     197     237     284     341     410 492      590     708     850     1020    1224    1469    1763    2116    2539    3047    3657    4388    5266    6319    7583    9100    10920   13104   15725   18870   20000
    do
    bedtools makewindows -g  /home/yli/ref/genome/male.mm10/male.mm10.genomeSize -w `expr 1000 \* ${i}` >male.mm10.${j}\.bed &
    j=$((j+1))

done
}

    mkdir counts
    for spl in `ls male.mm10*.bed`; do bedtools map -a <(sort -k1,1 -k2,2n ${spl}) -b ../D1_M_F.merge_MF_NPreplicates.DSS  -c 4,5,6,7,8,9,10,11 -o sum > counts/${spl}\.counts.bed&  done
    mkdir -p pixel
    cd counts
    for spl in `ls *.bed.counts.bed`; do bedtools intersect -wa -wb -a ../male.mm10.1.bed -b ${spl} | awk -v OFS="\t" '{if($7*$9*$11*$13>0)print $1,$2,$3,$10/$9-$8/$7,$14/$13-$12/$11;else print $1,$2,$3,"NA","NA"}' | sort -k1,1 -k2,2n|uniq  > ../pixel/${spl}\.meCG.bed& done
