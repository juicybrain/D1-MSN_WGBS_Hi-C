#
#  Author Yuxiang Li
#
### plot the meCG level in DHZs

M_P="D1M_WGBS_7_CG.DSS.gz.bw        D1M_WGBS_8_CG.DSS.gz.bw        D1M_WGBS_9_CG.DSS.gz.bw        D1M_WGBS_10_CG.DSS.gz.bw    D1M_WGBS_12_CG.DSS.gz.bw     "
M_N="D1M_WGBS_1_CG.DSS.gz.bw        D1M_WGBS_2_CG.DSS.gz.bw        D1M_WGBS_3_CG.DSS.gz.bw        D1M_WGBS_4_CG.DSS.gz.bw     D1M_WGBS_5_CG.DSS.gz.bw        D1M_WGBS_11_CG.DSS.gz.bw "
F_N="D1F_WGBS_1_CG.DSS.gz.bw        D1F_WGBS_2_CG.DSS.gz.bw        D1F_WGBS_4_CG.DSS.gz.bw        D1F_WGBS_8_CG.DSS.gz.bw     D1F_WGBS_12_CG.DSS.gz.bw       "
F_P="D1F_WGBS_3_CG.DSS.gz.bw        D1F_WGBS_5_CG.DSS.gz.bw        D1F_WGBS_6_CG.DSS.gz.bw        D1F_WGBS_7_CG.DSS.gz.bw     D1F_WGBS_10_CG.DSS.gz.bw       D1F_WGBS_11_CG.DSS.gz.bw   "

size=6
region=100000
bin=10000

computeMatrix scale-regions -S  ${M_P}  ${M_N}     -R D1_M_WTvsCre.${diff}.cDMR.bed     --regionBodyLength ${region}  --binSize ${bin} -b ${region} -a ${region} -o M_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.gz
computeMatrix scale-regions -S  ${F_P}  ${F_N}     -R D1_F_WTvsCre.${diff}.cDMR.bed     --regionBodyLength ${region}  --binSize ${bin} -b ${region} -a ${region} -o F_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.gz
computeMatrix scale-regions -S   ${M_N}  ${F_N}    -R D1_WT_Male_vs_Female.${diff}.cDMR.bed     --regionBodyLength ${region}  --binSize ${bin} -b ${region} -a ${region} -o WT.${CX}\_${diff}\.${bin}\.${region}\.gz
computeMatrix scale-regions -S  ${M_P} ${F_P}      -R D1_KO_Male_vs_Female.${diff}.cDMR.bed     --regionBodyLength ${region}  --binSize ${bin} -b ${region} -a ${region} -o KO.${CX}\_${diff}\.${bin}\.${region}\.gz

wait
for diff in hyper hypo
do
    for CX in CG
    do
    plotProfile   --perGroup -m M_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.gz  --colors  red red red red red blue blue blue blue blue blue  --samplesLabel M_KO_Rep1 M_KO_Rep2 M_KO_Rep3 M_KO_Rep4 M_KO_Rep5  M_WT_Rep1 M_WT_Rep2 M_WT_Rep3 M_WT_Rep4 M_WT_Rep5  M_WT_Rep6    --plotHeight ${size}  --plotWidth ${size} -out M_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.${size}.svg --startLabel ""  --endLabel ""
    
    plotProfile   --perGroup -m F_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.gz  --colors red red red red red red blue blue blue blue blue --samplesLabel  F_KO_Rep1 F_KO_Rep2 F_KO_Rep3 F_KO_Rep4 F_KO_Rep5 F_KO_Rep6 F_WT_Rep1 F_WT_Rep2 F_WT_Rep3 F_WT_Rep4 F_WT_Rep5    --plotHeight ${size} --plotWidth ${size} -out F_WTvsCre.${CX}\_${diff}\.${bin}\.${region}\.${size}.svg   --startLabel ""  --endLabel ""

        plotProfile   --perGroup -m WT.${CX}\_${diff}\.${bin}\.${region}\.gz  --colors  red red red red red red blue blue blue blue blue  --samplesLabel   M_WT_Rep1 M_WT_Rep2 M_WT_Rep3 M_WT_Rep4 M_WT_Rep5  M_WT_Rep6  F_WT_Rep1 F_WT_Rep2 F_WT_Rep3 F_WT_Rep4 F_WT_Rep5     --plotHeight ${size} --plotWidth ${size} -out WT.${CX}\_${diff}\.${bin}\.${region}\.${size}.svg  --startLabel ""  --endLabel ""

            plotProfile   --perGroup -m KO.${CX}\_${diff}\.${bin}\.${region}\.gz  --colors   red red red red red blue blue blue blue blue blue  --samplesLabel  M_KO_Rep1 M_KO_Rep2 M_KO_Rep3 M_KO_Rep4 M_KO_Rep5 F_KO_Rep1 F_KO_Rep2 F_KO_Rep3 F_KO_Rep4 F_KO_Rep5 F_KO_Rep6    --plotHeight ${size} --plotWidth ${size} -out KO.${CX}\_${diff}\.${bin}\.${region}\.${size}.svg  --startLabel ""  --endLabel ""
            
    done
done

### Count the intersection between non-CG DML and CG DHZs

  bedtools intersect -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4>0' ) -b M.hyper.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4<0' ) -b M.hyper.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4>0' ) -b M.hypo.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4<0' ) -b M.hypo.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4>0' ) -a M.hyper.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4<0' ) -a M.hyper.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4>0' ) -a M.hypo.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph | awk '$4<0' ) -a M.hypo.cDMR.bed| sort -k1,1 -k2,2n | uniq |wc -l  
  bedtools intersect -v -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph |awk '$4<0'  ) -b M.hyper.cDMR.bed M.hypo.cDMR.bed | sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -v -wa -a <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph |awk '$4>0'  ) -b M.hyper.cDMR.bed M.hypo.cDMR.bed | sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -v -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph   ) -a M.hyper.cDMR.bed  | sort -k1,1 -k2,2n | uniq |wc -l
  bedtools intersect -v -wa -b <(cat dml_D1_Smth_M_CtlvsKo_5p1e-4.nonCG.bedgraph   ) -a M.hypo.cDMR.bed | sort -k1,1 -k2,2n | uniq |wc -l






















