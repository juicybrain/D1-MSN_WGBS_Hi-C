for spl in  M  F
do
    computeMatrix reference-point -S ${spl}\N.CG.cov.bw  -R ${spl}\_D1_Sal_fosb_summits.valid.bed.top500 ${spl}\_D1_Sal_fosb_summits.valid.bed.bottom500  ${spl}\_D1_Coc_fosb_summits.valid.bed.top500 ${spl}\_D1_Coc_fosb_summits.valid.bed.bottom500   --binSize 50 -b 2000 -a 2000 -o ${spl}.ref.rp.gz
# computeMatrix reference-point -S *.bw -R /home/yli/ref/bed/Striatum.mm10.enhancer/Striatum.mm10.s.bed   --binSize 10 -b 5000 -a 5000 -o WGBS_EMseq_Str_EnHr.rp.gz
# computeMatrix scale-regions -S D1[FM]_delta.50bin.bw -R chr6_DHR.bed --binSize 1000 -b 10000000 -a 10000000 -o chr6_DHR.delta.sr.gz

#  plotProfile -m WGBS_EMseq_Str_EnHr.rp.gz --perGroup  --colors red red blue blue  --plotHeight 21  --plotWidth 40 -out WGBS_EMseq_EnHr.svg
#     plotProfile -m ${spl}.rp.gz --perGroup --colors orange orange purple purple  --startLabel DMR --endLabel DMR --plotHeight 21  --plotWidth 60 -out ${spl}\.scale.FosB.in.DMR.svg
    plotProfile -m ${spl}.ref.rp.gz   --startLabel DMR --endLabel DMR --plotHeight 21  --plotWidth 60 -out ${spl}\.ref.FosB_Peak.meCG.svg

     plotHeatmap -m  ${spl}.ref.rp.gz   -out ${spl}.ref.heatmap.Fosb.Peak.meCG.svg --colorMap hot  --zMin 0 --zMax 1

#  plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar'  --kmeans 4


done
