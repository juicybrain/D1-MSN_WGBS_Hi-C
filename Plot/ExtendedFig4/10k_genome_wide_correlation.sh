#
#   Author Yuxiang Li
#
#### Genome-wide CG and non-CG methyaltion level correlation

    multiBigwigSummary bins --numberOfProcessors 10 --bwfiles <bwfiles>  --chromosomesToSkip chrX chrY    --numberOfProcessors 16 -out CG.nonCG.10kbin.genomewide.npz --outRawCounts MF.CG_nonCG.10kbin.tab --labels <labels>
    plotCorrelation -in MF.CG_nonCG.1k_10kbin.genomewide.npz  --corMethod pearson --skipZeros --plotTitle "Pearson Correlation" --whatToPlot heatmap   -o MF.CG_nonCG.10kbin.npz.pearson.g.svg --outFileCorMatrix MF.CG_nonCG.10kbin.npz.tab  --plotHeight 10 --plotWidth 10


