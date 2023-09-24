### Processing and alignment of Hi-C sequencing reads
#
#  Author Yuxiang Li
#

# Trim adaptor
trim_galore --paired   ${spl}\_R1.fastq.gz ${spl}\_R2.fastq.gz

#  Align reads with Hi-C Pro with the following parameters in config-hicpro.txt
#  N_CPU =24
#  MIN_MAPQ = 10 
#  BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
#  BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder
#  REFERENCE_GENOME = male.mm10/female.mm10 (remove chrY for female samples)
#  LIGATION_SITE = AAGCTAGCTT
#  MIN_FRAG_SIZE = 100
#  MAX_FRAG_SIZE = 100000
#  MIN_INSERT_SIZE = 100
#  MAX_INSERT_SIZE = 1000
#  GET_ALL_INTERACTION_CLASSES = 1
#  RM_SINGLETON = 1
#  RM_DUP = 0
HiC-Pro -i /raw_data -o ./hicpro_out -c config-hicpro.txt 

# Create ".hic" files for JuicerTools HiC-Pro valid pairs 
hicpro2juicebox.sh -i ${spl} -g male.mm10.natural.genomeSize -j juicer_tools_1.22.01.jar
# identify TAD at 25k resolution for merged four groups Hi-C data: male Con, male KO, female Con, female KO
java -jar  juicer_tools_1.22.01.jar arrowhead -r 25000 -k KR --ignore-sparsity ${spl}\.hic ${spl}\.25k.KR

# Identify loop at 10000, 25000, 50000 resolution
java -jar juicer_tools_1.22.01.jar hiccups  -k KR -r 5000,10000,25000 -f 0.1 --ignore-sparsity ${spl}  ${spl}\_10.25.50k --cpu
java -jar juicer_tools_1.22.01.jar hiccups  -k KR -r 50000 -f 0.1 --ignore-sparsity ${spl}  ${spl}\_10.25.50k --cpu

# create Homer tag depository
  cut -f 1-7 validPairs > validPairs.homer
  makeTagDirectory ${name} -format HiCsummary ${spl}
# calculate the loop score for each Hi-C replicate
findTADsAndLoops.pl score  -loop ../D1_loop.bedpe  -o male.loop   -d <Hi-C replicates directory> -cpu 4 -res 25000  -window 50000 -normTotal 1e8 -raw

# convert ".hic" file to ".cool" file using Hicexplorer
hicConvertFormat   --matrices ${spl}  --outFileName ${spl}\.mcool --inputFormat hic --outputFormat cool

# get matrix of certain resolution
hicConvertFormat --matrices ${spl}\.mcool::/resolutions/25000 --outFileName ${spl}\_250k_${i}\.cool --inputFormat cool --outputFormat cool 

# mask bad regions (downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz and http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/genomicSuperDups.txt.gz )
hicAdjustMatrix -m ${spl} --action mask --regions   male.mm10.badrigion.GT25k.bed  -o mask/${spl}\.mask.cool

# normalize the reads of Hi-C replicates
current_time=$(date +%Y.%m.%d-%H.%M.%S)
hicNormalize -m  <all replicates>  --normalize smallest  -o <normalized files>  2>${current_time}\.normalTosmallest.log  

# Compare Tet1 KO to Tet1 Con Hi-C matrices
hicCompareMatrices -m ${spl}\_250k.cool.mask.nrl.cool F_N_250k.cool.mask.nrl.cool --outFileName ${spl}\_nrl_PosMinusNeg.log2ratio.cool --operation log2ratio

# KR correction of Hi-C matrices
hicCorrectMatrix correct --matrix  ${spl}\_25k.mask.cool   --correctionMethod  KR -o ${spl}\_25k.mask.KR.cool

# Hi-C A/B compartment analysis PC1 analysis with hicPCA, K27ac ChIPseq data as marker for active chromatin

hicPCA --matrix ../${spl}\.2Rep.allValidPairs.hic.mcool::/resolutions/1000000    --outputFileName  ${spl}\.distalN.1M.pca.txt  --whichEigenvectors 1  --format bedgraph --method dist_norm  --extraTrack  GSE197668_MSD1K27A.bw  --histonMarkType active
hicPCA --matrix ../${spl}\.2Rep.allValidPairs.hic.mcool::/resolutions/1000000    --outputFileName  ${spl}\.distalN.1M.pca.txt  --whichEigenvectors 1  --format bedgraph --method dist_norm  --extraTrack  GSE197668_FSD1K27A.bw  --histonMarkType active




