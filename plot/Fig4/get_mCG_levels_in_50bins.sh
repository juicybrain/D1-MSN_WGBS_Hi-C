
for index in 1 2 3 4 5 6 7 8 9 10

do

            mkdir -p promoter
            mkdir -p gene_body
            mkdir -p downstream
            for quantile in     chr8.DHZ chr14.DHZ
            do
   		bedtools map -a <(sort -k1,1 -k2,2n promoter/${quantile}\_promoter_${index}.bed) -b D1_all_sample.DSS -c 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52  -o sum|awk '{for (i=1; i<7;i++)printf("%s\t", $i)}{for (i = 7; i < NF; i += 2) {if ($(i+1) > 0) { result = $i / $(i+1);printf( "%s\t", result);} else { printf("NA\t", $i, $(i+1));} } printf("\n"); }'  >  promoter/${quantile}\_promoter_${index}.bedgraph &

   		bedtools map -a <(sort -k1,1 -k2,2n downstream/${quantile}\_downstream_${index}.bed) -b D1_all_sample.DSS -c 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52  -o sum | awk '{for (i=1; i<7;i++)printf("%s\t", $i)}{for (i = 7; i < NF; i += 2) {if ($(i+1) > 0) { result = $i / $(i+1);printf( "%s\t", result);} else { printf("NA\t", $i, $(i+1));} } printf("\n"); }'  >  downstream/${quantile}\_downstream_${index}.bedgraph &

             done
             wait
done


for index in {1..50}
do
            mkdir -p gene_body
            for quantile in chr8.cDMR chr14.cDMR
            do

           bedtools map -a <(sort -k1,1 -k2,2n gene_body/${quantile}\_gene_body_${index}.bed ) -b D1_all_sample.DSS -c 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52  -o sum|awk '{for (i=1; i<7;i++)printf("%s\t", $i)}{for (i = 7; i < NF; i += 2) {if ($(i+1) > 0) { result = $i / $(i+1);printf( "%s\t", result);} else { printf("NA\t", $i, $(i+1));} } printf("\n"); }'  >  gene_body/${quantile}\_gene_body_${index}.bedgraph &
            done
            wait
done
