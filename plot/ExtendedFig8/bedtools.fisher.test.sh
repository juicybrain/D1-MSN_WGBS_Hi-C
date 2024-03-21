#
# bedtool fisher test calculate the significance of overlop between D1-MSN CG DMRs and D1-MSN Hi-C loops, genome size adjusted for excluding uncovered regions (bad regions) by WGBS
#
bedtools fisher -a DMR.bed   -b  loop_anchor.bed -g WGBS.covered.genome.size