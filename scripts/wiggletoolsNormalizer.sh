#!/usr/bin/sh

# 0. Given files:
neg.bam
pos.bam
map.wig
neg.trimFill.bw
pos.trimFill.bw

# 1. Calculate fragment length.
wiggletools_max_shift.py

# 2. Shift read starts in 3' direction by Li/2 (or shift neg read starts by Li)
wiggletools write neg_shifted.wig shiftPos Li neg.bam

# 3. Compute read-start coverage for pos.bam, already have it for neg_shifted.wig
wiggletools write pos.wig pos.bam

# 4. Compute smooth weighted sum of read counts, using a window based on max. estimated frag. len. for that dataset.
#    Equivalent to approximate number of fragments overlapping each position.
wiggletools write pos_smooth.wig smooth 300 pos.wig
wiggletools write neg_smooth_shifted.wig smooth 300 neg_shifted.wig

# 5. Add together fragment counts from both strands.
wiggletools write sum.wig sum pos_smooth.wig neg_smooth_shifted.wig

# 6. Compute cumulative mappability.
wiggletools write map_smooth.wig smooth 300 map.wig

# 7. Compute expected fragment counts assuming uniform distribution over all uniquely mappable locations.
#    samtools view -c -F 4 bam for total_mapped_reads_in_bam???
#    wiggletools AUC map.wig for total_mappable_genome_size
wiggletools write expected.wig scale total_mapped_reads_in_bam/total_mappable_genome_size map_smooth.wig

# 8. Normalized signal.
wiggletools write norm.wig ratio strict sum.wig expected.wig

# 9. Trim regions less than 0.25 of max expected.
wiggletools write filtered_norm.wig trim gt maxI expected.wig scale 4 expected.wig norm.wig

# Can be summed up as:

wiggletools write sum.wig sum smooth 300 pos.bam smooth 300 shiftPos Li neg.bam
wiggletools write expected.wig scale total_mapped_reads_in_bam/total_mappable_genome_size smooth 300 map.wig
wiggletools write_bg norm.bg trim gt maxI expected.wig scale 4 expected.wig ratio sum.wig expected.wig

# Trim the output.
awk 'FNR==NR{a[$1]=$2;next}{if ($2 <= a[$1]) print($3 <= a[$1])?$0:($1"\t"$2"\t"a[$1]"\t"$4)}' ../../examples/yeast-genome.fa.fai norm.bg > norm_trimmed.bg
