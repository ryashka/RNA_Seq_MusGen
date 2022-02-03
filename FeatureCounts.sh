featureCounts -M -T 64 -p -a \
~/Documents/Programming/reference_files/mm10_star_reference_files/gencode.vM25.primary_assembly.annotation.gtf  \
-t exon -g gene_id -o ~/Documents/Programming/Keith_seq/processed_data/output_files/featureCounts_files/combinedFeatureCounts.txt \
~/Documents/Programming/Keith_seq/processed_data/output_files/bam_files/*.bam
