processData () {
local_sn=$1
directory_name=`sed "${local_sn}q;d" /Users/smolkinr/Documents/Programming/Keith_seq/input_files/directories.txt`

gatk FastqToSam \
-F1 input_data/raw_data/${directory_name}/*_1.fq.gz \
-F2 input_data/raw_data/${directory_name}/*_2.fq.gz \
-O processed_data/sam_files/"${directory_name}_fastqtosam.bam" \
-SM $sn

gatk --java-options "-Xmx300G" MarkIlluminaAdapters \
-I processed_data/sam_files/"${directory_name}_fastqtosam.bam" \
-O processed_data/sam_files/"${directory_name}_markilluminaadapters.bam" \
-M processed_data/sam_files/"${directory_name}_markilluminaadapters_metrics.txt"

gatk --java-options "-Xmx300G" ClipReads \
-I processed_data/sam_files/"${directory_name}_markilluminaadapters.bam" \
-O processed_data/sam_files/"${directory_name}_clipped.bam" \
-QT 10

gatk --java-options "-Xmx300G" SamToFastq \
-I processed_data/sam_files/"${directory_name}_clipped.bam" \
--FASTQ processed_data/fastq_files/"${directory_name}_samtofastq_R1.fq" \
--SECOND_END_FASTQ processed_data/fastq_files/"${directory_name}_samtofastq_R2.fq" \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE false --NON_PF true


STAR \
--genomeDir ~/Documents/Programming/reference_files/mm10_star_reference_files/star_genome \
--readFilesIn processed_data/fastq_files/"${directory_name}_samtofastq_R1.fq" processed_data/fastq_files/"${directory_name}_samtofastq_R2.fq" \
--outFileNamePrefix processed_data/sam_files/"${directory_name}_" \
--runThreadN 1000

gatk --java-options "-Xmx300G" MergeBamAlignment \
--ALIGNED_BAM processed_data/sam_files/"${directory_name}_Aligned.out.sam" \
--UNMAPPED_BAM processed_data/sam_files/"${directory_name}_fastqtosam.bam" \
--OUTPUT processed_data/sam_files/"${directory_name}_not_piped.bam" \
-R ~/Documents/Programming/reference_files/mm10_star_reference_files/GRCm38.primary_assembly.genome.fa --CREATE_INDEX true --ADD_MATE_CIGAR true \
--CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS

rm processed_data/sam_files/"${directory_name}_Aligned.out.sam"
rm processed_data/sam_files/"${directory_name}_fastqtosam.bam"

gatk --java-options "-Xmx300G" MarkDuplicatesSpark \
-I processed_data/sam_files/"${directory_name}_not_piped.bam" \
-O processed_data/output_files/"${directory_name}_dms_with_clipping.bam" \
--remove-sequencing-duplicates

rm processed_data/sam_files/*
rm processed_data/fastq_files/*

}
for sn in {1..24}
do
    processData $sn
done
