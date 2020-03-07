// GATK4 Variant Calling Pre Processing Pipeline
 
params.reads = "/data/input_fastq/HT27VBGXC_n0{1,2}_{5}*.fastq*"
params.ref = "/data/Protist/Plasmodium_vivax/PlasmoDB/Sal1/PlasmoDB-46_PvivaxSal1_Genome.fasta"
params.outdir = "/data/gatk4/out/"
params.tmpdir ="/tmp"

println "reads: $params.reads"
println "ref: $params.ref"
println "outdir: $params.outdir"

ref = file(params.ref)

// Prepare the fastq read pairs for input.
// Use the size parameter to not auto-group, and instead
// use the mapping through getBaseName() and subtract
// two regexs to get the ID.
// This is custom for NYU CGSB sequence data file naming format
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/_n0[12]/ - ~/.fastq/ }
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .set { read_pairs_ch }

process align {
    publishDir "${params.outdir}/aligned_reads", mode:'copy'

    input:
    set pair_id, file(reads) from read_pairs_ch

    output:
    set val(pair_id), file("${pair_id}_aligned_reads.sam") into aligned_reads_ch

    script:
    readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${pair_id}"
    """
    bwa mem -K 100000000 -v 3 -t ${task.cpus} -Y \
	-R \"${readGroup}\" ${ref} ${reads[0]} ${reads[1]} > ${pair_id}_aligned_reads.sam
    """
}

process markDuplicatesSpark {
    publishDir "${params.outdir}/dedup_sorted", mode:'copy'

    input:
    set val(pair_id), file(aligned_reads) from aligned_reads_ch

    output:
    set val(pair_id), \
	file("${pair_id}_sorted_dedup.bam") \
	into sorted_dedup_ch, sorted_dedup_ch_for_metrics
    set val(pair_id), \
	file ("${pair_id}_dedup_metrics.txt") \
	into dedup_metrics_ch

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.tmpdir}" \
	 MarkDuplicatesSpark \
	-I ${aligned_reads} \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam \
    """
}

process getMetrics {
    publishDir "${params.outdir}/metrics", mode:'copy'

    input:
    set val(pair_id), \
	file(sorted_dedup_reads) \
	from sorted_dedup_ch_for_metrics

    output:
    set val(pair_id),
            file("${pair_id}_alignment_metrics.txt"), \
            file("${pair_id}_insert_metrics.txt"), \
            file("${pair_id}_insert_size_histogram.pdf"), \
            file("${pair_id}_depth_out.txt") \
            into alignment_metrics_output

    script:
    """
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${pair_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf
    samtools depth -a ${sorted_dedup_reads} > ${pair_id}_depth_out.txt
    """
}

process haplotypeCaller {
    input:
    set val(pair_id), file(preprocessed_bam) from sorted_dedup_ch

    output:
    set val(pair_id), file("${pair_id}_raw_variants.vcf") into hc_output_ch

    script:
    """
    gatk HaplotypeCaller \
	-R $ref \
	-I $preprocessed_bam \
	-O ${pair_id}_raw_variants.vcf \
	-G StandardAnnotation -G StandardHCAnnotation
    """
}

process selectVariants {
    input:
    set val(pair_id), file(raw_variants) from hc_output_ch

    output:
    set val(pair_id), \
	file("${pair_id}_raw_snps.vcf") \
	into raw_snps_ch, raw_snps_metrics_ch
    set val(pair_id), \
	file("${pair_id}_raw_indels.vcf") \
	into raw_indels_ch

    script:
    """
    gatk SelectVariants \
	-R $ref \
	-V $raw_variants \
	-select-type SNP \
	-O ${pair_id}_raw_snps.vcf
    gatk SelectVariants \
        -R $ref \
        -V $raw_variants \
        -select-type INDEL \
        -O ${pair_id}_raw_indels.vcf
    """
}

process filterSnps {
    publishDir "${params.outdir}/filtered_snps", mode:'copy'

    input:
    set val(pair_id), file(raw_snps) from raw_snps_ch

    output:
    set val(pair_id), \
	file("${pair_id}_filtered_snps.vcf") \
	into filtered_snps_ch, filtered_snps_metrics_ch

    script:
    """
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${pair_id}_filtered_snps.vcf \
	-filter-name "QUAL_filter" -filter "QUAL < 30" \
	-filter-name "GQ_filter" -filter "GQ < 40" \
	-filter-name "AF_filter" -filter "AF < 1.00"
    """
}

process filterIndels {
    publishDir "${params.outdir}/filtered_indels", mode:'copy'
    input:
    set val(pair_id), file(raw_indels) from raw_indels_ch

    output:
    set val(pair_id), \
	file("${pair_id}_filtered_indels.vcf") \
	into filtered_indels_ch

    script:
    """
    gatk VariantFiltration \
        -R $ref \
        -V $raw_indels \
        -O ${pair_id}_filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process parseMetrics {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    set val(pair_id), \
	file("${pair_id}_alignment_metrics.txt"), \
	file("${pair_id}_insert_metrics.txt"), \
	file("${pair_id}_insert_size_histogram.pdf"), \
	file("${pair_id}_depth_out.txt") \
	from alignment_metrics_output
    set val(pair_id), \
	file("${pair_id}_dedup_metrics.txt") \
	from dedup_metrics_ch
    set val(pair_id), \
	file("${pair_id}_raw_snps.vcf") \
	from raw_snps_metrics_ch
    set val(pair_id), \
        file("${pair_id}_filtered_snps.vcf") \
        from filtered_snps_metrics_ch

    output:
    file("report.csv") into parse_metrics_output

    script:
    header = "ID,# reads,aligned reads,% aligned,aligned bases,read length,% paired, %dup, mean insert size,# SNPs,# SNPs filtered,average coverage"
    """
    # If report.csv doesn't exist, create it with the header
    [ -f report.csv ] || echo $header > report.csv

    # run the parse_metrics.sh and append output to report.csv
    parse_metrics.sh ${pair_id} >> report.csv
    """
}
