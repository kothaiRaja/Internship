nextflow.enable.dsl=2



process STAR_ALIGNMENT {
    tag { sample_id }

   container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/STAR", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path star_index_dir
    

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    # ngs-nf-dev Align reads to genome
  STAR --genomeDir $star_index_dir \
       --readFilesIn ${trimmed_r1} ${trimmed_r2}  \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)  
  STAR --genomeDir $star_index_dir \
       --readFilesIn ${trimmed_r1} ${trimmed_r2} \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
	   --outFileNamePrefix ${sample_id}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$sample_id LB:library PL:illumina PU:machine SM:GM12878
    """
}

process SAMTOOLS_SORT_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/sorted_bam", mode: "copy"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam ${bam}

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam
    """
}

process GATK_MARK_DUPLICATES {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/mark_duplicates", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(sample_id), 
          path("${sample_id}_marked_duplicates.bam"), 
          path("${sample_id}_marked_duplicates.bai"), 
          path("${sample_id}_dup_metrics.txt")

    script:
    """
    gatk MarkDuplicates \
        -I ${sorted_bam} \
        -O ${sample_id}_marked_duplicates.bam \
        -M ${sample_id}_dup_metrics.txt \
        --CREATE_INDEX true
    """
}

process SPLIT_NCIGAR_READS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/split_ncigar", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(genome_fasta), path (index), path (genome_dict)

    output:
    tuple val(sample_id), 
          path("${sample_id}_split.bam"), 
          path("${sample_id}_split.bai")

    script:
    """
    gatk SplitNCigarReads \
        -R ${genome_fasta} \
        -I ${bam} \
        -O ${sample_id}_split.bam \
        --create-output-bam-index true
    """
}

process GATK_RECALIBRATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org-singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/recalibrated_bams", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(genome_fasta), path(index), path(dict), path(known_variants), path(known_variants_index)

    output:
    tuple val(sample_id), 
          path("${sample_id}_recalibrated.bam"), 
          path("${sample_id}_recalibrated.bai"), 
          path("${sample_id}_recal_data.table")

    script:
    """
    # Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -R ${genome_fasta} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal_data.table

    # Step 2: ApplyBQSR
    gatk ApplyBQSR \
        -R ${genome_fasta} \
        -I ${bam} \
        --bqsr-recal-file ${sample_id}_recal_data.table \
        -O ${sample_id}_recalibrated.bam
    """
}

process BED_TO_INTERVAL_LIST {
    tag { bed_file.name }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/intervals", mode: "copy"

    input:
    path bed_file
    path genome_fasta
    path genome_dict

    output:
    path("${bed_file.baseName}.interval_list")

    script:
    """
    gatk BedToIntervalList \
        -I ${bed_file} \
        -O ${bed_file.baseName}.interval_list \
        -SD ${genome_dict}
    """
}

process SCATTER_INTERVAL_LIST {
    tag "Scatter interval list"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/scattered_intervals", mode: 'copy'

    input:
    path interval_list
    path genome_dict

    output:
    path "*.interval_list"

    script:
    """
    mkdir -p scattered_intervals
    gatk IntervalListTools \
        --INPUT ${interval_list} \
        --OUTPUT scattered_intervals \
        --SCATTER_COUNT ${params.scatter_count} \
        --UNIQUE true

    # Move and rename the output files to the working directory with .interval_list extension
    for f in scattered_intervals/*/*; do
        mv "\$f" "\$(dirname \$f)/\$(basename \$f).interval_list"
    done

    mv scattered_intervals/*/*.interval_list .
    rm -r scattered_intervals
    """
}


process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \
        -R ${genome} \
        -I ${bam} \
        -O ${sample_id}.g.vcf.gz \
        -L ${interval_list} \
        -ERC GVCF
    """
}





workflow {
    

    // Load trimmed reads
    trimmed_reads_ch = Channel.fromFilePairs(params.reads, flat: true)
        
    // Run STAR Alignment
    aligned_bams = STAR_ALIGNMENT(trimmed_reads_ch, params.star_index_dir)
	
	// Sort and index BAM files
    sorted_bams = SAMTOOLS_SORT_INDEX(aligned_bams)
	
	// Mark duplicates
    marked_bams = GATK_MARK_DUPLICATES(sorted_bams)
	
	// Split N CIGAR reads
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict) })
	
	// Recalibrate and Apply BQSR in one step
    recalibrated_bams = GATK_RECALIBRATION(
        split_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict, params.filtered_vcf, params.filtered_vcf_index) })
		
	// Convert BED to interval list
    interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist, params.genome, params.genome_dict)
	
	// Scatter the Interval List
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.genome_dict)
	
    //GATK HaplotypeCaller
	gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams, params.genome,params.fasta_index,params.genome_dict,  scattered_intervals_ch)


}
