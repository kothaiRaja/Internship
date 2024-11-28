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
}
