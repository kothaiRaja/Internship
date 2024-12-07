nextflow.enable.dsl=2



process STAR_ALIGNMENT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/STAR/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path star_index_dir

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    # First Pass
    STAR --genomeDir $star_index_dir \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN $task.cpus \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMattributes NH HI AS nM MD NM \
         --outFileNamePrefix ${sample_id}_pass1_

    # Second Pass
    STAR --genomeDir $star_index_dir \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN $task.cpus \
         --readFilesCommand zcat \
         --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMattributes NH HI AS nM MD NM \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 12000000000 \
         --outSAMattrRGline ID:$sample_id LB:library PL:illumina PU:machine SM:$sample_id \
         --outFileNamePrefix ${sample_id}_
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

process SAMTOOLS_FLAGSTAT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt")

    script:
    """
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt
    """
}


process GATK_MARK_DUPLICATES {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
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

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
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

process SAMTOOLS_CALMD {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), 
          path(bam), 
          path(bai) 
          path(genome_fasta) 
		  path (index)
		  
    output:
    tuple val(sample_id), 
          path("${sample_id}_calmd.bam"), 
          path("${sample_id}_calmd.bam.bai") 
         

    script:
    """
    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${sample_id}_calmd.bam

    # Index the updated BAM file
    samtools index ${sample_id}_calmd.bam

    """
}


process GATK_RECALIBRATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
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

	# Step 3: Validation
	gatk ValidateSamFile \
		-I ${sample_id}_recalibrated.bam \
		-MODE SUMMARY 


    """
}

process BED_TO_INTERVAL_LIST {
    tag { bed_file.name }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
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

     container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
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

     container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(table)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi") 

    script:
    """
    gatk HaplotypeCaller \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${genome} \
    --output output_${sample_id}.vcf.gz \
    -I $bam \
    --standard-min-confidence-threshold-for-calling 5.0 \
    --dont-use-soft-clipped-bases true \
    --min-base-quality-score 10 \
	--output-mode EMIT_VARIANTS_ONLY \
    --intervals $interval_list

    """
}

process GATK_MERGE_VCFS {
    tag "Combine GVCFs"
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/combined_gvcf", mode: 'copy'

    input:
    path genome
	path genome_index
	path genome_dict
    tuple val(sample_id), path (vcfs), path(vcf_index)  // List of GVCFs

    output:
    path "merged_output.vcf.gz", emit: vcf          // Merged VCF file
    path "merged_output.vcf.gz.tbi", emit: vcf_idx  // Index for the merged VCF file

    script:
    """
    gatk MergeVcfs \
        ${vcfs.collect { "-I $it" }.join(' ')} \
        --OUTPUT merged_output.vcf.gz
    """
}

process BCFTOOLS_STATS {
    tag "bcftools_stats"

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h8b25389_0" // Adjust with your container or module if needed
    publishDir "${params.outdir}/bcftools_stats/beforefilteration", mode: "copy"

    input:
    path(vcf_file)    // Input VCF file
    path(vcf_index)   // Input VCF index (.tbi) file

    output:
    path("stats.txt")     // Output stats file

    script:
    """
    # Generate stats
    bcftools stats ${vcf_file} > stats.txt

    """
}


process GATK_VARIANT_FILTER {
    tag "variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter", mode: "copy"

    input:
    path(vcf_file)
	path(vcf_index)  // Input VCF and index
    path genome 	// Reference genome
	path genome_index
	path genome_dict

    output:
          path("final.vcf.gz") 
          path("final.vcf.gz.tbi")  // Filtered VCF and index

    script:
    """
    gatk VariantFiltration \
    -R ${genome} \
    -V ${vcf_file} \
    --cluster-window-size 35 --cluster-size 3 \
    --filter-name "LowQual" --filter-expression "QUAL < 30.0" \
    --filter-name "LowQD" --filter-expression "QD < 2.0" \
    --filter-name "HighFS" --filter-expression "FS > 60.0" \
    --filter-name "LowMQ" --filter-expression "MQ < 40.0" \
    --filter-name "HighSOR" --filter-expression "SOR > 3.0" \
    --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < -2.0" \
    -O final.vcf.gz

    """
}


process ANNOTATE_VARIANTS {
    tag "Annotate variants"
    
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf	// Filtered VCF file
	path index
    path snpEffJar  // Path to the SnpEff JAR file
    path snpEffConfig  // Path to the SnpEff configuration file
    path snpEffDbDir	// Path to the SnpEff database directory
	val genomedb

    output:
    path "annotated.vcf"
	path "annotated.summary.html"

    script:
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated.vcf

    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated.summary.html \
        ${vcf} > /dev/null
    """
}

process MULTIQC_REPORT {
    tag "Generate MultiQC report"

    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.20--pyhdfd78af_0" 
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path results_dir 

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${results_dir} -o .
    """
}



workflow {
    

    // Load trimmed reads
    trimmed_reads_ch = Channel.fromFilePairs(params.reads, flat: true)
        
    // Run STAR Alignment
    aligned_bams = STAR_ALIGNMENT(trimmed_reads_ch, params.star_index_dir)
	
	// Sort and index BAM files
    sorted_bams = SAMTOOLS_SORT_INDEX(aligned_bams)
	
	// Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(sorted_bams)
	
	// Mark duplicates
    marked_bams = GATK_MARK_DUPLICATES(sorted_bams)
	
	// Split N CIGAR reads
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict) })
	
	//SAMTOOLS CALMD process
	calmd_ch = SAMTOOLS_CALMD(split_bams , params.genome, params.fasta_index )
	
	// Recalibrate and Apply BQSR in one step
    recalibrated_bams = GATK_RECALIBRATION(
        calmd_ch.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict, params.filtered_vcf, params.filtered_vcf_index) })
		
	// Convert BED to interval list
    interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist, params.genome, params.genome_dict)
	
	// Scatter the Interval List
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.genome_dict)
	
	//GATK HaplotypeCaller
	gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams, params.genome,params.fasta_index,params.genome_dict, scattered_intervals_ch)
	

	// Combine GVCFs
	merged_vcf = GATK_MERGE_VCFS(params.genome,params.fasta_index, params.genome_dict,gvcf_output)
	
	//provide stats
	bcftools_stats_ch = BCFTOOLS_STATS(merged_vcf)

	// Variant Filtering
    filtered_vcf = GATK_VARIANT_FILTER(merged_vcf,params.genome,params.fasta_index, params.genome_dict )
		
	// Pass all required inputs to ANNOTATE_VARIANTS
    annotated_vcf = ANNOTATE_VARIANTS(filtered_vcf, file('./data/test/snpEff/snpEff.jar'),
        file('./data/test/snpEff/snpEff.config'),
        file('./data/test/snpEff/snpEff/data') , params.genomedb)   
		
	//multiqc
	multiqc_results = MULTIQC_REPORT(Channel.fromPath("${params.outdir}"))

}
