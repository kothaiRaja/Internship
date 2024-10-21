nextflow.enable.dsl = 2

// Defining all the parameters
params.transcriptome_file = "$launchDir/data/ggal/transcriptome.fa"
params.accession = "SRR1777174"
params.with_fastqc = false
params.with_fastp = false
params.cache = "${launchDir}/cache"
params.out = "${launchDir}/out"

process prefetch {
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
    storeDir params.cache
    input:
        val accession
    output:
        path "${accession}"
    """
    prefetch ${accession}
    """
}

process fasterqdump {
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
    storeDir params.cache
    input:
        path sradir
    output:
        path "${sradir}.fastq"
    """
    fasterq-dump ${sradir} 
    """
} 

process fastqc {
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    input:
        path fastqfile
    output:
        path "fastqc_${fastqfile.getSimpleName()}"
    """
    mkdir fastqc_${fastqfile.getSimpleName()} 
    fastqc -o fastqc_${fastqfile.getSimpleName()} ${fastqfile}
    """
}

process trimReads {
    container "https://depot.galaxyproject.org/singularity/trim_galore%3A0.6.7--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path fastqfile
    output:
        path "${fastqfile.baseName}_trimmed.fastq.gz"
    """
    trim_galore --output_dir. --gzip ${fastqfile}
    """
}

process alignReads {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10a--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path trimmed_fastq
        path genome_index
    output:
        path "${trimmed_fastq.baseName}_aligned.bam"
        path "${trimmed_fastq.baseName}_aligned.bam.bai"
    script:
    """
    STAR --runThreadN 4 --genomeDir ${genome_index} --readFilesIn ${trimmed_fastq} --outFileNamePrefix ${trimmed_fastq.baseName} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15
    samtools index ${trimmed_fastq.baseName}_aligned.bam
    """
}

process markDuplicatesAndSplit {
    container "https://depot.galaxyproject.org/singularity/picard%3A2.27.4--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path bamfile
        path bai_file
    output:
        path "${bamfile.baseName}_dedupped.bam"
        path "${bamfile.baseName}_dedupped.bam.bai"
    script:
    """
    java -jar /usr/local/share/picard-tools/picard.jar AddOrReplaceReadGroups I=${bamfile} O=rg_added_sorted.bam SO=coordinate RGID=ID_NAME RGLB=library RGPL=illumina RGPU=identifier RGSM=sample_name
    java -jar /usr/local/share/picard-tools/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
    java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${params.genome_ref} -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -fixMisencodedQuals
    samtools index split.bam
    mv split.bam ${bamfile.baseName}_dedupped.bam
    mv split.bam.bai ${bamfile.baseName}_dedupped.bam.bai
    """
}

process baseRecalibration {
    container "https://depot.galaxyproject.org/singularity/gatk%3A4.2.6.1--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path bamfile
        path bai_file
    output:
        path "${bamfile.baseName}_recalibrated.bam"
        path "${bamfile.baseName}_recalibrated.bam.bai"
    script:
    """
    gatk --java-options '-Xmx${params.memory}' BaseRecalibrator -R ${params.genome_ref} -I ${bamfile} -knownSites ${params.known_sites} -O recalibration.table
    gatk --java-options '-Xmx${params.memory}' PrintReads -R ${params.genome_ref} -I ${bamfile} -BQSR recalibration.table -O recalibrated.bam
    samtools index recalibrated.bam
    mv recalibrated.bam ${bamfile.baseName}_recalibrated.bam
    mv recalibrated.bam.bai ${bamfile.baseName}_recalibrated.bam.bai
    """
}

process callVariantsDeepVariant {
    container "https://depot.galaxyproject.org/singularity/deepvariant%3A1.2.0--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path bamfile
        path genome_ref
    output:
        path "${bamfile.baseName}_variants.vcf"
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${genome_ref} \
        --reads=${bamfile} \
        --output_vcf=${bamfile.baseName}_variants.vcf \
        --output_gvcf=${bamfile.baseName}_variants.gvcf.gz
    """
}

process filterVariants {
    container "https://depot.galaxyproject.org/singularity/gatk%3A4.2.6.1--hdfd78af_0"
    storeDir params.cache
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path vcf_file
        path genome_ref
    output:
        path "${vcf_file.baseName}_filtered.vcf.gz"
    script:
    """
    gatk --java-options '-Xmx${params.memory}' VariantFiltration \
      -R ${genome_ref} \
      -V ${vcf_file} \
      -filterName QDFilter -filter "QD < 2.0" \
      -filterName FSFilter -filter "FS > 30.0" \
      -o ${vcf_file.baseName}_filtered.vcf.gz
    bgzip -c ${vcf_file.baseName}_filtered.vcf.gz > ${vcf_file.baseName}_filtered.vcf.gz
    tabix -p vcf ${vcf_file.baseName}_filtered.vcf.gz
    """
}

workflow {
    if(params.accession == null) {
        print("Please provide an accession, e.g. '--accession SRR1777174'.")
        System.exit(1)
    }
    accession_channel = Channel.from(params.accession)
    sra_channel = prefetch(accession_channel)
    fastq_channel = fasterqdump(sra_channel)
    quality_channel = fastqc(fastq_channel)
    trimmed_channel = trimReads(fastq_channel)
    aligned_channel = alignReads(trimmed_channel, params.genome_index)
    deduped_channel = markDuplicatesAndSplit(aligned_channel)
    recalibrated_channel = baseRecalibration(deduped_channel)
    variant_calling_channel = callVariantsDeepVariant(recalibrated_channel, params.genome_ref)
    filtered_channel = filterVariants(variant_calling_channel, params.genome_ref)
}
