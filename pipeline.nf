nextflow.enable.dsl = 2

import nextflow.Channel

// Processes
process download_reads {
	publishDir "data/reads", mode: "copy"
    tag { sample_id }

    input:
    tuple val(sample_id), val(url1), val(url2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fq.gz"), path("${sample_id}_R2.fq.gz")

    script:
    """
    wget -O ${sample_id}_R1.fq.gz ${url1}
    wget -O ${sample_id}_R2.fq.gz ${url2}
    """
}

process fastqc_raw {
    tag { sample_id }
	publishDir "${params.outdir}/fastqc/raw", mode: "copy"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")


    script:
    """
    fastqc ${r1} ${r2} --outdir .
    """
}

process trim_reads {
    tag { sample_id }
	container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_3"
    publishDir "data/trimmed_input", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), emit: trimmed_reads
	tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    

    script:
    """
    fastp -i ${r1} -I ${r2} \
          -o trimmed_${sample_id}_R1.fastq.gz \
          -O trimmed_${sample_id}_R2.fastq.gz \
          --detect_adapter_for_pe \
          --html ${sample_id}_fastp.html \
          --json ${sample_id}_fastp.json
    """
}

//Create a fasta Genome Index

process create_fasta_index {
	container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "data/", mode: "copy"
	
    input:
    path genome_fasta

    output:
    path("${genome_fasta}.fai")

    script:
    """
    samtools faidx ${genome_fasta}
    """
}

// Process: Create Genome Dictionary

process create_genome_dict {
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "data/", mode: "copy"
    input:
    path genome_fasta

    output:
    path("${genome_fasta.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary -R $genome_fasta -O ${genome_fasta.baseName}.dict
    """
}

//Process: STAR Genome Index

process create_star_index {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "data/", mode: "copy"

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR_index", type: 'dir'

    script:
    """
    mkdir -p STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir STAR_index \
         --genomeFastaFiles ${genome_fasta} \
         --sjdbGTFfile ${genome_gtf} \
         --runThreadN 8
    """
}






nextflow.enable.dsl = 2

workflow {
    // Step 1: Load and parse sample metadata
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.read1, row.read2) }

    // Step 2: Download raw reads
    reads_channel = download_reads(samples_channel)

    // Step 3: Perform FastQC on raw reads
    fastqc_raw_results = fastqc_raw(reads_channel)

    // Step 4: Trim reads with fastp
    trimmed_reads = trim_reads(reads_channel)

    // Step 6: Load genome files
    genome_fasta = Channel.value(file(params.genome_fasta))
    genome_gtf = Channel.value(file(params.genome_gtf))

    // Step 8: Create genome index files
    fasta_index = create_fasta_index(Channel.value(file(params.genome_fasta)))
    genome_dict = create_genome_dict(Channel.value(file(params.genome_fasta)))
    star_index = create_star_index(Channel.value(file(params.genome_fasta)), Channel.value(file(params.genome_gtf)))

}
