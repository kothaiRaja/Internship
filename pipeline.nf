nextflow.enable.dsl = 2



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
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), emit: trimmed_reads
	

    

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
    publishDir "${params.outdir}/genome/index", mode: "copy"
	
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
    publishDir "${params.outdir}/genome/dict", mode: "copy"
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
    publishDir "${params.outdir}/genome/starindex", mode: "copy"
    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR_index"

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

process star_mapping {
    tag { sample_id }
	container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/star", mode: "copy"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR --runThreadN 8 \
         --genomeDir ${star_index} \
         --readFilesIn ${r1} ${r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${sample_id}.
    """
}



workflow {
    // Step 1: Load and parse the CSV file into a channel
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.read1, row.read2) }
    samples_channel.view() // Debug: View the parsed samples

    // Step 2: Download reads
    reads_channel = download_reads(samples_channel)
    reads_channel.view() // Debug: View the downloaded reads

    // Step 3: Perform trimming with Fastp
	fastp_results = trim_reads(reads_channel)

	// Reshape the output for downstream processes
	trimmed_reads = fastp_results.map { sample_id, r1, r2 ->
    tuple(sample_id, r1, r2) // Separate r1 and r2 for STAR mapping
	}
	trimmed_reads.view() // Debug: Verify reshaped output


    // Step 5: Input genome files
    genome_channel = Channel.fromPath(params.genome_fasta)
    gtf_channel = Channel.fromPath(params.genome_gtf)
    genome_channel.view() // Debug: View the genome FASTA
    gtf_channel.view() // Debug: View the GTF file

    // Step 6: Create FASTA index
    fasta_index = create_fasta_index(genome_channel)

    // Step 7: Create genome dictionary
    genome_dict = create_genome_dict(genome_channel)

    // Step 8: Generate STAR genome index
    star_index = create_star_index(genome_channel, gtf_channel)
    star_index.view() // Debug: View the STAR index directory

    // Step 9: Perform STAR mapping
    aligned_reads = star_mapping(trimmed_reads, star_index)
    aligned_reads.view() // Debug: View the STAR mapping outputs
}

