nextflow.enable.dsl = 2

// Process to download the test genome
process DOWNLOAD_TEST_GENOME {
    tag "Download test genome"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "genome.fa"

    script:
    """
    wget -q -O genome.fa ${params.test_data_url}/genome.fa
    """
}

// Process to download the test variants VCF
process DOWNLOAD_TEST_VARIANTS {
    tag "Download test variants VCF"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "variants.vcf"

    script:
    """
    wget -q -O variants.vcf ${params.test_data_url}/subset_chr22.vcf.gz
    """
}

// Process to download the test denylist BED
process DOWNLOAD_TEST_DENYLIST {
    tag "Download test denylist BED"
	container null	
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "denylist.bed"

    script:
    """
    wget -q -O denylist.bed ${params.test_data_url}/denylist_chr22_to_22.bed
    """
}

// Process to download gtf files
process DOWNLOAD_TEST_GTF {
    tag "Download test GTF"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "annotations.gtf"

    script:
    """
    wget -q -O annotations.gtf ${params.test_data_url}/annotations.gtf
    """
}


process CHECK_JAVA {
    tag "Check Java"
	container null
    output:
    path "java_check.txt"

    script:
    """
    if ! java -version &>/dev/null; then
        echo "Java is not installed or not in PATH." > java_check.txt
        exit 1
    else
        echo "Java is available." > java_check.txt
    fi
    """
}


process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
	publishDir "${params.test_data_dir}", mode: 'copy'
	container null
    output:
    path "${params.snpeff_jar_dir}/snpEff.jar"
	path "${params.snpeff_jar_dir}/snpEff.config"

    script:
    """
    mkdir -p ${params.snpeff_jar_dir}
    wget -q -O snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -j snpEff_latest_core.zip -d ${params.snpeff_jar_dir}
    rm snpEff_latest_core.zip
    """
}

process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
	publishDir "${params.test_data_dir}/snpEff", mode: 'copy'
	container null
    input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"

    script:
    """
    # Ensure the output directory exists first
    mkdir -p ${params.snpeff_db_dir}

    # Use an absolute path for the data directory
    data_dir=\$(realpath ${params.snpeff_db_dir})

    # Download the database
    java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -v
    """
}

// Processes
process download_reads {
	publishDir "${params.test_data_dir}", mode: "copy"
	container null
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
    publishDir "${params.outdir}/fastp", mode: "copy"
	storeDir "${params.test_data_dir}"
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
    publishDir "${params.test_data_dir}", mode: 'copy'
	
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
    publishDir "${params.test_data_dir}", mode: 'copy'
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
    publishDir "${params.test_data_dir}", mode: 'copy'

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




workflow {
    // Step 1: Download reference files
    def genome = DOWNLOAD_TEST_GENOME()
    def variants = DOWNLOAD_TEST_VARIANTS()
    def denylist = DOWNLOAD_TEST_DENYLIST()
    def genome_gtf = DOWNLOAD_TEST_GTF()

    // Step 2: Check Java installation
    CHECK_JAVA()

    // Step 3: Download and configure SnpEff
    def snpeff_jar_and_config = DOWNLOAD_SNPEFF_TOOL()
    def snpeff_db = DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_and_config[0])

    // Step 4: Load and parse sample metadata
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.read1, row.read2) }

    // Step 5: Download raw reads
    reads_channel = download_reads(samples_channel)

    // Step 6: Perform FastQC on raw reads (depends on downloaded reads)
    fastqc_raw_results = fastqc_raw(reads_channel)

    // Step 7: Trim reads with fastp (depends on FastQC)
    trimmed_reads = trim_reads(reads_channel)

    // Step 8: Create genome index files
    // Fasta index depends on genome download
    fasta_index = create_fasta_index(genome)
    // Genome dictionary depends on genome download
    genome_dict = create_genome_dict(genome)
    // STAR index depends on genome and GTF downloads
    star_index = create_star_index(genome, genome_gtf)
}
