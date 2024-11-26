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



workflow {
    // Load and parse the CSV file into a channel
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row -> [ row.sample_id, row.read1, row.read2 ] }

    // Download reads
    reads_channel = download_reads(samples_channel)

    // Perform FastQC on the downloaded reads
    fastqc_raw(reads_channel)
}

