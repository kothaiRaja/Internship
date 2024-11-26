nextflow.enable.dsl = 2



process parse_csv {
	publishDir "data/reads", mode: "copy"
    input:
    path csv_file // Accepts the CSV file as input

    output:
    path "samples.txt" // Produces a text file listing samples

    script:
    """
    # Extract sample names and URLs, removing empty lines
    tail -n +2 $csv_file | tr -d '\r' | while IFS=, read -r sample_name r1_url r2_url
    do
        # Only include lines with all fields populated
        if [[ -n "\$sample_name" && -n "\$r1_url" && -n "\$r2_url" ]]; then
            echo "\$sample_name \$r1_url \$r2_url" >> samples.txt
        fi
    done
    """
}



process download_reads {
    
	publishDir "data/reads", mode: "copy", saveAs: { filename -> filename }
    input:
    path sample_list // Takes samples.txt as input

    output:
    path "*.fastq.gz"

    script:
    """
    tr -d '\r' < $sample_list | while read sample_name r1_url r2_url; do
        # Skip empty lines
        if [[ -z "\$sample_name" || -z "\$r1_url" || -z "\$r2_url" ]]; then
            continue
        fi
        echo "Processing Sample: \$sample_name"
        wget -O \${sample_name}_R1.fastq.gz \$r1_url
        wget -O \${sample_name}_R2.fastq.gz \$r2_url
    done
    """
}





process fastqc_raw {
    publishDir "${params.outdir}/fastqc/raw", mode: "copy"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"

    input:
    path fastq_file

    output:
    path "*.html" // FastQC HTML reports
    path "*.zip"  // FastQC ZIP files

    script:
    """
    echo "Running FastQC for: ${fastq_file}"
    fastqc --outdir . ${fastq_file}
    """
}








workflow {
    // Create a channel for the CSV file
    csv_channel = Channel.fromPath(params.csv_file)

    // Pass the CSV file to parse_csv process
    parsed_output = parse_csv(csv_channel)

    // Pass the output of parse_csv to download_reads process
    downloaded_reads = download_reads(parsed_output)

    // Perform FastQC on each downloaded read
    fastqc_raw(downloaded_reads)
}




