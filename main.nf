nextflow.enable.dsl=2

params.input_dir = "./trimmed_input"  // Directory containing trimmed reads
params.output_dir = "./test_output"  // Output directory for test files

workflow {
    Channel
        .fromFilePairs("${params.input_dir}/*_R{1,2}.fastq.gz", flat: true)
        .set { sample_pairs }

    sample_pairs | testProcess
}

process testProcess {
    input:
    tuple val(sample_id), path(reads)

    output:
    path("${params.output_dir}/")

    script:
    """
    mkdir -p ${params.output_dir}
    echo "Processing sample ${sample_id}" > ${params.output_dir}/${sample_id}.txt
    """
}
