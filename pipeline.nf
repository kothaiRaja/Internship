nextflow.enable.dsl = 2

//defining all the parameters

//params.reads = "$launchDir/data/ggal/gut_{1,2}.fq"
//params.transcriptome_file = "$launchDir/data/ggal/transcriptome.fa"
//params.multiqc = "$launchDir/multiqc"
params.accession= "SRR16641606"
params.cache = "$launchDir/output"

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

workflow {
  if(params.accession == null) {
    print("Please provide an accession, e.g. '--accession SRR1777174'.")
    System.exit(1)
  }
  accession_channel = Channel.from(params.accession)
  sra_channel = prefetch(accession_channel)
}
