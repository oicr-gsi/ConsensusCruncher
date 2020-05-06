version 1.0

workflow consensusCruncher {
  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
  }

  parameter_meta {
    fastqR1: "First Fastq File"
    fastqR2: "Second Fastq File"
    outputFileNamePrefix: "Prefix to use for output file"
  }

  call alignConsensus { input: fastqR1 = fastqR1, fastqR2 = fastqR2, outputFileNamePrefix = outputFileNamePrefix}

  meta {
    author: "Alexander Fortuna"
    email: "alexander.fortuna@oicr.on.ca"
    description: "Workflow to run extract UMI from fastq and generate consensus Bams"
    dependencies: [
     {
      name: "hg19-bwa-index/0.7.12",
      url: "http://bio-bwa.sourceforge.net/"
     },
     {
      name: "samtools/1.9",
      url: "http://www.htslib.org/"
     },
     {
      name: "python/3.6",
      url: "https://www.python.org/downloads/"
     },
     {
      name: "picard/2.21.2",
      url: "https://broadinstitute.github.io/picard/"
     },
     {
      name: "rstats/3.6",
      url: "https://www.r-project.org/"
     },
     {
      name: "consensuscruncer-5.0",
      url: "https://github.com/pughlab/ConsensusCruncherb"
     }
    ]
  }

  output {
    File dcsScBam = alignConsensus.dcsScBam
    File dcsScBamIndex = alignConsensus.dcsScBamIndex
    File allUniqueBam = alignConsensus.allUniqueBam
    File allUniqueBamIndex = alignConsensus.allUniqueBamIndex
    File sscsScBam = alignConsensus.sscsScBam
    File sscsScBamIndex = alignConsensus.sscsScBamIndex
    File ccFolder = alignConsensus.ccFolder

  }
}

task alignConsensus {
  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
    String consensusCruncher = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String modules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12"
    String bwa = "$BWA_ROOT/bin/bwa"
    String bwaref = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String blist = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"
    String genome   = "hg19"
    String cytoband = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt"
    String name = "_R"
    Float cutoff  = 0.7
    Int threads = 4
    Int jobMemory = 64
    Int timeout = 72
  }

  parameter_meta {
    fastqR1: "path to left fastq file"
    fastqR2: "path to right fastq file"
    outputFileNamePrefix: "file name prefix"
    consensusCruncher: "path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
    bwa: "path to bwa binary"
    bwaref: "path to bwa index"
    samtools: "path to samtools binary"
    blist: "path to blacklist for barcodes"
    name: "file name taken to the left of this string"
    threads: "number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    genome: "which genome version to use"
    cytoband: "path to cytoband for genome"
    cutoff: "cutoff to use to call a consenus of reads"
  }

  command <<<

  set -euo pipefail

  mkdir consensuscruncher

  mkdir outputBam

  ~{consensusCruncher} fastq2bam \
         --fastq1 ~{fastqR1} \
         --fastq2 ~{fastqR2}\
         --output 5_consensuscruncher \
         --name "~{name}" \
         --bwa ~{bwa} \
         --ref ~{bwaref} \
         --samtools ~{samtools} \
         --skipcheck \
         --blist ~{blist}

   ~{consensusCruncher} consensus \
         --input 5_consensuscruncher/bamfiles/~{outputFileNamePrefix}.sorted.bam \
         --output 5_consensuscruncher/consensus \
         --samtools ~{samtools} \
         --cutoff ~{cutoff} \
         --genome ~{genome} \
         --bedfile ~{cytoband} \
         --bdelim '|'

   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/dcs_sc/~{outputFileNamePrefix}.sorted.dcs.sc.sorted.bam outputBam/
   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/dcs_sc/~{outputFileNamePrefix}.sorted.dcs.sc.sorted.bam.bai outputBam/
   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/~{outputFileNamePrefix}.sorted.all.unique.dcs.sorted.bam outputBam/
   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/~{outputFileNamePrefix}.sorted.all.unique.dcs.sorted.bam.bai outputBam/
   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/sscs_sc/~{outputFileNamePrefix}.sorted.sscs.sc.sorted.bam outputBam/
   mv consensuscruncher/consensus/~{outputFileNamePrefix}.sorted/sscs_sc/~{outputFileNamePrefix}.sorted.sscs.sc.sorted.bam.bai outputBam/

   gzip consensuscruncher/

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
      File dcsScBam = "outputBam/~{outputFileNamePrefix}.sorted.dcs.sc.sorted.bam"
      File dcsScBamIndex = "outputBam/~{outputFileNamePrefix}.sorted.dcs.sc.sorted.bam.bai"
      File allUniqueBam = "outputBam/~{outputFileNamePrefix}.sorted.all.unique.dcs.sorted.bam"
      File allUniqueBamIndex = "outputBam/~{outputFileNamePrefix}.sorted.all.unique.dcs.sorted.bam.bai"
      File sscsScBam = "outputBam/~{outputFileNamePrefix}.sorted.sscs.sc.sorted.bam"
      File sscsScBamIndex = "outputBam/~{outputFileNamePrefix}.sorted.sscs.sc.sorted.bam.bai"
      File ccFolder = "consensuscruncher.gz"
  }

  meta {
    output_meta: {
      dcsScBam: "DCS generated from SSCS + SC",
      dcsScBamIndex: "Index for DCS SC Bam",
      allUniqueBam: "DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons",
      allUniqueBamIndex: "Index for All Unique Bam",
      sscsScBam: "SSCS combined with corrected singletons (from both rescue strategies)",
      sscsScBamIndex: "Index for SSCS SC Bam",
      ccFolder: "output folder containing files not needed for downstream analysis; info on family size, QC metrics"
    }
  }

}
