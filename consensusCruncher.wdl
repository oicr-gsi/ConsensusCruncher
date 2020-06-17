version 1.0

workflow consensusCruncher {
  input {
    File? fastqR1
    File? fastqR2
    File? sortedBam
    File? sortedBai
    String outputFileNamePrefix
  }

  parameter_meta {
    fastqR1: "First Fastq File"
    fastqR2: "Second Fastq File"
    sortedBam: "Bam file from bwamem"
    sortedBai: "Bai file from bwamem"
    outputFileNamePrefix: "Prefix to use for output file"
  }

  if (!(defined(sortedBam)) && defined(fastqR1) && defined(fastqR2)) {
    call align {
      input:
        fastqR1 = select_first([fastqR1]),
        fastqR2 = select_first([fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  call consensus {
    input:
      inputBam = select_first([sortedBam, align.sortedBam]),
      inputBai = select_first([sortedBai, align.sortedBai]),
      basePrefix = outputFileNamePrefix
  }

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
      url: "https://github.com/pughlab/ConsensusCruncher"
     }
    ]
  }

  output {
    File dcsScBam = consensus.dcsScBam
    File dcsScBamIndex = consensus.dcsScBamIndex
    File allUniqueBam = consensus.allUniqueBam
    File allUniqueBamIndex = consensus.allUniqueBamIndex
    File sscsScBam = consensus.sscsScBam
    File sscsScBamIndex = consensus.sscsScBamIndex
    File ccFolder = consensus.ccFolder
  }
}

task align {
  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
    String modules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String bwa = "$BWA_ROOT/bin/bwa"
    String bwaref = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String blist = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"
    Int threads = 4
    Int jobMemory = 16
    Int timeout = 72
  }

  parameter_meta {
    fastqR1: "Path to left fastq file"
    fastqR2: "Path to right fastq file"
    outputFileNamePrefix: "File name prefix"
    consensusCruncherPy: "Path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
    bwa: "Path to bwa binary"
    bwaref: "Path to bwa index"
    samtools: "Path to samtools binary"
    blist: "Path to blacklist for barcodes"
    name: "File name taken to the left of this string"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    ~{consensusCruncherPy} fastq2bam \
         --fastq1 ~{fastqR1} \
         --fastq2 ~{fastqR2}\
         --output . \
         --bwa ~{bwa} \
         --ref ~{bwaref} \
         --samtools ~{samtools} \
         --skipcheck \
         --blist ~{blist}

    # Necessary for bam files to be named according to merged library name
    # Additionally if ".sorted" isn't omitted here, file names from align include ".sorted" twice
    mv bamfiles/*.bam bamfiles/"~{outputFileNamePrefix}.bam"
    mv bamfiles/*.bai bamfiles/"~{outputFileNamePrefix}.bam.bai"
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File? sortedBam = "bamfiles/~{outputFileNamePrefix}.bam"
    File? sortedBai = "bamfiles/~{outputFileNamePrefix}.bam.bai"
  }
}

task consensus {
  input {
    File inputBam
    File inputBai
    String basePrefix
    String modules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 samtools/1.9"
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String cytoband = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt"
    String genome   = "hg19"
    Float cutoff  = 0.7
    Int threads = 8
    Int jobMemory = 16
    Int timeout = 72
  }

  parameter_meta {
    inputBam: "Bam file either from bwamem or ConsensusCruncher align."
    inputBai: "Bai file either from bwamem or ConsensusCruncher align."
    consensusCruncherPy: "Path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
    samtools: "Path to samtools binary"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    genome: "Which genome version to use"
    cytoband: "Path to cytoband for genome"
    cutoff: "Cutoff to use to call a consenus of reads"
  }

  String ccDir = basePrefix + ".consensuscruncher"

  command <<<
  set -euo pipefail

   ~{consensusCruncherPy} consensus \
         --input ~{inputBam} \
         --output . \
         --samtools ~{samtools} \
         --cutoff ~{cutoff} \
         --genome ~{genome} \
         --bedfile ~{cytoband} \
         --bdelim '|'

   tar cf - ~{basePrefix} | gzip --no-name > ~{ccDir}.tar.gz
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File dcsScBam = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam"
    File dcsScBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam.bai"
    File allUniqueBam = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam"
    File allUniqueBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam.bai"
    File sscsScBam = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam"
    File sscsScBamIndex = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam.bai"
    File ccFolder = "~{ccDir}.tar.gz"
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
