#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/nascent
========================================================================================
 nf-core/nascent Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/nascent
----------------------------------------------------------------------------------------
 #### Authors
 Ignacio Tripodi <ignacio.tripodi@colorado.edu>
 Margaret Gruca <margaret.gruca@colorado.edu>
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        1a. SRA tools -- fastq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files
        1c. Gzip fastq -- compress fastq files for storage

    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)

    3. Mapping w/ HISAT2 -- map to genome reference file

    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM

    5. Quality control
        5a. preseq -- estimate library complexity
        5b. RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        5c. Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio

    6. Coverage files
        6a. deepTools : normalized bigwigs
        6b. BEDTools and kentUtils : 5' bigwigs for dREG
        6c. deepTools : normalized bedgraphs
        6d. BEDTools : non-normalized bedgraphs

    7. IGV Tools : bedGraph --> tdf

    8. MultiQC : generate QC report for pipeline

    9. Pipeline report

=======
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/nascent -profile slurm --reads '/project/*_{R1,R2}*.fastq' --outdir '/project/'
    nextflow run nf-core/nascent --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Required arguments:
         -profile                      Configuration profile to use. <base, fiji>
         --reads                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --reads not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --saveAll                      Compresses and saves all fastq reads.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --saveReference               Save the generated reference files the the Results directory.

    QC Options:
        --skipMultiQC                  Skip running MultiQC report.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/assets/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.workdir = "./nextflowTemp"


// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// Validate inputs

if ( params.fasta ){
   // genome_fasta = file(params.fasta)
   // if( !genome_fasta.exists() ) exit 1, "Genome directory not found: ${params.fasta}"
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
           .into { genome_fasta; ch_fasta_for_hisat_index}
}
else {
    params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
}

if( params.chrom_sizes ){
    Channel
        .fromPath(params.chrom_sizes, checkIfExists: true)
        .ifEmpty { exit 1, "Chrom sizes file not found: ${params.chrom_sizes}" }
        .into { chrom_sizes_for_bed;
                chrom_sizes_for_bigwig;
                chrom_sizes_for_igv }
}
else {
    params.chrom_sizes = null
}

if ( params.bbmap_adapters ){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}
else {
    hisat2_indices = null
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}
else {
    genome_refseq = null
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
 // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

/*
 * Create a channel for input read files
 */
if (params.reads) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.baseName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.baseName, file) }
        fastq_reads_gzip = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.baseName, file) }
    } else {
        Channel
            .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    params.reads = null
}

if (params.sras) {
    if (params.singleEnd) {
    println("Pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
    } else {
         Channel
             .fromFilePairs( params.sras, size: params.singleEnd ? 1 : 2 )
             .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
             .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    }
}

else {
    read_files_sra = Channel.empty()
 }


/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_files_fastqc; read_files_trimming }
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
if(params.reads) summary['Fastqs']           = params.reads
if(params.sras) summary['SRAs']             = params.sras
summary['Genome Ref']       = params.fasta
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-nascent-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/nascent Workflow Summary'
    section_href: 'https://github.com/nf-core/nascent'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    preseq > v_preseq.txt
    seqkit version > v_seqkit.txt
    bedtools --version > v_bedtools.txt
    export LC_ALL=C
    igvtools version > v_igv-tools.txt

    # Can't call this before running MultiQC or it breaks it
    read_distribution.py --version > v_rseqc.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1a -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"

    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_gzip_sra

    script:
    prefix = reads.baseName
    if (!params.threadfqdump) {
        """
        echo ${prefix}

        fastq-dump ${reads}
        """
    } else if (!params.singleEnd) {
         """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads ${task.cpus} \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}

        fastq-dump --split-3 ${reads}
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads ${task.cpus} \
            --sra-id ${reads}
        """
    }
}


/*
 * PREPROCESSING - Build HISAT2 index (borrowed from nf-core/rnaseq)
 */
// TODO: do we need --ss and --exon? probably not, need to check what was the actual hisat2-builder arguments used to generate the indices we have on fiji
if(!params.hisat2_indices && params.fasta){
    process make_hisat_index {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_hisat_index

        output:
        file "*.ht2" into hisat2_indices

        script:
        if( !task.memory ){
            log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
            avail_mem = 0
        } else {
            log.info "[HISAT2 index build] Available memory: ${task.memory}"
            avail_mem = task.memory.toGiga()
        }
        """
        hisat2-build -p ${task.cpus} ${fasta} ${fasta.baseName}-hisat2_index
        """
    }
}


/*
 * STEP 1b - FastQC
 */

process fastqc {
    tag "$prefix"
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    prefix = reads.baseName
    """
#    echo `which gunzip`

    fastqc $reads
		#extract_fastqc_stats.sh --srr=${prefix} > ${prefix}_stats_fastqc.txt
#    GC=\$(gunzip -c "\$(find . -name *_fastqc.zip)" "${prefix}"_fastqc/fastqc_data.txt \
#            | grep "%GC" | grep -o "[0-9]*")
#    SEQ=\$(gunzip -c "\$(find . -name *_fastqc.zip)" "${prefix}"_fastqc/fastqc_data.txt | \
#              grep "Total Sequences" | \
#              grep -o "[0-9]*")
#    DEDUP=\$(gunzip -c "\$(find . -name *_fastqc.zip)" "${prefix}"_fastqc/fastqc_data.txt | \
#                grep "#Total Deduplicated Percentage" | \
#                grep -o "[0-9,.]*")
#
#    echo -e "SRR\t%GC\tTotal_Sequences\t%Total_Deduplicated" > ${prefix}_stats_fastqc.txt
#    echo -e "${prefix}""\$(printf "\\t")""\$GC""\$(printf "\\t")""\$SEQ""\$(printf "\\t")""\$DEDUP" >> ${prefix}_stats_fastqc.txt
    """
}


/*
 *STEP 1c - Compress fastq files for storage
 */

process gzip_fastq {
    tag "$name"
    publishDir "${params.outdir}/fastq", mode: 'copy'

    when:
    params.savefq || params.saveAllfq

    input:
    set val(name), file(fastq_reads) from fastq_reads_gzip.mix(fastq_reads_gzip_sra)

    output:
    set val(name), file("*.gz") into compressed_fastq

    script:
    """
    gzip -c ${name}.fastq > ${name}.fastq.gz
    """
 }


/*
 * STEP 2a - Trimming
 */

process bbduk {
    validExitStatus 0,1
    tag "$name"
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file ("*.trim.fastq") into trimmed_reads_fastqc, trimmed_reads_hisat2, trimmed_reads_gzip
    file "*.txt" into trim_stats

    script:
//    prefix = fastq.baseName
    bbduk_mem = task.memory.toGiga()
    if (!params.singleEnd && params.flip) {
        """
        echo ${name}

        seqkit seq -j 16 -r -p \
                  ${name}_R1.flip.fastq \
                  -o ${name}.flip.fastq
                  
        seqkit seq -j 16 -r -p \
                 ${name}_R2.flip.fastq \
                 -o ${name}.flip.fastq

        

        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${name}_R1.flip.fastq \
                  in2=${name}_R2.flip.fastq \
                  out=${name}_R1.flip.trim.fastq \
                  out2=${name}_R2.flip.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt \
                  ehist=${name}.ehist.txt
        """
    } else if (params.flip) {
        """
        echo ${name}


        seqkit seq -j 16 -r -p \
                  ${name}.fastq \
                  -o ${name}.flip.fastq

        
        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${name}.flip.fastq \
                  out=${name}.flip.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt \
                  ehist=${name}.ehist.txt
        """
    }
        else if (!params.singleEnd) {
        """
        echo ${name}      

        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${name}_R1.fastq \
                  in2=${name}_R2.fastq \
                  out=${name}_R1.trim.fastq \
                  out2=${name}_R2.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt \
                  ehist=${name}.ehist.txt
        """
    } else {
        """
        echo ${name}
        echo ${bbduk_mem}
        
        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${name}.fastq \
                  out=${name}.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt \
                  ehist=${name}.ehist.txt
        """
    }
}


/*
 * STEP 2b - Trimmed FastQC
 */

process fastqc_trimmed {
    validExitStatus 0,1
    tag "$prefix"
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    prefix = trimmed_reads.baseName
    """
    echo ${prefix}

    fastqc $trimmed_reads
		extract_fastqc_stats.sh --srr=${prefix} > ${prefix}_stats_fastqc.txt
    """
}

/*
 *STEP 2c - Compress trimmed fastq files for storage
 */

process gzip_trimmed {
    tag "$prefix"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    when:
    params.saveTrim || params.saveAllfq

    input:
    file(trimmed_reads) from trimmed_reads_gzip

    output:
    set val(prefix), file("*.gz") into trimmed_gzip

    script:
    prefix = trimmed_reads.baseName
    """
    gzip -c $trimmed_reads > ${prefix}.fastq.gz
    """
 }


/*
 * STEP 3 - Map reads to reference genome
 */

process hisat2 {
    // NOTE: this tool sends output there even in successful (exit code 0) 
    // termination, so we have to ignore errors for now, and the next 
    // process will blow up from missing a SAM file instead.
    //errorStrategy 'ignore'
    tag "$name"
    validExitStatus 0,143

    input:
    val(indices) from hisat2_indices.first()
    set val(name), file(trimmed_reads) from trimmed_reads_hisat2

    output:
    set val(name), file("*.sam") into hisat2_sam

    script:
    index_base = indices[0].toString() - ~/.\d.ht2/
    if (!params.singleEnd) {
        """
        echo ${name}
    
        hisat2  -p ${task.cpus} \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${index_base} \
                -1 ${name}_R1.trim.fastq \
                -2 ${name}_R2.trim.fastq
                > ${name}.sam
        """
    } else {
        """
        echo ${name}
    
        hisat2  -p ${task.cpus} \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${index_base}\
                -U ${trimmed_reads} \
                > ${name}.sam
        """
    }
}


/*
 * STEP 4 - Convert to BAM format and sort
 */

/*
 * STEP 4 - Convert to BAM format and sort
 */

process samtools {
    tag "$prefix"
    publishDir "${params.outdir}/mapped/bams", mode: 'copy', pattern: "${prefix}.sorted.bam"
    publishDir "${params.outdir}/mapped/bams", mode: 'copy', pattern: "${prefix}.sorted.bam.bai"
    publishDir "${params.outdir}/qc/mapstats", mode: 'copy', pattern: "${prefix}.sorted.bam.flagstat"
    publishDir "${params.outdir}/qc/mapstats", mode: 'copy', pattern: "${prefix}.sorted.bam.millionsmapped"

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${prefix}.sorted.bam") into sorted_bam_ch
    set val(name), file("${prefix}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${prefix}.sorted.bam.flagstat") into bam_flagstat
    set val(name), file("${prefix}.sorted.bam.millionsmapped") into bam_milmapped_bedgraph

    script:
    prefix = mapped_sam.baseName
// Note that the millionsmapped arugments below are only good for SE data. When PE is added, it will need to be changed to:
    // -F 0x40 rootname.sorted.bam | cut -f1 | sort | uniq | wc -l  > rootname.bam.millionsmapped
    if (!params.singleEnd) {
    """

    samtools view -@ ${task.cpus} -bS -o ${prefix}.bam ${mapped_sam}
    samtools sort -@ ${task.cpus} ${prefix}.bam > ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools view -@ ${task.cpus} -F 0x40 ${prefix}.sorted.bam | cut -f1 | sort | uniq | wc -l > ${prefix}.sorted.bam.millionsmapped
    samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
    """
    } else {
    """

    samtools view -@ ${task.cpus} -bS -o ${prefix}.bam ${mapped_sam}
    samtools sort -@ ${task.cpus} ${prefix}.bam > ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools view -@ ${task.cpus} -F 0x904 -c ${prefix}.sorted.bam > ${prefix}.sorted.bam.millionsmapped
    samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
    """
    }
}

sorted_bam_ch
   .into {sorted_bams_for_bedtools_bedgraph; sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup}

sorted_bam_indices_ch
    .into {sorted_bam_indices_for_bedtools_bedgraph; sorted_bam_indices_for_bedtools_normalized_bedgraph; sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc}

/*
 *STEP 5a - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """

    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}

    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 5b - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    validExitStatus 0,143
    publishDir "${params.outdir}/qc/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else filename
        }

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc
    file(bam_indices) from sorted_bam_indices_for_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """

    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_dist.txt

    read_duplication.py -i ${bam_file} \
                        -o ${name}.read_duplication

    infer_experiment.py -i ${bam_file} \
                        -r ${genome_refseq} \
                        > ${name}.infer_experiment.txt
    """
 }



/*
 *STEP 5c - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    pileup_mem = task.memory.toGiga()
    """

    pileup.sh -Xmx${pileup_mem}g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }

/*
 *STEP 6a - Create non-normalzied bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "*{neg,pos}.bedGraph"
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"
    publishDir "${params.outdir}/mapped/rcc_bedgraphs", mode: 'copy', pattern: "${name}.rcc.bedGraph"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
    set val(name), file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("*.bedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    """

    genomeCoverageBed \
                     -bg \
                     -strand + \
                     -g hg38 \
                     -ibam ${bam_file} \
                     > ${name}.pos.bedGraph

    genomeCoverageBed \
                     -bg \
                     -strand - \
                     -g hg38 \
                     -ibam ${bam_file} \
                     > ${name}.tmp.neg.bedGraph

     awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}.tmp.neg.bedGraph \
        > ${name}.neg.bedGraph
        rm ${name}.tmp.neg.bedGraph

    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph

    sortBed \
             -i ${name}.unsorted.bedGraph \
             > ${name}.bedGraph

    rm ${name}.unsorted.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph \

    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph

    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph
    rm ${name}.unsorted.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph

    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    rm ${name}.unsorted.neg.rcc.bedGraph

    """
 }


/* Idea borrowed from the nf-core/atacseq workflow:
 * Just generate the chromosome sizes file using samtools, if not provided.
  */
if(!params.chrom_sizes) {
  process make_chromosome_sizes {
      tag "$fasta"
      publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'

      input:
      file fasta from genome_fasta

      output:
      file("${fasta}.sizes") into chrom_sizes_ch

      script:
      """
      samtools faidx $fasta
      cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
      """
  }
}

chrom_sizes_ch.into{chrom_sizes_for_bed; chrom_sizes_for_bigwig; chrom_sizes_for_igv}


/*
 *STEP 6b - Create bedGraphs and bigwigs for dREG
 */

process dreg_prep {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    publishDir "${params.outdir}/mapped/dreg_input", mode: 'copy', pattern: "*.bw"

    input:
    set val(name), file(bam_file) from sorted_bams_for_dreg_prep
    file(chr_sizes) from chrom_sizes_for_bed

    output:
        set val(name), file("*.bw") into dreg_bigwig

    script:
    """

    echo "Creating BigWigs suitable as inputs to dREG"

    bedtools bamtobed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sortBed -i ${name}.dreg.bed > ${name}.dreg.sort.bed

    echo positive strand processed to bedGraph

    bedtools genomecov -bg -i ${name}.dreg.sort.bed -g ${chr_sizes} -strand + > ${name}.pos.bedGraph
    sortBed -i ${name}.pos.bedGraph > ${name}.pos.sort.bedGraph
    bedtools genomecov -bg -i ${name}.dreg.sort.bed -g ${chr_sizes} -strand - \
    | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' > ${name}.neg.bedGraph
    sortBed -i ${name}.neg.bedGraph > ${name}.neg.sort.bedGraph

    echo negative strand processed to bedGraph

    ${params.bedGraphToBigWig} ${name}.pos.sort.bedGraph ${chr_sizes} ${name}.pos.bw
    ${params.bedGraphToBigWig} ${name}.neg.sort.bedGraph ${chr_sizes} ${name}.neg.bw

    echo bedGraph to bigwig done
    """
 }

/*
 *STEP 7 - Normalize bigWigs by millions of reads mapped for visualization on nascent2.0
 */

process normalized_bigwigs {
    validExitStatus 0
    tag "$name"
    publishDir "${params.outdir}/mapped/rcc_bigwig", mode: 'copy'

    input:
    set val(name), file(neg_bedgraph) from bedgraph_bigwig_neg
    set val(name), file(pos_bedgraph) from bedgraph_bigwig_pos
    file chrom_sizes from chrom_sizes_for_bigwig

    output:
    set val(name), file("*.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.bedGraphToBigWig} ${pos_bedgraph} ${chrom_sizes} ${name}.pos.rcc.bw
    ${params.bedGraphToBigWig} ${neg_bedgraph} ${chrom_sizes} ${name}.neg.rcc.bw

    """
}

/*
 *STEP 8 - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    tag "$name"
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file chrom_sizes from chrom_sizes_for_igv

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    export LC_ALL=C
    igvtools toTDF ${normalized_bg} ${name}.rcc.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 9 - MultiQC
 */
process multiqc {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

//TO DO : Need to build a new multiqc container for the newest version

    """
    export PATH=~/.local/bin:$PATH

    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/nascent] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/nascent] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/nascent] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/nascent] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/nascent] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/nascent] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/nascent]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/nascent]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/nascent v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
