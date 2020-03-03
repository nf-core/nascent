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
        SRA tools -- fasterq-dump sra to generate fastq file
        FastQC (pre-trim) -- perform pre-trim FastQC on fastq files
    2. Trimming & Mapping
        BBDuk -- trim fastq files for quality and adapters
        FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
        HISAT2 -- Map reads to a reference genome
    3. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM, BAM --> CRAM, index CRAM
    4. Quality control
        preseq -- estimate library complexity
        RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio
        Picard -- Mark duplicates, GC content
        NQC -- Perform nascent-specific quality control analysis
    5. Mapping Visualization
        BEDTools : non-normalized & nornmalized bedgraphs
        BEDTools and kentUtils : 5' bigwigs for dREG & normalized bigwigs for genome browser     
        IGV Tools : bedGraph --> tdf
    6. MultiQC : generate QC report for pipeline
    7. Read Counting : BEDTools multicov to count reads over genes
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/nascent -profile singularity --reads '/project/*_{R1,R2}*.fastq' --outdir '/project/'
    or:
    nextflow run nf-core/nascent --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Required arguments:
         -profile                      Configuration profile to use. <base, singularity>
         --genomeid                    Genome ID (e.g. hg38, mm10, rn6, etc.).
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --single_end                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.        
        
    Strandedness:
        --forwardStranded              The library is forward stranded
        --reverseStranded              The library is reverse stranded
        --unStranded                   The default behaviour        
        
    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Saves compressed raw fastq reads.
        --saveTrim                     Saves compressed trimmed fastq reads.
        --saveBAM                      Save BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Saves all compressed fastq reads.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --saveReference               Save the generated reference files the the Results directory.

    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipFastQC                   Skip running FastQC.
        --skipRSeQC                    Skip running RSeQC.
        --skippicard                   Skip running picard.        
        --skippreseq                   Skip running preseq.
        --skippileup                   Skip running pileup.sh.
        --skipAllQC                    Skip running all QC.
        
    Analysis Options:
        --noTrim                       Skip trimming and map only. Will also skip flip/flipR2 (any BBMap) steps.
        --counts                       Run BEDTools mutlicov for each sample to obtain gene counts over the RefSeq annotation.
        --dreg                         Produce bigwigs formatted for input to dREG.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.name = false
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/assets/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.merge_counts = "$baseDir/bin/merge_counts.py"
params.workdir = "./nextflowTemp"
output_docs = file("$baseDir/docs/output.md")
params.extract_fastqc_stats = "$baseDir/bin/extract_fastqc_stats.sh"

// Validate inputs

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

if ( params.fasta ){
   // genome_fasta = file(params.fasta)
   // if( !genome_fasta.exists() ) exit 1, "Genome directory not found: ${params.fasta}"
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
           .into { genome_fasta; ch_fasta_for_hisat_index; ch_fasta_for_samtools; ch_fasta_for_picard}
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

chrom_sizes_ch.into{chrom_sizes_for_bed; chrom_sizes_for_dreg; chrom_sizes_for_bigwig; chrom_sizes_for_igv }


if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

//if ( params.picard_path ){
//    picard_path = file(params.picard_path)
//}

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

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (params.reads) {
    if (params.single_end) {
        fastq_reads_qc = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_hisat2_notrim = Channel
                            .fromPath(params.reads)
                            .map { file -> tuple(file.simpleName, file) }        
    } else {
        Channel
            .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
if(params.reads) summary['Reads']            = params.reads
if(params.sras) summary['SRAs']             = params.sras
summary['Genome Ref']       = params.genome
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Strandedness']     = (params.unStranded ? 'None' : params.forwardStranded ? 'Forward' : params.reverseStranded ? 'Reverse' : 'None')
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save BAM']         = params.saveBAM ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run Multicov']     = params.counts ? 'YES' : 'NO'
summary['Skip Trimming']    = params.noTrim ? 'NO' : 'YES'
summary['Run FastQC']       = params.skipFastQC ? 'NO' : 'YES'
summary['Run preseq']       = params.skippreseq ? 'NO' : 'YES'
summary['Run pileup']       = params.skippileup ? 'NO' : 'YES'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Skip All QC']      = params.skipAllQC ? 'YES' : 'NO'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['dREG']             = params.dreg ? 'YES' : 'NO'
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile']   = workflow.profile

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
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """

    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    preseq 2> v_preseq.txt
    bedtools --version > v_bedtools.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    picard CollectGcBiasMetrics --version &> v_collectgcbiasmetrics.txt  || true
    export LC_ALL=C
    igvtools version > v_igv-tools.txt
    infer_experiment.py --version > v_rseqc.txt
    multiqc --version > v_multiqc.txt
    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1 -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus ${task.cpus} }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_hisat2_notrim_sra
   

    script:
    prefix = reads.baseName
    if (!params.threadfqdump) {
        """
        echo ${prefix}
        fastq-dump ${reads} --gzip
        """
    } else if (!params.single_end) {
         """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads ${task.cpus} \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.single_end) {
        """
        echo ${prefix}
        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads ${task.cpus} \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * PREPROCESSING - Build HISAT2 index (borrowed from nf-core/rnaseq)
 */
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

process fastQC {
    validExitStatus 0,1
    tag "$prefix"
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)     "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)    "qc/fastqc/$filename"
        else if (filename.indexOf("txt") > 0)     "qc/fastqc_stats/$filename"
        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html}" into fastqc_results
    file "*.fastqc_stats.txt" into fastqc_stats

    script:
    """
    echo ${prefix}
    fastqc $reads
    
    ${params.extract_fastqc_stats} \
        --srr=${prefix} \
        > ${prefix}.fastqc_stats.txt    
    """
}



/*
 * STEP 2 - Trimming & Mapping
 */

process bbduk_hisat2 {
    validExitStatus 0
    tag "$name"
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.{refstats,trimstats}.txt"
    publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*hisat2_mapstats.txt"    
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }
    
    when:
    !params.noTrim

    input:
    val(indices) from hisat2_indices.first()
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file("*.trim.fastq.gz") into trimmed_reads_fastqc
    file "*.{refstats,trimstats}.txt" into trim_stats
    set val(name), file("*.sam") into hisat2_sam
    file("*hisat2_mapstats.txt") into hisat2_mapstats       

    script:
    index_base = indices[0].toString() - ~/.\d.ht2/
    bbduk_mem = task.memory.toGiga()
    prefix_pe = reads[0].toString() - ~/(_1\.)?(_R1)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
    prefix_se = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    
    def rnastrandness = ''
    if (params.forwardStranded && !params.unStranded){
        rnastrandness = params.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (params.reverseStranded && !params.unStranded){
        rnastrandness = params.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }
    
    if (!params.single_end && params.flip) {
        """
        echo ${name}         
        reformat.sh -Xmx${bbduk_mem}g \
                t=${task.cpus} \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}_1.flip.fastq.gz \
                out2=${prefix_pe}_2.flip.fastq.gz \
                rcomp=t
                
        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${prefix_pe}_1.flip.fastq.gz \
                  in2=${prefix_pe}_2.flip.fastq.gz \
                  out=${prefix_pe}_1.flip.trim.fastq.gz \
                  out2=${prefix_pe}_2.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_pe}.trimstats.txt \
                  refstats=${prefix_pe}.refstats.txt
                  
        hisat2 -p ${task.cpus} \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${index_base} \
               -1 ${prefix_pe}_1.flip.trim.fastq.gz \
               -2 ${prefix_pe}_2.flip.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                  
        """
    } else if (params.single_end && params.flip) {
        """
        echo ${name}        
        reformat.sh -Xmx${bbduk_mem}g \
                t=${task.cpus} \
                in=${reads} \
                out=${prefix_se}.flip.fastq.gz \
                rcomp=t
        
        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${prefix_se}.flip.fastq.gz \
                  out=${prefix_se}.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_se}.trimstats.txt \
                  refstats=${prefix_se}.refstats.txt
                  
        hisat2  -p ${task.cpus} \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${index_base}\
                -U ${prefix_se}.flip.trim.fastq.gz  \
                $rnastrandness \
                --new-summary \
                > ${prefix_se}.sam \
                2> ${prefix_se}.hisat2_mapstats.txt                  
        """
    } else if (!params.single_end && params.flipR2) {
                """
        echo ${prefix_pe}

        reformat.sh -Xmx${bbduk_mem}g \
                t=${task.cpus} \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}.flip.fastq.gz \
                out2=${prefix_pe}.flip.fastq.gz \
                rcompmate=t

        bbduk.sh -Xmx${bbduk_mem}g \
                t=${task.cpus} \
                in=${prefix_pe}.flip.fastq.gz \
                in2=${prefix_pe}.flip.fastq.gz \
                out=${prefix_pe}.flip.trim.fastq.gz \
                out2=${prefix_pe}.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_pe}.trimstats.txt \
                refstats=${prefix_pe}.refstats.txt
                
        hisat2 -p ${task.cpus} \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${index_base} \
               -1 ${prefix_pe}_1.flip.trim.fastq.gz \
               -2 ${prefix_pe}_2.flip.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                   
        """
    } else if (!params.single_end) {
        """
        echo ${prefix_pe}

        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${reads[0]} \
                  in2=${reads[1]} \
                  out=${prefix_pe}_1.trim.fastq.gz \
                  out2=${prefix_pe}_2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_pe}.trimstats.txt \
                  refstats=${prefix_pe}.refstats.txt
                  
        hisat2 -p ${task.cpus} \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${index_base} \
               -1 ${prefix_pe}_1.trim.fastq.gz \
               -2 ${prefix_pe}_2.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                          
        """
    } else {
        """
        echo ${prefix_se}

        bbduk.sh -Xmx${bbduk_mem}g \
                  t=${task.cpus} \
                  in=${reads} \
                  out=${prefix_se}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_se}.trimstats.txt \
                  refstats=${prefix_se}.refstats.txt
                  
        hisat2  -p ${task.cpus} \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${index_base}\
                -U ${prefix_se}.trim.fastq.gz \
                $rnastrandness \
                --new-summary \
                > ${prefix_se}.sam \
                2> ${prefix_se}.hisat2_mapstats.txt                  
        """
    } 
}


/*
 * STEP 2+ - Trimmed FastQC
 */

process fastQC_trim {
    validExitStatus 0,1
    tag "$name"
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)   "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)  "qc/fastqc/$filename"
        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC && !noTrim

    input:
    set val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    script:
    """
    echo ${name}
    
    fastqc --quiet ${trimmed_reads}
    """
}


/*
 * STEP 2+ - Map reads to reference genome w/o trimming
 */

if (params.noTrim) {
    process hisat2 {
        tag "$name"
        validExitStatus 0
        publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*.txt"
    
        input:
        val(indices) from hisat2_indices.first()
        set val(name), file(reads) from fastq_reads_hisat2_notrim_sra.mix(fastq_reads_hisat2_notrim)
    
        output:
        set val(name), file("*.sam") into hisat2_sam
        file("*.txt") into hisat2_mapstats
    
        script:
        index_base = indices[0].toString() - ~/.\d.ht2/
        prefix_pe = trimmed_reads[0].toString() - ~/(_1\.)?(_R1)?(flip)?(trim)?(\.flip)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
        prefix_se = trimmed_reads[0].toString() - ~/(\.flip)?(\.fq)?(\.fastq)?(\.gz)?$/
        
        def rnastrandness = ''
        if (params.forwardStranded && !params.unStranded){
            rnastrandness = params.single_end ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (params.reverseStranded && !params.unStranded){
            rnastrandness = params.single_end ? '--rna-strandness R' : '--rna-strandness RF'
        }
        
        if (!params.single_end) {
            """
            echo ${prefix_pe}
        
            hisat2 -p ${task.cpus} \
                   --very-sensitive \
                   --no-spliced-alignment \
                    $rnastrandness \
                   -x ${index_base} \
                   -1 ${reads[0]} \
                   -2 ${reads[1]} \
                   --new-summary \
                   > ${prefix_pe}.sam \
                   2> ${prefix_pe}.hisat2_mapstats.txt                
            """
        }
        else {
            """
            echo ${prefix_se}
        
            hisat2  -p ${task.cpus} \
                    --very-sensitive \
                    --no-spliced-alignment \
                    $rnastrandness \
                    -x ${index_base} \
                    -U ${reads} \
                    --new-summary \
                    > ${prefix_se}.sam \
                    2> ${prefix_se}.hisat2_mapstats.txt                
            """
        }
    }
}

/*
 * STEP 3 - Convert to BAM/CRAM format and sort
 */

process samtools {
    tag "$name"
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if ((filename.indexOf("sorted.bam") > 0) & params.saveBAM)     "mapped/bams/$filename"
        else if ((filename.indexOf("sorted.bam.bai") > 0) & params.saveBAM) "mapped/bams/$filename"
        else if (filename.indexOf("flagstat") > 0)                          "qc/mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)                    "qc/mapstats/$filename"
        else if (filename.indexOf("sorted.cram") > 0)                       "mapped/crams/$filename"
        else if (filename.indexOf("sorted.cram.crai") > 0)                  "mapped/crams/$filename"
        else null            
    }

    input:
    set val(name), file(mapped_sam) from hisat2_sam
    file fasta from ch_fasta_for_samtools

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${name}.flagstat") into bam_flagstat
    set val(name), file("${name}.millionsmapped") into bam_milmapped_bedgraph
    set val(name), file("${name}.sorted.cram") into cram_ch
    set val(name), file("${name}.sorted.cram.crai") into cram_index_ch

    script:
    if (!params.single_end) {
    """
    samtools view -@ ${task.cpus} -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ ${task.cpus} ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ ${task.cpus} -F 0x40 ${name}.sorted.bam | cut -f1 | sort | uniq | wc -l > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ ${task.cpus} -C -T ${fasta} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ ${task.cpus} -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    } else {
    """
    samtools view -@ ${task.cpus} -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ ${task.cpus} ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ ${task.cpus} -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ ${task.cpus} -C -T ${fasta} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ ${task.cpus} -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    }
}

sorted_bam_ch
   .into { sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup; sorted_bams_for_picard; sorted_bam_for_bedgraph }

sorted_bam_indices_ch
    .into { sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc; sorted_bam_indices_for_picard; sorted_bam_index_for_bedgraph }

cram_ch
    .into { cram_for_counts; cram_dreg_prep }

cram_index_ch
    .into { cram_index_for_counts; cram_index_dreg_prep} 

/*
 *STEP 4+ - Picard tools
 */

process picard {
    tag "$name"
    errorStrategy 'ignore'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("marked_dup_metrics.txt") > 0)            "qc/picard/dups/$filename"
        else if (filename.indexOf("gc_bias_metrics.pdf") > 0)               "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("gc_bias_metrics.txt") > 0)               "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("summary_metrics.txt") > 0)               "qc/picard/gc_bias/$filename"
        else null            
    }
    
    when:
    !params.skippicard && !params.skipAllQC    

    input:
    set val(name), file(bam_file) from sorted_bams_for_picard
    file(bam_indices) from sorted_bam_indices_for_picard
    file fasta from ch_fasta_for_picard

    output:
    file "*.{txt,pdf}" into picard_stats_multiqc
    
    script:
    markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
    """
    picard ${markdup_java_options} MarkDuplicates \
         I=${bam_file} \
         O=${name}.marked_duplicates.bam \
         M=${name}.marked_dup_metrics.txt             
         
    picard ${markdup_java_options} CollectGcBiasMetrics \
          I=${bam_file} \
          O=${name}.gc_bias_metrics.txt \
          CHART=${name}.gc_bias_metrics.pdf \
          S=${name}.summary_metrics.txt \
          R=${fasta}    
    """
}

/*
 *STEP 4+ - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"
    
    when:
    !params.skippreseq && !params.skipAllQC    

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
 *STEP 4+ - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    validExitStatus 0,1,143
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
            else null
        }
    
    when:
    !params.skipRSeQC && !params.skipAllQC

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
 *STEP 4+ - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"
    
    when:
    !params.skippileup && !params.skipAllQC    

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
 *STEP 5 - Create non-normalzied bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"      

    input:
    set val(name), file(bam_file) from sorted_bam_for_bedgraph
    set val(name), file(bam_indices) from sorted_bam_index_for_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph
    file chrom_sizes from chrom_sizes_for_bed

    output:
    set val(name), file("${name}.pos.bedGraph") into pos_non_normalized_bedgraphs
    set val(name), file("${name}.neg.bedGraph") into neg_non_normalized_bedgraphs
    set val(name), file("${name}.bedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    if (params.single_end) {
    """    
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        > ${name}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.neg.bedGraph

    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    """
    } else {
    """   
    samtools view \
        -h -b -f 0x0040 \
        ${bam_file} \
        > ${name}.first_pair.bam
        
    samtools view \
        -h -b -f 0x0080 \
        ${bam_file} \
        > ${name}.second_pair.bam
        
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        > ${name}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.first_pair.neg.bedGraph
                     
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | sortBed \
        > ${name}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${name}.second_pair.neg.bedGraph
                     
    unionBedGraphs \
        -i ${name}.first_pair.pos.bedGraph ${name}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.pos.bedGraph 
        
    unionBedGraphs \
        -i ${name}.first_pair.neg.bedGraph ${name}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.neg.bedGraph 
    
    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph     
    """
    }
 }


/*
 *STEP 5+ - Create bedGraphs and bigwigs for dREG
 */

process dreg_prep {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    publishDir "${params.outdir}/mapped/dreg_input", mode: 'copy', pattern: "*.bw"
    
    when:
    params.dreg

    input:
    set val(name), file(cram_file) from cram_dreg_prep
    set val(name), file(cram_index) from cram_index_dreg_prep
    file(chrom_sizes) from chrom_sizes_for_dreg

    output:
    set val(name), file("*.bw") into dreg_bigwig

    script:
    """
    echo "Creating BigWigs suitable as inputs to dREG"
    
    export CRAM_REFERENCE=${genome}    
    
    bamToBed -i ${cram_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sortBed -i ${name}.dreg.bed > ${name}.dreg.sort.bed
    
    echo positive strand processed to bedGraph
    
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand + \
            > ${name}.pos.bedGraph
    sortBed \
            -i ${name}.pos.bedGraph \
            > ${name}.pos.sort.bedGraph
            
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${name}.neg.bedGraph
    sortBed \
            -i ${name}.neg.bedGraph \
            > ${name}.neg.sort.bedGraph
            
    echo negative strand processed to bedGraph
    
    ${params.bedGraphToBigWig} ${name}.pos.sort.bedGraph ${chrom_sizes} ${name}.pos.bw
    ${params.bedGraphToBigWig} ${name}.neg.sort.bedGraph ${chrom_sizes} ${name}.neg.bw
    
    echo bedGraph to bigwig done
    """
 }

/*
 *STEP 5+ - Normalize bigWigs by millions of reads mapped for visualization
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
 *STEP 5+ - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
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
 * STEP 6 - MultiQC
 */

process multiQC {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC && !params.skipAllQC

    input:
    file multiqc_config from ch_multiqc_config.collect()
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from ch_software_versions_yaml
    file ('qc/hisat2_mapstats/*') from hisat2_mapstats.collect()
    file ('qc/picard/*') from picard_stats_multiqc.collect()    

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}

/*
 * STEP 10 - Counts -- BEDTools multicov
 */

process multicov {
    tag "$name"
    validExitStatus 0
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "*.bed"
    
    when:
    params.counts
    
    input:
    set val(name), file (cram_count) from cram_for_counts
    set val(name), file (cram_index) from cram_index_for_counts
       
    output:
    file ("*.bed") into counts_bed_out
        
    script:
        """
        export CRAM_REFERENCE=${genome}
        
        bedtools multicov \\
            -bams ${cram_count} \\
            -s \\
            -bed ${genome_refseq} \\
            > ${name}_counts.bed
        """
}

/*
 * STEP 10+ - Merge Counts
 */

process merge_multicov {
    validExitStatus 0
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "merged_gene_counts.bed"
    
    when:
    params.counts
    
    input:
    file (multicov:'counts/*') from counts_bed_out.collect()
       
    output:
    file ("merged_gene_counts.bed") into merged_counts_bed_out
        
    script:
        """
        python3 ${params.merge_counts} \\
            -b './counts/' \\
            -o merged_gene_counts.bed
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

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/nascent] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/nascent] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/nascent]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/nascent]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/nascent v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
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
