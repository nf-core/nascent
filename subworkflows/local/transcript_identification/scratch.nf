workflow {
    def chr = Channel.from( 1..4 )
    def samples = Channel.from( "jurkat", "hepg2" )

    RUN_FOR_EACH_CHR (
        samples,
        chr
    )

    RUN_FOR_EACH_CHR.out.beds
    .groupTuple(by: [0])
    .map { meta, beds -> [meta, beds.flatten()] }
    .view()

}

process RUN_FOR_EACH_CHR {
    input:
    val(meta)
    each chr_name // optional

    output:
    tuple val(meta), path("*.bed")     , optional:true, emit: beds

    script:
    def prefix = "${meta}" + (chr_name ? '_' + chr_name : '_all')
    """
    echo "Processing ${prefix} 1" > ${prefix}_1.bed
    echo "Processing ${prefix} 2" > ${prefix}_2.bed
    echo "Processing ${prefix} 3" > ${prefix}_3.bed
    """
}
