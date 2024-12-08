process CONCAT_TRANSCRIPT_BEDS {
    tag "${meta.id}"

    input:
    // grohmm, homer, pints
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("*combined_transcripts.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    """
    # Add prefix to GRO-HMM transcripts
    if [ -f "${beds[0]}" ]; then
        awk 'BEGIN{OFS="\t"} {\$4="grohmm_"\$4; print}' ${beds[0]} > grohmm_prefixed.bed
    else
        touch grohmm_prefixed.bed
    fi

    # Add prefix to HOMER transcripts
    if [ -f "${beds[1]}" ]; then
    awk 'BEGIN{OFS="\t"} {\$4="homer_"\$4; print}' ${beds[1]} > homer_prefixed.bed
    else
        touch homer_prefixed.bed
    fi

    # Add prefix to PINTS transcripts
    if [ -f "${beds[2]}" ]; then
    awk 'BEGIN{OFS="\t"} {\$4="pints_"\$4; print}' ${beds[2]} > pints_prefixed.bed
    else
        touch pints_prefixed.bed
    fi

    # Concatenate all files
    cat grohmm_prefixed.bed homer_prefixed.bed pints_prefixed.bed > ${meta.id}_combined_transcripts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
