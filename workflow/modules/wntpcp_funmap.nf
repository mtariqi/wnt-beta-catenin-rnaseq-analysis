process WNTPCP_FUNMAP {
    tag "${prefix}"

    container params.container_image

    publishDir params.outdir,
        mode: 'copy',
        overwrite: true

    input:
    path pipeline_script
    path deseq
    path orthology
    path funmap
    path string_file
    path literature
    path gene_file
    val string_species
    val prefix
    val title
    val funmap_score_min
    val string_score_min
    val max_edges
    val seed

    output:
    path "${prefix}.png", emit: png
    path "${prefix}.pdf", emit: pdf
    path "${prefix}.svg", emit: svg
    path "${prefix}_edge_audit.csv", emit: edge_audit
    path "${prefix}_node_audit.csv", emit: node_audit
    path "${prefix}_run_provenance.json", emit: provenance

    script:
    """
    python ${pipeline_script} \
        --deseq ${deseq} \
        --orthology ${orthology} \
        --funmap ${funmap} \
        --string ${string_file} \
        --string-species ${string_species} \
        --literature ${literature} \
        --gene-file ${gene_file} \
        --strict-orthology \
        --funmap-score-min ${funmap_score_min} \
        --string-score-min ${string_score_min} \
        --max-edges ${max_edges} \
        --seed ${seed} \
        --outdir . \
        --prefix ${prefix} \
        --title "${title}"
    """
}
