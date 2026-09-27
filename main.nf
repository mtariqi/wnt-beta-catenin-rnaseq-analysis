nextflow.enable.dsl = 2

include { WNTPCP_FUNMAP } from './workflow/modules/wntpcp_funmap'


workflow {
    WNTPCP_FUNMAP(
        file(params.pipeline_script, checkIfExists: true),
        file(params.deseq, checkIfExists: true),
        file(params.orthology, checkIfExists: true),
        file(params.funmap, checkIfExists: true),
        file(params.string_file, checkIfExists: true),
        file(params.literature, checkIfExists: true),
        file(params.gene_file, checkIfExists: true),
        params.string_species,
        params.prefix,
        params.title,
        params.funmap_score_min,
        params.string_score_min,
        params.max_edges,
        params.seed
    )
}
