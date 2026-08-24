"""Top-level ``snakemake`` file that runs analysis.

The core analysis is done by the `seqneut-pipeline` submodule; the additional analyses
specific to this project are in the `.smk` files in `rules/`, each of which defines the
final outputs it contributes to `rule all`.

"""


configfile: "config.yml"


include: "seqneut-pipeline/seqneut-pipeline.smk"
include: "rules/library_qc.smk"
include: "rules/analyze_titers.smk"
include: "rules/trees.smk"


rule all:
    input:
        seqneut_pipeline_outputs,
        library_qc_outputs,
        analyze_titers_outputs,
        trees_outputs,
