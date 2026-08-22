"""Validation of the viral libraries."""


rule validate_viral_library:
    """Validate a viral library."""
    input:
        csv=lambda wc: config["viral_libraries"][wc.viral_library],
    output:
        validation="results/validate_viral_library/{viral_library}_validation.txt",
    log:
        "results/logs/validate_viral_library_{viral_library}.txt",
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        validations=config["viral_library_validations"],
    script:
        "../scripts/validate_viral_library.py"


validate_viral_library_outputs = expand(
    rules.validate_viral_library.output.validation,
    viral_library=config["viral_libraries"],
)
