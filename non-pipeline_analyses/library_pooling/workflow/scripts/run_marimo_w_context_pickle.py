"""Run a marimo notebook non-interactively with snakemake context via a pickle.

Adapted from ``seqneut-pipeline/scripts/run_marimo_w_context_pickle.py``. Invoked
as a Snakemake ``script:``; it pickles the relevant fields of the ``snakemake``
object, exports the notebook to a self-contained HTML with that pickle passed on
the command line, and post-processes the HTML to hide the app code.

The paired notebook cell detects the ``--context-pickle`` argument to decide
whether it is running under Snakemake (this driver) or interactively via
``marimo edit`` -- see ``notebooks/example_analysis.py``.
"""

# `snakemake` is a global injected by Snakemake when this file is run as a
# `script:`; ruff can't see it, so silence undefined-name for the whole file.
# ruff: noqa: F821

import os
import pickle
import subprocess
import sys

# Redirect stderr and stdout to the Snakemake log.
log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

# Build a plain-dict context (picklable) with a flat input/output/params layout
# matching the seqneut-pipeline notebooks. Snakemake disallows function-based
# rule outputs, so the notebook's own output files can't live in
# `snakemake.output`; they are passed via `params.notebook_output` and merged
# into context["output"] here. Likewise config params come via
# `params.notebook_params`.
#
# context["input"] excludes the notebook file itself (`marimo_nb`) so it holds
# only the analysis inputs the notebook reads by key.
nb_input = {k: v for k, v in dict(snakemake.input).items() if k != "marimo_nb"}
nb_output = dict(snakemake.params.notebook_output)
nb_params = dict(snakemake.params.notebook_params)

context = {
    "workdir": os.getcwd(),
    "input": nb_input,
    "output": nb_output,
    "params": nb_params,
    "threads": snakemake.threads,
    "resources": dict(snakemake.resources),
    "wildcards": dict(snakemake.wildcards),
}

marimo_nb = snakemake.input.marimo_nb
marimo_html = snakemake.output.marimo_html
context_pickle = snakemake.output.context_pickle

print(f"Running marimo notebook: {marimo_nb=}")
print(f"Using context pickle: {context_pickle=}")
print(f"Context keys: {list(context.keys())=}")
for key in context:
    if context[key]:
        print(f"  For {key}, context is: {context[key]}")
print(f"Output HTML: {marimo_html=}")

# Write context to the pickle the notebook will read.
with open(context_pickle, "wb") as f:
    pickle.dump(context, f, protocol=pickle.HIGHEST_PROTOCOL)

# Export the notebook to a self-contained HTML, passing the pickle to the
# notebook via marimo's post-`--` script args.
cmd = [
    "marimo",
    "export",
    "html",
    marimo_nb,
    "-o",
    marimo_html,
    "--",
    "--context-pickle",
    context_pickle,
]

print(f"Running command: {' '.join(cmd)}")
print("=" * 80)

try:
    result = subprocess.run(
        cmd,
        stdout=log_file,
        stderr=subprocess.STDOUT,
        check=False,
    )
    print("=" * 80)
    if result.returncode != 0:
        print(f"ERROR: marimo export failed with return code {result.returncode}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    print("marimo export completed successfully")
except Exception as e:
    print(f"ERROR: Exception running marimo export: {e}")
    raise

# Post-process HTML to hide the app code in the rendered output.
print("Post-processing HTML to set showAppCode to false...")
with open(marimo_html, "r") as f:
    html_content = f.read()

old_string = '"view": {"showAppCode": true}'
new_string = '"view": {"showAppCode": false}'

occurrence_count = html_content.count(old_string)
if occurrence_count != 1:
    error_msg = (
        f"ERROR: Expected exactly 1 occurrence of '{old_string}' in HTML, "
        f"but found {occurrence_count} occurrences"
    )
    print(error_msg)
    raise ValueError(error_msg)

html_content = html_content.replace(old_string, new_string)

with open(marimo_html, "w") as f:
    f.write(html_content)

print("Successfully updated showAppCode setting in HTML")
print(f"HTML output written to: {marimo_html}")

log_file.close()
