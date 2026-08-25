"""Count public HA sequences by mutation-count category vs. library strains."""

import io
import multiprocessing
import sys
import tarfile
from collections import defaultdict

import pandas as pd
import regex
from Bio import SeqIO
from mutation_categories import category_label

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

CANONICAL_AAS = set("ACDEFGHIKLMNPQRSTVWY")

max_mutations_computed = snakemake.params.max_mutations_computed
time_bin_freq = snakemake.params.time_bin_freq
regions = snakemake.params.regions
gisaid_fasta_members = snakemake.params.gisaid_fasta_members
n_workers = snakemake.threads

_region_patterns = None
_max_mutations_computed = None


def _init_worker(region_patterns, max_mutations):
    global _region_patterns, _max_mutations_computed
    _region_patterns = region_patterns
    _max_mutations_computed = max_mutations


def _match_one(seq_and_period):
    """Match one query sequence against every region's library patterns.

    Module-level (not nested) so it and its `_region_patterns`/
    `_max_mutations_computed` globals (set once per worker via the pool
    initializer) are usable by `multiprocessing.Pool` workers.
    """
    seq, period_start = seq_and_period
    result = {}
    for region, patterns in _region_patterns.items():
        best_ndiff = None
        for ref, pattern in patterns.values():
            if ref in seq:
                best_ndiff = 0
                break  # can't do better than an exact match
            m = pattern.search(seq)
            if m is None:
                continue
            ndiff = sum(m.fuzzy_counts)
            if best_ndiff is None or ndiff < best_ndiff:
                best_ndiff = ndiff
        result[region] = category_label(best_ndiff, _max_mutations_computed)
    return period_start, result


lib = pd.read_csv(
    snakemake.input.library_csv,
    usecols=["strain", "subtype", "protein_sequence_HA_ectodomain"],
    dtype=str,
)

unexpected_subtypes = set(lib["subtype"]) - set(gisaid_fasta_members)
if unexpected_subtypes:
    raise ValueError(
        f"library has subtypes not in config['regions']: {unexpected_subtypes}"
    )

rows = []

for subtype in gisaid_fasta_members:
    print(f"\n=== {subtype} ===")

    # one reference sequence per named strain (collapses barcode replicates)
    strain_seqs = {}
    for strain, group in lib.loc[lib["subtype"] == subtype].groupby("strain"):
        seqs = group["protein_sequence_HA_ectodomain"].str.upper().unique()
        if len(seqs) != 1:
            raise ValueError(f"strain {strain!r} has inconsistent library sequences")
        (seq,) = seqs
        if set(seq) - CANONICAL_AAS:
            raise ValueError(f"library strain {strain!r} has non-canonical amino acids")
        strain_seqs[strain] = seq
    print(f"{len(strain_seqs)} unique library strains")

    # per region: trimmed reference + a `regex` fuzzy pattern per strain. The `(?b)`
    # (BESTMATCH) mode finds the lowest-edit-distance match of the reference as a
    # substring anywhere in the (untrimmed, unaligned) query sequence, so no prior
    # alignment of query to library numbering is needed.
    region_patterns = {}
    for region, bounds_by_subtype in regions.items():
        start, end = bounds_by_subtype[subtype]
        assert (
            1 <= start < end
        ), f"invalid region bounds for {subtype} {region}: {start}-{end}"
        patterns = {}
        for strain, seq in strain_seqs.items():
            assert len(seq) >= end, f"{strain!r} shorter than region {region} end {end}"
            ref = seq[start - 1 : end]
            pattern = regex.compile(
                f"(?b)(?:{regex.escape(ref)}){{e<={max_mutations_computed}}}"
            )
            patterns[strain] = (ref, pattern)
        region_patterns[region] = patterns

    n_invalid_aa = 0
    n_dup_accession = 0
    seen_accessions = set()
    records_to_match = []  # (seq, period_start) per valid, deduplicated record

    member = gisaid_fasta_members[subtype]
    with (
        tarfile.open(snakemake.input.gisaid_tar, mode="r:xz") as tf,
        io.TextIOWrapper(tf.extractfile(member)) as text_stream,
    ):
        for record in SeqIO.parse(text_stream, "fasta"):
            seq = str(record.seq).upper().rstrip("*")
            if set(seq) - CANONICAL_AAS:
                n_invalid_aa += 1
                continue

            fields = record.description.split("|")
            if len(fields) != 4:
                raise ValueError(f"expected 4 pipe-delimited header fields: {fields}")
            accession, _, _, date_str = fields
            date = pd.Timestamp(date_str)  # raises on a malformed/incomplete date

            if accession in seen_accessions:
                n_dup_accession += 1
                continue
            seen_accessions.add(accession)

            period_start = date.to_period(time_bin_freq).start_time.date()
            records_to_match.append((seq, period_start))

    n_kept = len(records_to_match)
    print(f"dropped {n_invalid_aa} sequences with non-canonical amino acids")
    print(f"dropped {n_dup_accession} duplicate accessions")
    print(f"kept {n_kept} sequences for matching, using {n_workers} worker processes")

    # the matching loop is embarrassingly parallel across query sequences. Use an
    # explicit "fork" context (not this Python's "forkserver" default): workers
    # must inherit the already-built `region_patterns`/data via a copy-on-write
    # fork, rather than reimporting/re-running this whole top-level script from
    # scratch to reconstruct that state, which is what "forkserver"/"spawn" would
    # otherwise require.
    counts = defaultdict(int)  # (region, period_start, category) -> n_sequences
    with multiprocessing.get_context("fork").Pool(
        processes=n_workers,
        initializer=_init_worker,
        initargs=(region_patterns, max_mutations_computed),
    ) as pool:
        for period_start, region_categories in pool.imap_unordered(
            _match_one, records_to_match, chunksize=200
        ):
            for region, category in region_categories.items():
                counts[(region, period_start, category)] += 1

    # every kept sequence must be assigned to exactly one category per region
    for region in region_patterns:
        region_total = sum(n for (r, _, _), n in counts.items() if r == region)
        assert region_total == n_kept, (
            f"{subtype} {region}: counts sum to {region_total}, "
            f"expected {n_kept} kept sequences"
        )

    for (region, period_start, category), n_sequences in counts.items():
        rows.append(
            {
                "subtype": subtype,
                "region": region,
                "week_start": period_start,
                "category": category,
                "n_sequences": n_sequences,
            }
        )

df = (
    pd.DataFrame(rows)
    .sort_values(["subtype", "region", "week_start", "category"])
    .reset_index(drop=True)
)
df.to_csv(snakemake.output.counts_csv, index=False, float_format="%.5g")
