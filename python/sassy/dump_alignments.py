"""Dump all search_all_alignments() results for a pattern+FASTA to TSV.

Usage::

    python dump_alignments.py --pattern ACGT --fasta refs.fa --k 3 --output out.tsv
"""

import time
from pathlib import Path

import defopt
import pysam

from fgpyo.util.metric import MetricWriter

from sassy import Searcher
from sassy.alignment_record import AlignmentRecord, ContigMetrics

_STRAND_STR: dict[str, str] = {"+": "Fwd", "-": "Rev"}


def dump(
    *,
    pattern: str,
    fasta: Path,
    k: int,
    output: Path = Path("alignments.tsv"),
    prune_suboptimal: bool = False,
) -> None:
    """Dump all alignments for a pattern against every reference in a FASTA file.

    :param pattern: DNA pattern to search for (forward strand only).
    :param fasta: Path to input FASTA file.
    :param k: Maximum edit distance.
    :param output: Output TSV file path.
    :param prune_suboptimal: If set, prune suboptimal alignments before writing.
    """
    # _STRAND_STR maps raw sassy strand symbols to human-readable names;
    # the "-"/"Rev" entry is reserved for when rc=True is enabled.
    searcher = Searcher("iupac")
    records: list[AlignmentRecord] = []
    metrics_path = output.with_name(output.stem).with_suffix(".metrics.tsv")

    with MetricWriter(metrics_path, ContigMetrics) as metrics_writer:
        with pysam.FastxFile(str(fasta)) as fh:
            for ref in fh:
                t0 = time.perf_counter()
                groups = searcher.search_all_alignments(
                    pattern.encode(),
                    ref.sequence.encode(),
                    k,
                    # margin=k, # jf_ branch only
                    prune_suboptimal=prune_suboptimal, # td_ branch only
                )

                for group in groups:
                    for m in group:
                        raw_strand = str(m.strand)
                        records.append(
                            AlignmentRecord(
                                ref_name=ref.name,
                                text_start=m.text_start,
                                text_end=m.text_end,
                                cost=m.cost,
                                cigar=str(m.cigar),
                                strand=_STRAND_STR[raw_strand],
                            )
                        )
                elapsed = time.perf_counter() - t0
                metrics_writer.write(ContigMetrics(
                    contig_name=ref.name,
                    contig_length=len(ref.sequence),
                    num_alignments=sum(len(g) for g in groups),
                    runtime_seconds=elapsed,
                ))

    AlignmentRecord.write(output, *records)


if __name__ == "__main__":
    defopt.run(dump)
