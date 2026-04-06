"""Dump all search_all_alignments() results for a pattern+FASTA to TSV.

Usage::

    python dump_alignments.py --pattern ACGT --fasta refs.fa --k 3 --output out.tsv
"""

from pathlib import Path

import defopt
import pysam

from sassy import Searcher
from sassy.alignment_record import AlignmentRecord

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

    with pysam.FastxFile(str(fasta)) as fh:
        for ref in fh:
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

    AlignmentRecord.write(output, *records)


if __name__ == "__main__":
    defopt.run(dump)
