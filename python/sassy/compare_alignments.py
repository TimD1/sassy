"""Compare two alignment TSV dump files group-by-group.

Usage::

    python compare_alignments.py --file-a td.tsv --file-b jf.tsv
"""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

import defopt

from sassy.alignment_record import AlignmentRecord

# (ref_name, text_end) -> set of (text_start, cigar, strand)
Groups = Dict[Tuple[str, int], Set[Tuple[int, str, str]]]

# (key, alignments_only_in_a, alignments_only_in_b, alignments_in_common)
DiffEntry = Tuple[
    Tuple[str, int],
    Set[Tuple[int, str, str]],
    Set[Tuple[int, str, str]],
    Set[Tuple[int, str, str]],
]


@dataclass
class CompareResult:
    # group-level counts
    identical_groups: List[Tuple[str, int]]
    only_a_groups: List[Tuple[str, int]]
    only_b_groups: List[Tuple[str, int]]
    different_groups: List[DiffEntry]
    # alignment-level counts
    alignments_in_common: int
    alignments_only_a: int
    alignments_only_b: int


def build_groups(path: Path) -> Groups:
    """Read a TSV and group alignments by (ref_name, text_end)."""
    groups: Groups = defaultdict(set)
    for r in AlignmentRecord.read(path):
        groups[(r.ref_name, r.text_end)].add((r.text_start, r.cigar, r.strand))
    return dict(groups)


def compare_groups(groups_a: Groups, groups_b: Groups) -> CompareResult:
    """Diff two group dicts and classify every (ref, text_end) key."""
    all_keys = sorted(set(groups_a) | set(groups_b))
    identical_groups: List[Tuple[str, int]] = []
    only_a_groups: List[Tuple[str, int]] = []
    only_b_groups: List[Tuple[str, int]] = []
    different_groups: List[DiffEntry] = []
    alignments_in_common = 0
    alignments_only_a = 0
    alignments_only_b = 0

    for key in all_keys:
        in_a = key in groups_a
        in_b = key in groups_b
        if in_a and not in_b:
            only_a_groups.append(key)
            alignments_only_a += len(groups_a[key])
        elif in_b and not in_a:
            only_b_groups.append(key)
            alignments_only_b += len(groups_b[key])
        elif groups_a[key] == groups_b[key]:
            identical_groups.append(key)
            alignments_in_common += len(groups_a[key])
        else:
            common = groups_a[key] & groups_b[key]
            a_only = groups_a[key] - groups_b[key]
            b_only = groups_b[key] - groups_a[key]
            different_groups.append((key, a_only, b_only, common))
            alignments_in_common += len(common)
            alignments_only_a += len(a_only)
            alignments_only_b += len(b_only)

    return CompareResult(
        identical_groups=identical_groups,
        only_a_groups=only_a_groups,
        only_b_groups=only_b_groups,
        different_groups=different_groups,
        alignments_in_common=alignments_in_common,
        alignments_only_a=alignments_only_a,
        alignments_only_b=alignments_only_b,
    )


def _print_report(result: CompareResult, label_a: str, label_b: str) -> None:
    W = 38
    print(f"{'Groups':<{W}} {'Count':>6}")
    print("─" * (W + 8))
    print(f"{'  Identical':<{W}} {len(result.identical_groups):>6}")
    print(f"{'  Different':<{W}} {len(result.different_groups):>6}")
    print(f"{'  Only in ' + label_a:<{W}} {len(result.only_a_groups):>6}")
    print(f"{'  Only in ' + label_b:<{W}} {len(result.only_b_groups):>6}")
    print()
    print(f"{'Alignments':<{W}} {'Count':>6}")
    print("─" * (W + 8))
    print(f"{'  In common':<{W}} {result.alignments_in_common:>6}")
    print(f"{'  Only in ' + label_a:<{W}} {result.alignments_only_a:>6}")
    print(f"{'  Only in ' + label_b:<{W}} {result.alignments_only_b:>6}")

    for key, a_only, b_only, common in result.different_groups:
        ref, end = key
        print(
            f"\nref={ref}  text_end={end}  "
            f"[common={len(common)}  "
            f"{label_a}_only={len(a_only)}  {label_b}_only={len(b_only)}]"
        )
        for ts, cig, strand in sorted(a_only):
            print(f"  only in {label_a}: text_start={ts} cigar={cig} strand={strand}")
        for ts, cig, strand in sorted(b_only):
            print(f"  only in {label_b}: text_start={ts} cigar={cig} strand={strand}")

    for ref, end in result.only_a_groups:
        print(f"\nref={ref}  text_end={end}  [group only in {label_a}]")

    for ref, end in result.only_b_groups:
        print(f"\nref={ref}  text_end={end}  [group only in {label_b}]")


def compare(
    *,
    file_a: Path,
    file_b: Path,
    label_a: str = "A",
    label_b: str = "B",
) -> None:
    """Compare two alignment dump TSV files and report differences.

    :param file_a: First TSV file (output of dump_alignments.py).
    :param file_b: Second TSV file (output of dump_alignments.py).
    :param label_a: Label for file_a in the report.
    :param label_b: Label for file_b in the report.
    """
    groups_a = build_groups(file_a)
    groups_b = build_groups(file_b)
    result = compare_groups(groups_a, groups_b)
    _print_report(result, label_a, label_b)


if __name__ == "__main__":
    defopt.run(compare)
