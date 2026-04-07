# max_gaps parameter for iterate_all_alignments

## Context

`iterate_all_alignments` enumerates all alignment paths of cost ≤ `k` via DFS over the DP cost
matrix. With no gap limit, the number of paths can be large: any combination of insertions and
deletions within budget is explored. Adding a `max_gaps: Option<usize>` parameter lets callers
bound the number of individual gap operations (each `I` or `D` CIGAR op counts as one gap),
pruning paths early and reducing total paths checked.

## Definition

A **gap** is a single `I` (insertion) or `D` (deletion) CIGAR operation. Consecutive indels each
count separately. `max_gaps: None` means no limit (current behaviour).

## API Changes

### `iterate_all_alignments` (public)

```rust
pub fn iterate_all_alignments(
    &self,
    pattern: &[u8],
    text: &[u8],
    k: usize,
    matches: &[Match],
    partial_matches: bool,
    prune_suboptimal: bool,
    max_gaps: Option<usize>,   // ← new
    callback: &mut impl Callback,
)
```

### `search_all_alignments` (public)

```rust
pub fn search_all_alignments(
    &mut self,
    pattern: &[u8],
    text: &[u8],
    k: usize,
    prune_suboptimal: bool,
    max_gaps: Option<usize>,   // ← new
) -> Vec<Vec<Match>>
```

### `iterate_one_strand` (private)

Gains the same `max_gaps: Option<usize>` parameter, propagated from `iterate_all_alignments`.

### `Context` struct

Two new fields:

```rust
max_gaps: Option<usize>,
gaps: usize,   // gap count on the active DFS path; reset to 0 per text_end iteration
```

## DFS Changes (`dfs()`)

**Guard in the edge loop** (alongside the existing leading/trailing deletion guard):

```rust
if (op == CigarOp::Ins || op == CigarOp::Del)
    && self.max_gaps.is_some_and(|mg| self.gaps >= mg)
{
    continue;
}
```

**Increment/decrement around the recursive call** (alongside the existing cost mutation/revert):

```rust
let is_gap = matches!(op, CigarOp::Ins | CigarOp::Del);
if is_gap { self.gaps += 1; }
let continuation = self.dfs::<P>();
if is_gap { self.gaps -= 1; }
```

`gaps` is reset to `0` at the start of each `text_end` iteration in `iterate_one_strand`,
consistent with how `m.cost`, `m.pattern_start`, etc. are reset.

## Existing Callers

All existing callers of `iterate_all_alignments` and `search_all_alignments` pass `None` for
`max_gaps` to preserve current behaviour.

## Testing

New test in `src/search.rs`:

```rust
fn max_gaps_excludes_high_gap_alignments() {
    // pattern=ACGT text=AACGT: exact match at [1..5] (0 gaps),
    // plus cost-1 alignments with 1 gap (e.g. 1I3=).
    // max_gaps=0 should keep only gap-free alignments.
    let pattern = b"ACGT";
    let text = b"AACGT";
    let k = 1;
    let groups_limited =
        Searcher::<Dna>::new(false, None).search_all_alignments(pattern, text, k, false, Some(0));
    let groups_unlimited =
        Searcher::<Dna>::new(false, None).search_all_alignments(pattern, text, k, false, None);

    for group in &groups_limited {
        for m in group {
            let cigar = m.cigar.to_string();
            assert!(!cigar.contains('I') && !cigar.contains('D'),
                "gap found with max_gaps=0: {cigar}");
        }
    }
    assert!(
        groups_limited.iter().map(|g| g.len()).sum::<usize>()
            < groups_unlimited.iter().map(|g| g.len()).sum::<usize>(),
        "max_gaps=0 should exclude at least one alignment"
    );
}
```
