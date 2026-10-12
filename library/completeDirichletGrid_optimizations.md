# `completeDirichletGrid` performance work — changelog

This documents all changes made in this session related to `AdaptiveSparseGrid::completeDirichletGrid()`
and its subfunctions. The goal throughout was to add **timing visibility** and **faster drop-in
alternatives** for the slow steps, without ever modifying the original functions in place — so the
original pipeline still exists unchanged as a correctness/performance baseline to compare against.

`completeDirichletGrid()` (`source/sgrid/sparseGrid.cc`) does:
```cpp
completeGridWithoutBoundary();
CompleteToLocalTensorProductGrid();
completeGridWithoutBoundary();
addHangingNodes();
addPoints2(*this);
```

## 1. Timing instrumentation

**Added** `AdaptiveSparseGrid::completeDirichletGrid_timed()`
(`source/sgrid/sparseGrid.h`, `source/sgrid/sparseGrid.cc`)

Identical in behavior to `completeDirichletGrid()`, but wraps each of the five internal calls with
`std::chrono::high_resolution_clock` timestamps and prints a per-step breakdown (ms) plus a total. This
was the tool used to find every bottleneck described below. `#include <chrono>` was added to
`sparseGrid.cc` for this.

## 2. `CompleteToLocalTensorProductGrid` — fixed an O(N) full-grid rescan

**Problem found:** `LocalForTensor::findAnichtWithDepth` (`source/sgrid/localTensorProduct.cc`) scanned
the *entire* grid (`for i = iNow+1 .. maximalOccupiedSecondTable`), recomputing `Depth(grid->getIndexOfTable(i))`
for every candidate, just to find the next point at a specific depth `T`. It's called once per connected
component at every distinct depth in the grid, so cost was roughly O(depths × components × N).

**Fix:** the grid already maintains a `MultiDepthHashGrid` that buckets every point by exact depth
(`getMultiDepthHashGrid()->getGridForDepth(T)._mapPosToGridPos`, always sorted ascending as a side
effect of how points are appended). Added:

- `LocalForTensor::findAnichtWithDepthFast(iNow, T)` (`localTensorProduct.h`/`.cc`) — binary-searches
  (`std::upper_bound`) into that depth's own bucket instead of scanning the whole grid, and needs no
  `Depth` recomputation per candidate (every entry in a depth's bucket already shares that depth).
  Special-cases the `iNow == (unsigned long)-1` "start of search" sentinel explicitly (the original
  relied on unsigned-integer wraparound for this, which `upper_bound` can't reproduce).
- `AdaptiveSparseGrid::CompleteToLocalTensorProductGrid_fast()` (`sparseGrid.h`/`localTensorProduct.cc`)
  — exact copy of `CompleteToLocalTensorProductGrid()`, calling the fast finder instead.

Complexity per search call: O(N) → O(log N_T + bucket length), where N_T is the (typically much
smaller) point count at that one depth.

## 3. `addHangingNodes` — fixed a genuine D× redundant-loop bug

**Problem found:** `addHangingNodes()` (`source/sgrid/sparseGrid.cc`) had:
```cpp
for (int d = 0; d < DimensionSparseGrid; ++d) {
    for (MultiDimFiveCompass mc; mc.goon(); ++mc) { ... }
}
```
`MultiDimFiveCompass` already enumerates all `5^DimensionSparseGrid` neighbor offsets across every
dimension at once — the `nextFive(mc, Tlocal)` overload used here takes no dimension argument, so `d` is
never used inside the loop body. The entire `5^D` inner sweep was executed `D` times, doing identical
work each time. Confirmed (by grepping every other `MultiDimFiveCompass` call site in the codebase) that
no other function has this pattern — isolated bug.

**Fix:** added `AdaptiveSparseGrid::addHangingNodes_fast()` (`sparseGrid.h`/`sparseGrid.cc`) — identical
logic with the outer `d` loop removed, leaving a single `5^D` sweep per active point.

## 4. Fast pipeline wiring

Added `AdaptiveSparseGrid::completeDirichletGrid_fast()` and `completeDirichletGrid_fast_timed()`
(`sparseGrid.h`/`sparseGrid.cc`) — drop-in replacements for `completeDirichletGrid()` /
`completeDirichletGrid_timed()` that use `CompleteToLocalTensorProductGrid_fast()` and
`addHangingNodes_fast()` in place of the originals. `addPoints2(*this)` is unchanged (still called
directly — it's a free function, not a member).

## 5. `addPoints2` — attempted, no real gain found, reverted

**Problem:** same shape as before — `addPoints2()` (`source/sgrid/GridGeneration.cpp`) loops over every
pair of distinct depths present in the grid, and for pairs differing in ≥2 dimensions, a further nested
loop over every active point at each of those two depths (O(depth-pairs × N₁×N₂)), calling `refill()`
for each pair, all inside a `while(changed)` fixed-point loop that reruns until a full pass adds nothing.
On the current D=6 test grid (729 points, 64 depths), this took ~23-26s — completely dominating
everything else (which totals well under 300ms combined).

**Two ideas tried, both added as `_fast` siblings (`refill_fast`, `addPoints2_fast`) and measured — then
reverted when neither helped:**

1. *Hoist the depth-pair/active-point-pair enumeration out of the `while(changed)` loop.* Verified that
   `refill()` only ever adds new points as hanging (`AddPoint(I, false)`) and never flips an existing
   point's active flag, so the candidate list is provably invariant across passes and can be built once.
   **Measured result: no meaningful difference** (25.2s vs. 24.8s baseline, within run-to-run noise) —
   meaning the `while(changed)` loop converges in very few passes on this grid, so there was little
   redundant enumeration to actually save.
2. *Merge `refill()`'s two `RectangularIteratorExact` box scans into one*, collecting missing points into
   a `std::vector` instead of scanning the box twice. **Measured result: net slower** (38s vs. 26s) —
   the per-call heap allocation for the collection vector cost more than the redundant-but-allocation-free
   second array scan it was meant to eliminate, since most boxes are tiny.

**Conclusion:** the real cost is the sheer number of `(I1, I2)` pairs being enumerated and passed through
`overlap()` — O(N₁×N₂) summed over qualifying depth pairs — which neither change actually reduced.
`addPoints2()` was left completely unmodified; `refill_fast`/`addPoints2_fast` were added, measured,
found not to help, and removed again. Since the pair count itself couldn't be cut without a riskier
rewrite of the core overlap logic, the next attempt used more hardware instead of doing less work — see
§6.

## 6. `addPoints2` — parallelized with OpenMP, ~4.6× measured

Since serial micro-optimizations didn't reduce the O(N₁×N₂) pair-enumeration cost, and 8 cores were
available, the actual win came from parallelizing that enumeration rather than trying to shrink it.

`refill()`'s read-side work (`overlap()`, the box scan, `occupied()` lookups) is embarrassingly parallel
and touches no shared mutable state — the only unsafe part is `AddPoint()`, which mutates the grid's
shared hash tables. The codebase already half-attempted this (`addPoints2_NEU`/`refill_NEU`, never wired
into any pipeline): it defers `AddPoint` into a `pointsToAdd` vector instead of calling it inside the
parallel region — the right idea — but its `#pragma omp parallel for` writes to that vector via
`pointsToAdd.push_back(I)` from every thread with **no synchronization**, a genuine data race.

**Added** `addPoints2_omp()` (`source/sgrid/GridGeneration.cpp`/`.h`), reusing the existing (correct in
isolation) `refill_NEU` but fixing the race: builds the candidate `{I1, I2, T, p, q}` list once (same
invariant-enumeration hoisting as §5.1 — harmless on its own, and here it's also what's needed to get a
flat list to parallelize over), then for each `changed` pass runs a `#pragma omp parallel` region where
each thread accumulates into its **own indexed slot** of a
`std::vector<std::vector<IndexDimension>> perThreadPoints(omp_get_max_threads())` — no shared vector, no
critical section — over a `#pragma omp for schedule(dynamic)` loop calling `refill_NEU`. All `AddPoint`
calls happen only after the parallel region closes, applied serially from the merged per-thread buffers.

**Behavioral nuance:** deferring every `AddPoint` until after a full parallel sweep changes the iteration
from Gauss-Seidel-like (original: a `refill` call can see points added earlier in the *same* pass) to
Jacobi-like (a pass only sees points added in *previous* passes). Since the algorithm is monotonic
(points only ever added, never removed; a box's "already has a depth-T point" check only becomes more
true over time, never less), the final fixed point is provably the same — this can only mean a few more
`changed`-passes to converge, not a different final grid. Confirmed empirically: `grid_omp` matches
`grid` exactly (`gridsAreSimilar` true, same 729 active / 0 hanging point counts).

Added `AdaptiveSparseGrid::completeDirichletGrid_omp()` / `completeDirichletGrid_omp_timed()`
(`sparseGrid.h`/`.cc`) — identical to `completeDirichletGrid_fast()`/`_fast_timed()` (so the comparison
isolates just this change on top of the already-validated `CompleteToLocalTensorProductGrid_fast` and
`addHangingNodes_fast`) with `addPoints2_omp(*this)` swapped in for `addPoints2(*this)`.

**Measured result (8 cores, same D=6 729-point test grid): `addPoints2_omp` 4.94s vs. `addPoints2`
22.8-23.4s — ~4.6× speedup.** Total `completeDirichletGrid_omp_timed`: 5.2s vs. ~23.5-23.7s for the other
two pipelines.

## What was intentionally left unchanged

- `completeDirichletGrid()`, `CompleteToLocalTensorProductGrid()`, `addHangingNodes()` — untouched, kept
  as the baseline for correctness (`gridsAreSimilar`, already in `tests/gridgeneration.cpp`) and
  performance comparisons.
- `addPoints2()` (`source/sgrid/GridGeneration.cpp`) — untouched, kept as the baseline. Its parallel
  sibling `addPoints2_omp()` (§6) is ~4.6× faster and is what `completeDirichletGrid_omp[_timed]` use;
  `completeDirichletGrid_fast[_timed]` still call the original serial `addPoints2()`.

## Measured results

On the D=2 test grid (2 seed points, depth 6 in 2 dimensions):
- `CompleteToLocalTensorProductGrid_fast` vs. original: ~35–46% faster.

On the D=6 test grid (2 seed points, depth 6 in dimensions 0/1 — an earlier version of the test file):
- `addHangingNodes` 37.5 ms → `addHangingNodes_fast` 6.0 ms (~6.25×, matching the `D=6` factor).
- Total `completeDirichletGrid`: 54.2 ms → 20.6 ms (~2.6×).

On a larger D=6 test grid (2 seed points, depth 6 in dimensions 0/1, giving 119 active + 210 hanging
points):
- `addHangingNodes` 105.4 ms → `addHangingNodes_fast` 18.3 ms (~5.8×).
- Total: 194.0 ms → 110.7 ms (~1.75×), with `addPoints2` (unoptimized) now the largest remaining cost.

Note: a subsequent test-file variant using seed points refined in *all six* dimensions simultaneously
(`leftCenter`/`rightCenter`, `for d in 0..5: leftSon(d)` repeated) turned out to be dramatically more
expensive — `completeDirichletGrid` did not finish within 10 minutes on that case. This is consistent
with `addPoints2`'s complexity being sensitive to how many distinct depths and how many dimensions are
simultaneously refined, and is a strong signal that `addPoints2` is worth investigating next if such
grids need to be practical.

## 7. `refill`/`refill_NEU` — merged the two box scans, per-depth hash lookup (~7% on `_omp`)

Follow-up session, still targeting `addPoints2`/`addPoints2_omp` as the dominant cost (this
is exactly what `literature/DissertationScherner.pdf` formalizes as Definition 3, "lokal
volles Prewavelet-Gitter" / `D_(p,q)^full`, and Definition 4's fixed-point completion
criterion — the dissertation states the definition but not an implementation algorithm, so
`refill`/`refill_NEU`'s brute-force per-pair box enumeration is this library's own choice).

**Profiling finding:** instrumenting `refill()` (temporary, removed again) on the test grid
that was current at that point (single `leftSon`/`rightSon` per dimension, 729 points) showed
the cost is dominated by the two `RectangularIteratorExact` box scans inside `refill()`, not
by the number of `(I1,I2)` pairs or `overlap()` calls: 220,408 `refill()` calls but 711M
iterations in the "does this box already contain an occupied point" check-loop and 2.7B
iterations in the "add every missing point" loop, with `addNodes` true on *every single
call* — meaning the entire check-loop was pure duplicate work of what the add-loop redid
anyway. Each iteration calls `AdaptiveSparseGrid_Base::occupied()`, a lookup into a
~1.5M-slot table (`sparseGrid.cc:46`), consistent with the ~26s observed for `addPoints2` on
that grid.

**Added** `refill_v2()`/`refill_NEU_v2()` (`source/sgrid/GridGeneration.cpp`/`.h`), wired
into new `addPoints2_v2()`/`addPoints2_omp_v2()` and
`completeDirichletGrid_v2[_timed]`/`completeDirichletGrid_omp_v2[_timed]`
(`sparseGrid.h`/`.cc`), following the same never-touch-the-baseline convention as §5/§6.
Two changes, both semantics-preserving (verified via `gridsAreSimilar`):

1. **Single-pass box scan.** Instead of a check-loop then a separate add-loop, one
   `RectangularIteratorExact` pass buffers missing cells into a local `pending` vector until
   the first occupied cell is seen; from then on missing cells are added immediately. If the
   pass ends without ever finding an occupied cell, `pending` is discarded and nothing is
   added — identical result to the original two-pass version, computed in one walk of the
   box instead of two.
2. **Per-depth hash lookup.** `occupied()` checks go through
   `grid.getMultiDepthHashGrid()->getGridForDepth(T)` (fetched once per `refill_v2` call,
   since `T` is fixed for the whole box) instead of the whole-grid table. This is exact,
   not approximate: every candidate cell already satisfies `Tneu == T` before the call, and
   every point added to the grid is unconditionally registered in its own depth's
   `SingleDepthHashGrid` too (`setIndexInTable()` → `multiDepthHashGrid->addPoint()`,
   `sparseGrid.cc:494`). Checked first (via a throwaway diagnostic over
   `MultiDepthHashGrid::getAllGrids()`/`SingleDepthHashGrid::getNumberOfEntries()`) that this
   wouldn't backfire: `SimpleMultiHash`'s primary table is a **fixed 100 slots, never
   resized** (`simpleMultiHash.h:231`), so a depth bucket bigger than that would see long
   collision chains; measured max 32 / avg 18.4 entries per depth bucket (256 depths) on the
   current test grid, comfortably clear.
   `addPoints2_omp_v2()` additionally pre-creates every candidate's `SingleDepthHashGrid`
   *serially*, before the parallel region — `getGridForDepth()` inserts into a shared
   `unordered_multimap` on a miss, which is not safe to call concurrently, and the fixed
   candidate list (built once, same invariant as §6) already contains every `T` any pass will
   ever need.

**Measured result — smaller than the profiling above suggested, and workload-dependent.**
The test file's seed points changed (uncommitted, independent of this work) between the
profiling run and this measurement — from one `leftSon(d)`/`rightSon(d)` per dimension to
three, giving a materially different grid (1169 active + 3552 hanging over 256 depths,
instead of 729 active over 64 depths). On *this* grid:
- `addPoints2_omp_v2` 2645 ms vs. `addPoints2_omp` 2830-2870 ms — **~7% faster**.
- `addPoints2_v2` 3416 ms vs. `addPoints2` 3298-3321 ms — **no improvement** (within
  run-to-run noise, possibly marginally slower).

Both `grid_v2` and `grid_omp_v2` are confirmed `gridsAreSimilar` to the untouched baseline.

**Why smaller than expected:** the single-pass merge only saves work on calls where
`addNodes` ends up `true` (original: check-loop + add-loop, both full box scans; merged: one
scan). When `addNodes` is `false`, both versions do exactly one full box scan — there is
nothing to merge away. The profiled 711M/2.7B numbers came from a grid where `addNodes` was
`true` 100% of the time; this grid's `(T1,T2)` depth-pair/point-pair mix apparently produces
far fewer confirmed hits relative to attempted ones, so the merge has much less redundant
work to eliminate. The per-depth hash lookup change is a straight, low-risk win (one fewer
`Depth` recomputation per lookup than the global table's `hash()`, which recomputes
`Depth(index)` internally — `sparseGrid.h:416` vs. `SimpleMultiHash` using its stored `depth`
member) but a small one when the global 1.5M-slot table was already near-collision-free at
these point counts.

**Left as a follow-up, not implemented:** the same bucket-size measurement that cleared
Step 2 (max 32 entries per depth bucket) suggests a bigger, more targeted fix for the
`addNodes == false`-heavy regime this grid apparently hits: determine "does this box contain
an occupied depth-`T` point" by iterating the small depth-`T` bucket itself and testing each
entry against `[Imin, Imax]` (`O(bucket size)`), instead of walking the (potentially much
larger) box and hash-probing every cell (`O(box volume)`) as `refill_v2` still does. Since
`T` is fixed across an entire `(T1, T2)` depth-pair's whole `i,j` double loop, the same small
bucket would be reused across every pair in that loop. This changes the check phase's
algorithmic shape (not just its constant factor) and wasn't attempted this session.

## Files touched

- `source/sgrid/sparseGrid.h` — declarations for `completeDirichletGrid_timed`,
  `CompleteToLocalTensorProductGrid_fast`, `addHangingNodes_fast`, `completeDirichletGrid_fast`,
  `completeDirichletGrid_fast_timed`, `completeDirichletGrid_omp`, `completeDirichletGrid_omp_timed`.
- `source/sgrid/sparseGrid.cc` — `#include <chrono>`, plus definitions of all the above.
- `source/sgrid/localTensorProduct.h` — declaration of `findAnichtWithDepthFast`.
- `source/sgrid/localTensorProduct.cc` — `#include <algorithm>`, plus definitions of
  `findAnichtWithDepthFast` and `CompleteToLocalTensorProductGrid_fast`.
- `source/sgrid/GridGeneration.h` — declaration of `addPoints2_omp`.
- `source/sgrid/GridGeneration.cpp` — definition of `addPoints2_omp` (reuses the existing `refill_NEU`
  unchanged).
- `tests/gridgeneration.cpp` — demo/benchmark harness exercising all three pipelines side by side via
  `gridsAreSimilar` for correctness and the `_timed` variants for performance (also independently edited
  by the user during this session to change seed-point placement and dimensionality).

§7 additions:
- `source/sgrid/GridGeneration.h` — declarations for `refill_v2`, `refill_NEU_v2`,
  `addPoints2_v2`, `addPoints2_omp_v2`.
- `source/sgrid/GridGeneration.cpp` — definitions of the above; `refill`, `refill_NEU`,
  `addPoints2`, `addPoints2_omp` remain unchanged.
- `source/sgrid/sparseGrid.h`/`.cc` — `completeDirichletGrid_v2[_timed]`,
  `completeDirichletGrid_omp_v2[_timed]`, mirroring `_fast`/`_omp`.
- `tests/gridgeneration.cpp` — added `grid_v2`/`grid_omp_v2` blocks alongside the existing
  three, same `gridsAreSimilar` + `_timed` pattern.
