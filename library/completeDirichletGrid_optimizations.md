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

## What was intentionally left unchanged

- `completeDirichletGrid()`, `CompleteToLocalTensorProductGrid()`, `addHangingNodes()` — untouched, kept
  as the baseline for correctness (`gridsAreSimilar`, already in `tests/gridgeneration.cpp`) and
  performance comparisons.
- `addPoints2()` (`source/sgrid/GridGeneration.cpp`) — **not yet optimized**. It has a nested loop over
  every pair of distinct depths present in the grid, and for depth pairs that differ in ≥2 dimensions, a
  further nested loop over every point at each of those two depths (O(depth-pairs × N₁×N₂)). On the
  D=6 test grid in `tests/gridgeneration.cpp`, this became the single largest cost in the fast pipeline
  (~73 ms out of ~111 ms total) once the two bugs above were fixed. This is the next candidate if
  further optimization is wanted.

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

## Files touched

- `source/sgrid/sparseGrid.h` — declarations for `completeDirichletGrid_timed`,
  `CompleteToLocalTensorProductGrid_fast`, `addHangingNodes_fast`, `completeDirichletGrid_fast`,
  `completeDirichletGrid_fast_timed`.
- `source/sgrid/sparseGrid.cc` — `#include <chrono>`, plus definitions of all the above.
- `source/sgrid/localTensorProduct.h` — declaration of `findAnichtWithDepthFast`.
- `source/sgrid/localTensorProduct.cc` — `#include <algorithm>`, plus definitions of
  `findAnichtWithDepthFast` and `CompleteToLocalTensorProductGrid_fast`.
- `tests/gridgeneration.cpp` — demo/benchmark harness exercising both pipelines side by side via
  `gridsAreSimilar` for correctness and the `_timed` variants for performance (also independently edited
  by the user during this session to change seed-point placement and dimensionality).
