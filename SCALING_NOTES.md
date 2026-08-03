# Scaling notes

This file records reproducible performance observations from the first search beyond
the 200-digit challenge entry.  Search manifests and candidate files are intentionally
not committed; the commands below are sufficient to reproduce their structure.

## Test system

- Apple M4, 10 physical cores, 16 GB RAM
- PARI/GP 2.17.4 with `seadata`
- GMP-ECM 7.0.6
- macOS, 3 August 2026

The target used below is the least prime greater than `10^210`, the first line of
`targets.txt`.

## Smooth-part extraction

`smoothpart` now caches the square-free product of every prime through the certificate
bound and repeatedly takes gcds with the order.  On ten representative 701-bit orders,
the old prime-by-prime loop took 53 ms and the cached-gcd path took 4 ms, a 13.2x
screening speedup.  Building the 701-bit target's cache took 128 ms; its product was
701,527 bits (about 86 KiB).  The old and new paths returned identical smooth and rough
parts.

The cache is shared by every root checkpoint attempted in one GP process.  It is
therefore especially useful with `--resume-per-worker`, which avoids rebuilding both
the cache and PARI's private factor table between candidate orders.

## Two-pass CM search

Before invoking Cornacchia, the scanner now tests the necessary quadratic-residue
condition.  A representation `4p=t^2+Dv^2` implies that `-D` is a square modulo
`p` (and the corresponding statement holds for the even-discriminant form), so a
Kronecker-symbol rejection is exact.  On the first 50,000 odd/even discriminant
tests at the 210-digit target, this reduced Cornacchia calls to 15,229 and wall
time from 22.65 seconds to 10.03 seconds, a 2.26x speedup, while preserving all
239 represented orders.

The production enumeration used no SEA calls and no factoring:

```sh
read -r p < targets.txt
python3 parallel_short.py "$p" -j 10 --seed 2113000 \
  --cm-bound 100000000 --cm-screen-only \
  --candidate-bits 40 --candidate-dir candidates/210 \
  --manifest search-runs-210.jsonl
```

It completed in 1,264.6 seconds of wall time and 9,634.9 aggregate child CPU seconds.
The run saved 1,095 distinct live orders before later factoring and tombstones changed
the frontier, reached an 86-bit smooth part, and made zero SEA calls.  Equal-width
discriminant partitioning kept the finite scan bounded and avoided the severe tail seen
when factor work and representation density were both assigned by square-root ranges.

The expensive second stage ranked the saved orders and assigned twenty to each worker:

```sh
python3 parallel_short.py "$p" -j 10 --seed 2115000 \
  --factor-seconds 20,20,20,20,20,20,5,5,5,5 \
  --factor-flags 0,0,0,0,0,0,9,9,9,9 \
  --factor-rounds 3,3,3,3,3,3,2,2,2,2 \
  --prefactor-seconds 2 --pm1-bound 1000000 --pp1-bound 1000000 \
  --ecm-bound 0,0,0,0,0,0,0,0,100000,100000 \
  --ecm-curves 0,0,0,0,0,0,0,0,8,8 \
  --resume-candidates candidates/210 --resume-top 200 \
  --resume-per-worker 20 --resume-only --candidate-bits 40 \
  --candidate-dir candidates/210 --curve-families 5 \
  --branch-curves 96 --manifest search-runs-210.jsonl -o cert210.txt
```

This pass exhausted its 200 assigned checkpoints in 564.6 seconds and 2,387.0 child
CPU seconds without a certificate.  It nevertheless demonstrated why checkpoints must
retain partial factor progress: one cofactor fell from 502 to 397 bits, and another
from 512 to 444 bits.  Re-ranking these reduced cofactors made them the first two items
in the next portfolio.  In the fast worker profiles, P-1 and P+1 found factors in many
orders; the two GMP-ECM profiles reported 17 factor-bearing batches in 34 attempts.

The practical workflow is consequently an adaptive frontier, not a fixed seed list:

1. enumerate many exact orders with no factor budget;
2. apply cheap P-1, P+1, and short partial-factor passes to the ranked frontier;
3. persist every smaller residual and re-rank;
4. spend independent ECM curves and longer complete-factor limits only on the new
   frontier; and
5. report both wall time and aggregate core-seconds, plus actual `resume_attempts`.

The frontier also applies an exact impossibility bound.  Every unresolved prime factor
is larger than `n^2`; for a composite residual `R`, a future child prime is therefore at
most `R/(n^2+1)`.  Once that upper bound times the known smooth part cannot clear the
certificate window, GP writes a tombstone and the Python loader suppresses every older
checkpoint for the same order.  A focused ECM pass triggered this rule by reducing a
397-bit residual to 310 bits: another split would leave at most 291 bits, below that
order's 299-bit minimum.

## Torsion-conditioned curves

Curve family 5 ports the optimized `X_1(27)` construction from OneShotSEA.  The retained
side has full rational 2-torsion and a point of order four, so its group order is
divisible by 216; the twist remains available for the same point count.  Generation at
701 bits takes a fraction of a second, whereas a PARI SEA call on an unconditioned curve
at this size takes roughly three to four minutes on one core.  A complete two-worker
20-digit smoke search found and verified a certificate with family 5 in 0.3 seconds.

## OneShotSEA boundary at 701 bits

OneShotSEA commit `c0bbf81` was built with its native C++20/GMP path.  Its generic
Weber-compatible curve at seed 2115000, index 0 was counted with one SEA thread.

With the checked-in 77 levels through 401, the run used 210.66 seconds wall and 179.95
seconds user CPU, then correctly stopped at the level limit with about 7.95e22 compatible
traces.  An arbitrary Montgomery `X_1(27)` curve was not eligible for this path because it
had no rational Weber lift.

The pinned archive was then materialized through level 601 (108 authenticated tables,
53 MB).  The second run stopped after level 461 in 557.07 seconds wall and 317.26 seconds
user CPU.  It soundly enumerated the complete remaining set of 67,306 traces, but did not
produce a unique point count.  Thus the current generic Weber path is not a practical
replacement for PARI SEA at 701 bits on this machine.

The useful longer-term route is narrower:

- generate `X_1(27)` curves together with their authenticated Weber coordinate, rather
  than converting an unrelated Montgomery curve after generation;
- feed the exact 216/432 torsion trace prior into SEA;
- apply exact smooth-part rejection to the complete partial trace set before paying for
  more levels; and
- replace finite high-level tables with direct finite-field specialization as the input
  grows.

The CM two-pass search is currently the better practical engine for the 210--300 digit
targets.  The Weber route remains the asymptotically motivated next-round design, but the
measurements above identify the missing integration points instead of assuming a
crossover that has not occurred.
