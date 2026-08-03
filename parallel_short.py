#!/usr/bin/env python3
"""Run independent short.gp searches in parallel and verify the first result."""

import argparse
from datetime import datetime, timezone
import fcntl
import json
import math
import os
from pathlib import Path
import re
import resource
import signal
import shutil
import subprocess
import sys
import tempfile
import time
import uuid


ROOT = Path(__file__).resolve().parent
SHORT_GP = ROOT / "short.gp"
VERIFIER = ROOT / "vsmallECPP.py"
COUNTERS_RE = re.compile(
    r"curves=(\d+)\s+sea_calls=(\d+)(?:\s+sea_aborts=(\d+))?"
    r"(?:\s+factor_attempts=(\d+)\s+factor_aborts=(\d+)"
    r"\s+descents=(\d+)\s+backtracks=(\d+)"
    r"(?:\s+ecm_attempts=(\d+)\s+ecm_factors=(\d+)"
    r"(?:\s+q_candidates=(\d+)\s+viable_q=(\d+)\s+window_candidates=(\d+)"
    r"\s+residual_bits=(\d+)\s+max_smooth_bits=(\d+)"
    r"\s+smooth_rejects=(\d+))?)?)?"
    r"(?:\s+deep_factor_attempts=(\d+))?"
    r"(?:\s+pm1_attempts=(\d+)\s+pm1_factors=(\d+)"
    r"\s+pp1_attempts=(\d+)\s+pp1_factors=(\d+))?"
    r"(?:\s+msieve_attempts=(\d+)\s+msieve_factors=(\d+))?"
    r"(?:\s+exhausted_orders=(\d+))?"
    r"(?:\s+factor_recoveries=(\d+))?"
    r"(?:\s+cm_tests=(\d+)\s+cm_last_D=(\d+))?"
    r"(?:\s+resume_attempts=(\d+))?"
)
RESUME_EXHAUSTED_RE = re.compile(r"resume_exhausted=(\d+)")
REAL_TIME_RE = re.compile(r"^real ([0-9.]+)$", re.MULTILINE)
USER_TIME_RE = re.compile(r"^user ([0-9.]+)$", re.MULTILINE)
SYS_TIME_RE = re.compile(r"^sys ([0-9.]+)$", re.MULTILINE)


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def positive_int_list(value):
    try:
        values = [positive_int(item) for item in value.split(",")]
    except (ValueError, argparse.ArgumentTypeError) as error:
        raise argparse.ArgumentTypeError("must be comma-separated positive integers") from error
    if not values:
        raise argparse.ArgumentTypeError("must not be empty")
    return values


def nonnegative_int(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return value


def nonnegative_int_list(value):
    try:
        values = [nonnegative_int(item) for item in value.split(",")]
    except (ValueError, argparse.ArgumentTypeError) as error:
        raise argparse.ArgumentTypeError(
            "must be comma-separated nonnegative integers"
        ) from error
    if not values:
        raise argparse.ArgumentTypeError("must not be empty")
    return values


def int_list(value):
    try:
        values = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError("must not be empty")
    return values


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search for a short ECPP with parallel PARI/GP workers."
    )
    parser.add_argument("p", type=positive_int, help="probable prime to certify")
    parser.add_argument(
        "-j",
        "--workers",
        type=positive_int,
        default=max(1, os.cpu_count() or 1),
        help="number of GP workers (default: number of logical CPUs)",
    )
    parser.add_argument(
        "--gp-threads",
        type=positive_int,
        default=1,
        help="PARI threads allowed per worker (default: 1)",
    )
    parser.add_argument(
        "--factor-seconds",
        type=nonnegative_int_list,
        default=[20, 20, 20, 20, 20, 20, 20, 2, 5, 20],
        help=(
            "comma-separated per-order factorization limits assigned cyclically "
            "to workers; 0 disables PARI factorint without disabling external ECM "
            "(default: 20,20,20,20,20,20,20,2,5,20)"
        ),
    )
    parser.add_argument(
        "--factor-flags",
        type=nonnegative_int_list,
        default=[0, 0, 0, 0, 0, 0, 0, 9, 8, 0],
        help=(
            "comma-separated PARI factorint flags assigned cyclically "
            "to workers (default: 0,0,0,0,0,0,0,9,8,0; "
            "9 is fast partial factoring)"
        ),
    )
    parser.add_argument(
        "--factor-rounds",
        type=positive_int_list,
        default=[3],
        help=(
            "comma-separated maximum timed factor passes; another pass runs "
            "only after the preceding one shrinks its residual (default: 3)"
        ),
    )
    parser.add_argument(
        "--deep-factor-bits",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated residual bit thresholds for adaptive longer "
            "factorization; 0 disables it (default: 0)"
        ),
    )
    parser.add_argument(
        "--deep-factor-seconds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated adaptive factorization limits used below the "
            "corresponding bit threshold, even when --factor-seconds is zero; "
            "0 disables it (default: 0)"
        ),
    )
    parser.add_argument(
        "--prefactor-seconds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated saved partial-factor limits assigned cyclically "
            "to workers; 0 disables the preliminary pass (default: 0)"
        ),
    )
    parser.add_argument(
        "--prefactor-flags",
        type=nonnegative_int_list,
        default=[8],
        help=(
            "comma-separated factorint flags for saved partial-factor passes "
            "(default: 8)"
        ),
    )
    parser.add_argument(
        "--pm1-bounds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated saved GMP P-1 stage-1 bounds assigned "
            "cyclically to workers; 0 disables P-1 (default: 0)"
        ),
    )
    parser.add_argument(
        "--pp1-bounds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated saved GMP P+1 stage-1 bounds assigned "
            "cyclically to workers; 0 disables P+1 (default: 0)"
        ),
    )
    parser.add_argument(
        "--ecm-bounds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated GMP-ECM stage-1 bounds assigned cyclically to "
            "workers; 0 disables external ECM (default: 0)"
        ),
    )
    parser.add_argument(
        "--ecm-curves",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated GMP-ECM curve counts assigned cyclically to "
            "workers; 0 disables external ECM (default: 0)"
        ),
    )
    parser.add_argument(
        "--ecm-rounds",
        type=positive_int_list,
        default=[1],
        help=(
            "comma-separated maximum saved-factor rounds assigned cyclically "
            "to GMP-ECM workers (default: 1)"
        ),
    )
    parser.add_argument(
        "--msieve-bits",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated residual bit thresholds for bounded msieve "
            "passes; 0 disables msieve (default: 0)"
        ),
    )
    parser.add_argument(
        "--msieve-seconds",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated wall-time limits for one-thread msieve passes; "
            "0 disables msieve (default: 0)"
        ),
    )
    parser.add_argument(
        "--curve-families",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated curve families assigned cyclically to workers: "
            "0 random; 1/2/3 known 4/8/16-torsion; 4 known 3-torsion; "
            "5 X_1(27) with a rational point of order 4 (default: 0)"
        ),
    )
    parser.add_argument(
        "--root-smooth-bits",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated minimum root smooth-part bit lengths assigned "
            "cyclically to workers; 0 disables filtering (default: 0)"
        ),
    )
    parser.add_argument(
        "--cm-start",
        type=nonnegative_int,
        default=3,
        help="first CM discriminant magnitude to scan (default: 3)",
    )
    parser.add_argument(
        "--cm-bound",
        type=nonnegative_int,
        default=0,
        help=(
            "scan fundamental CM discriminants through this bound before the "
            "random-SEA search; 0 disables CM presearch (default: 0)"
        ),
    )
    parser.add_argument(
        "--cm-smooth-bits",
        type=nonnegative_int,
        default=40,
        help=(
            "minimum smooth-part bit length before factoring a CM order "
            "(default: 40)"
        ),
    )
    parser.add_argument(
        "--cm-kind",
        choices=("both", "odd", "even"),
        default="both",
        help="CM discriminant types to scan (default: both)",
    )
    parser.add_argument(
        "--cm-partition",
        choices=("sqrt", "width"),
        default="sqrt",
        help=(
            "worker partition: sqrt balances expected represented orders when "
            "factoring inline; width balances discriminants in a screen-only "
            "enumeration pass (default: sqrt)"
        ),
    )
    parser.add_argument(
        "--cm-only",
        action="store_true",
        help=(
            "stop successfully after the finite CM range is exhausted instead "
            "of falling back to an unbounded random-SEA search"
        ),
    )
    parser.add_argument(
        "--cm-screen-only",
        action="store_true",
        help=(
            "enumerate and checkpoint CM orders without factoring; implies "
            "--cm-only, --cm-partition width, and zero factoring limits"
        ),
    )
    parser.add_argument(
        "--candidate-bits",
        type=nonnegative_int_list,
        default=[0],
        help=(
            "comma-separated minimum root smooth-part bit lengths for saving "
            "resumable curve orders; 0 disables saving (default: 0)"
        ),
    )
    parser.add_argument(
        "--candidate-dir",
        type=Path,
        help=(
            "directory for per-worker resumable root-order logs; required when "
            "--candidate-bits is nonzero"
        ),
    )
    parser.add_argument(
        "--resume-candidates",
        type=Path,
        help=(
            "checkpoint file or directory of per-worker checkpoint files; complete "
            "root levels and ranked saved orders are assigned before fresh search"
        ),
    )
    parser.add_argument(
        "--resume-top",
        type=positive_int,
        help="use only the top N ranked saved orders with --resume-candidates",
    )
    parser.add_argument(
        "--resume-offset",
        type=nonnegative_int,
        default=0,
        help="skip this many top-ranked saved orders before assignment (default: 0)",
    )
    parser.add_argument(
        "--resume-per-worker",
        type=positive_int,
        default=1,
        help=(
            "number of ranked checkpoints each worker tries before fallback; "
            "checkpoints are distributed round-robin in rank order (default: 1)"
        ),
    )
    parser.add_argument(
        "--resume-only",
        action="store_true",
        help=(
            "stop successfully after assigned resumed orders are exhausted "
            "instead of starting another CM or random-SEA search"
        ),
    )
    parser.add_argument(
        "--sea-torsions",
        type=int_list,
        default=[0],
        help=(
            "comma-separated ellsea torsion filters assigned cyclically to workers; "
            "0 uses unfiltered ellcard (default: 0)"
        ),
    )
    parser.add_argument(
        "--branch-curves",
        type=nonnegative_int,
        default=64,
        help=(
            "curve-order budget for each child subtree before backtracking "
            "(default: 64; 0 disables backtracking)"
        ),
    )
    parser.add_argument(
        "--branch-curves-list",
        type=nonnegative_int_list,
        help=(
            "comma-separated child-branch budgets assigned cyclically to workers; "
            "overrides --branch-curves"
        ),
    )
    parser.add_argument(
        "--seed",
        type=positive_int,
        default=1,
        help="base random seed; worker i uses seed+i (default: 1)",
    )
    parser.add_argument(
        "--seeds",
        type=positive_int_list,
        help=(
            "explicit comma-separated worker seeds; when supplied, its length "
            "determines the worker count and --seed/-j are ignored"
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="append start and completion records to a JSON Lines run manifest",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="write the verified certificate to this file"
    )
    return parser.parse_args()


def worker_input(
    p,
    factor_seconds,
    factor_flags,
    factor_rounds,
    deep_factor_bits,
    deep_factor_seconds,
    prefactor_seconds,
    prefactor_flags,
    pm1_bound,
    pp1_bound,
    ecm_bound,
    ecm_curves,
    ecm_rounds,
    msieve_bits,
    msieve_seconds,
    curve_family,
    root_smooth_bits,
    candidate_bits,
    candidate_file,
    sea_torsion,
    branch_curves,
    seed,
    gp_threads,
    resume_orders,
    cm_start,
    cm_bound,
    cm_smooth_bits,
    cm_slot,
    cm_slots,
    cm_kind,
    cm_partition,
    cm_only,
    resume_only,
):
    search = "c=0; "
    initialized = False
    for resume_order in resume_orders:
        reset = 0 if initialized else 1
        if resume_order["kind"] == "level":
            search += (
                f"if(c==0,c=shortcertfromlevel({p},{resume_order['A']},"
                f"{resume_order['x']},{resume_order['o']},{resume_order['q']},"
                f"{reset})); "
            )
        elif resume_order["kind"] == "cm_screen":
            search += (
                f"if(c==0,c=shortcertcmfromscreen({p},-{resume_order['D']},"
                f"{resume_order['N']},{resume_order['o']},{resume_order['q']},"
                f"{reset})); "
            )
        elif resume_order["xden"] == "0":
            search += (
                f"if(c==0,c=shortcertcmfromorder({p},-{resume_order['A']},"
                f"{resume_order['N']},{resume_order['residual']},{reset})); "
            )
        else:
            search += (
                f"if(c==0,c=shortcertfromorder({p},{resume_order['A']},"
                f"{resume_order['xden']},{resume_order['N']},0,"
                f"{resume_order['residual']},{reset})); "
            )
        search += (
            "if(c==0,warning(\"resume_exhausted=\",SC_exhaustedorders)); "
        )
        initialized = True
    if cm_bound and not resume_only:
        search += (
            f"if(c==0,c=shortcertcm({p},{cm_slot},{cm_slots},"
            f"{cm_start},{cm_bound},{cm_smooth_bits},{cm_kind},"
            f"{0 if initialized else 1},{cm_partition})); "
        )
        initialized = True
    if not cm_only and not resume_only:
        search += f"if(c==0,c=shortcert({p},{0 if initialized else 1})); "
    return (
        f"default(nbthreads,{gp_threads}); SC_tlim={factor_seconds}; "
        f"SC_factorflags={factor_flags}; "
        f"SC_factorrounds={factor_rounds}; "
        f"SC_deepfactorbits={deep_factor_bits}; "
        f"SC_deepfactortlim={deep_factor_seconds}; "
        f"SC_prefactortlim={prefactor_seconds}; "
        f"SC_prefactorflags={prefactor_flags}; "
        f"SC_pm1b1={pm1_bound}; "
        f"SC_pp1b1={pp1_bound}; "
        f"SC_ecmb1={ecm_bound}; "
        f"SC_ecmcurves={ecm_curves}; "
        f"SC_ecmrounds={ecm_rounds}; "
        f"SC_msievebits={msieve_bits}; "
        f"SC_msieveseconds={msieve_seconds}; "
        f"SC_workerkey={seed}; "
        f"SC_curvefamily={curve_family}; "
        f"SC_rootsmoothbits={root_smooth_bits}; "
        f"SC_candidatebits={candidate_bits}; "
        f"SC_candidatefile={json.dumps(candidate_file)}; "
        f"SC_levelfile={json.dumps(candidate_file)}; "
        f"SC_seators={sea_torsion}; "
        f"SC_branchcurves={branch_curves}; "
        "SC_progress=1; "
        f"setrand({seed}); "
        + search
        +
        'for(i=1,#c,printf("%d%s",c[i],if(i<#c," ","\\n"))); '
        'warning("curves=",SC_curves," sea_calls=",SC_seacalls,'
        '" sea_aborts=",SC_seaaborts," factor_attempts=",SC_factorattempts,'
        '" factor_aborts=",SC_factoraborts," descents=",SC_descents,'
        '" backtracks=",SC_backtracks," ecm_attempts=",SC_ecmattempts,'
        '" ecm_factors=",SC_ecmfactors," q_candidates=",SC_qcandidates,'
        '" viable_q=",SC_viableqcandidates," window_candidates=",SC_windowcandidates,'
        '" residual_bits=",SC_residualbits," max_smooth_bits=",SC_maxsmoothbits,'
        '" smooth_rejects=",SC_smoothrejects,'
        '" deep_factor_attempts=",SC_deepfactorattempts,'
        '" pm1_attempts=",SC_pm1attempts," pm1_factors=",SC_pm1factors,'
        '" pp1_attempts=",SC_pp1attempts," pp1_factors=",SC_pp1factors,'
        '" msieve_attempts=",SC_msieveattempts,'
        '" msieve_factors=",SC_msievefactors,'
        '" exhausted_orders=",SC_exhaustedorders,'
        '" factor_recoveries=",SC_factorrecoveries,'
        '" cm_tests=",SC_cmtests," cm_last_D=",SC_cmlastD,'
        '" resume_attempts=",SC_resumeattempts);\n'
    )


def smooth_least_prime(value, bound):
    """Return the least prime divisor when value is bound-smooth, else zero."""
    least = 0
    divisor = 2
    while divisor <= bound and divisor * divisor <= value:
        if value % divisor == 0:
            if not least:
                least = divisor
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        if value > bound:
            return 0
        if not least:
            least = value
    return least


def load_resume_candidates(path, p):
    files = sorted(path.glob("*.txt")) if path.is_dir() else [path]
    if not files:
        raise SystemExit(f"no candidate files found in {path}")
    candidates = {}
    exhausted_orders = set()
    exhausted_cm_screens = set()
    cm_screens_by_order = {}
    root = math.isqrt(p)
    lower_bound = root + 1 + math.isqrt(4 * root)
    trial_bound = p.bit_length() ** 2
    for candidate_file in files:
        try:
            lines = candidate_file.read_text(encoding="utf-8").splitlines()
        except OSError as error:
            raise SystemExit(f"cannot read candidate file {candidate_file}: {error}") from error
        for line_number, line in enumerate(lines, 1):
            fields = line.split()
            if len(fields) not in (5, 6, 7) or any(
                not field.isdigit() for field in fields
            ):
                raise SystemExit(
                    f"invalid candidate record at {candidate_file}:{line_number}"
                )
            if len(fields) == 5:
                cp, A, x, o, q = map(int, fields)
                if cp != p:
                    raise SystemExit(
                        f"candidate prime mismatch at {candidate_file}:{line_number}"
                    )
                if q == 1:
                    m = o
                elif o % q or q <= trial_bound or q >= root:
                    continue
                else:
                    m = o // q
                least = smooth_least_prime(m, trial_bound)
                if (
                    m <= 1
                    or not least
                    or o <= lower_bound
                    or o >= least * lower_bound
                ):
                    continue
                key = ("level", A, x, o, q)
                candidates[key] = {
                    "kind": "level",
                    "A": str(A),
                    "x": str(x),
                    "o": str(o),
                    "q": str(q),
                    "rank": 0,
                    "smooth_bits": 0,
                    "residual_bits": 0,
                    "minimum_q_bits": 0,
                    "cofactor_gap_bits": 0,
                    "source": str(candidate_file.resolve()),
                    "line": line_number,
                }
                continue
            if len(fields) == 7:
                cp, D, N, o, q, reserved1, reserved2 = map(int, fields)
                if cp != p:
                    raise SystemExit(
                        f"candidate prime mismatch at {candidate_file}:{line_number}"
                    )
                key = ("cm_screen", D, N, o, q)
                cm_order = (D, N)
                cm_screens_by_order.setdefault(cm_order, set()).add(key)
                if reserved1 == 1 and reserved2 == 0:
                    exhausted_cm_screens.add(key)
                    candidates.pop(key, None)
                    continue
                if reserved1 or reserved2:
                    raise SystemExit(
                        f"invalid screened-CM record at {candidate_file}:{line_number}"
                    )
                if key in exhausted_cm_screens:
                    continue
                # The certificate model has the rational 2-torsion point
                # (0,0), so an odd CM curve order cannot be reconstructed.
                if N % 2:
                    continue
                if q == 1:
                    m = o
                elif o % q or q <= trial_bound or q >= root:
                    continue
                else:
                    m = o // q
                least = smooth_least_prime(m, trial_bound)
                if (
                    m <= 1
                    or not least
                    or o <= lower_bound
                    or o >= least * lower_bound
                ):
                    continue
                candidates[key] = {
                    "kind": "cm_screen",
                    "D": str(D),
                    "N": str(N),
                    "o": str(o),
                    "q": str(q),
                    "rank": 1,
                    "smooth_bits": 0,
                    "residual_bits": 0,
                    "minimum_q_bits": 0,
                    "cofactor_gap_bits": 0,
                    "source": str(candidate_file.resolve()),
                    "line": line_number,
                }
                continue
            cp, A, xden, N, smooth, residual = map(int, fields)
            if cp != p:
                raise SystemExit(
                    f"candidate prime mismatch at {candidate_file}:{line_number}"
                )
            if xden == 0 and N % 2:
                continue
            key = ("order", A, xden, N)
            # A resumed order can only descend through a prime q with
            # smooth*q > lower_bound.  Once saved factoring has reduced the
            # entire composite residual below that exact threshold, no future
            # factor of the residual can work.  A reduced residual is known to
            # be composite, and all of its prime factors exceed trial_bound;
            # this gives the sharper upper bound residual/(trial_bound+1) on
            # every remaining candidate q.  Tombstone the whole order so an
            # earlier, larger residual is not replayed either: every prime
            # removed on the way to this residual was already tested.
            minimum_q = lower_bound // smooth + 1
            reduced_residual = residual != N // smooth
            if residual < minimum_q or (
                reduced_residual
                and residual // (trial_bound + 1) < minimum_q
            ):
                candidates.pop(key, None)
                exhausted_orders.add(key)
                continue
            if key in exhausted_orders:
                continue
            minimum_q_bits = minimum_q.bit_length()
            record = {
                "kind": "order",
                "A": str(A),
                "xden": str(xden),
                "N": str(N),
                "residual": str(residual),
                "smooth_bits": smooth.bit_length(),
                "residual_bits": residual.bit_length(),
                "minimum_q_bits": minimum_q_bits,
                "rank": 2,
                "cofactor_gap_bits": max(
                    0, residual.bit_length() - minimum_q_bits
                ),
                "source": str(candidate_file.resolve()),
                "line": line_number,
            }
            previous = candidates.get(key)
            if previous is None or record["residual_bits"] < previous["residual_bits"]:
                candidates[key] = record
    # Once every known viable screen for a CM order is tombstoned, repeating
    # its rough-order factoring can only reproduce the same exhausted choices.
    for (D, N), screens in cm_screens_by_order.items():
        if screens and screens <= exhausted_cm_screens:
            candidates.pop(("order", D, 0, N), None)
    if not candidates:
        raise SystemExit(f"no candidate records found in {path}")
    return sorted(
        candidates.values(),
        key=lambda candidate: (
            candidate["rank"],
            candidate["cofactor_gap_bits"],
            -candidate["smooth_bits"],
            candidate["residual_bits"],
        ),
    )


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def child_cpu_seconds():
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def append_manifest(path, record):
    if path is None:
        return
    with path.open("a", encoding="utf-8") as manifest:
        fcntl.flock(manifest, fcntl.LOCK_EX)
        json.dump(record, manifest, sort_keys=True)
        manifest.write("\n")
        manifest.flush()
        os.fsync(manifest.fileno())
        fcntl.flock(manifest, fcntl.LOCK_UN)


def read_worker_log(log):
    log.seek(0)
    return log.read().strip()


def latest_counters(detail):
    matches = COUNTERS_RE.findall(detail)
    if not matches:
        return None
    return tuple(int(counter or 0) for counter in matches[-1])


def worker_records(workers, configs):
    records = []
    for (process, log), config in zip(workers, configs):
        detail = read_worker_log(log)
        record = {**config, "returncode": process.returncode}
        counters = latest_counters(detail)
        if counters:
            for name, counter in zip(
                (
                    "curves",
                    "sea_calls",
                    "sea_aborts",
                    "factor_attempts",
                    "factor_aborts",
                    "descents",
                    "backtracks",
                    "ecm_attempts",
                    "ecm_factors",
                    "q_candidates",
                    "viable_q",
                    "window_candidates",
                    "residual_bits",
                    "max_smooth_bits",
                    "smooth_rejects",
                    "deep_factor_attempts",
                    "pm1_attempts",
                    "pm1_factors",
                    "pp1_attempts",
                    "pp1_factors",
                    "msieve_attempts",
                    "msieve_factors",
                    "exhausted_orders",
                    "factor_recoveries",
                    "cm_tests",
                    "cm_last_D",
                    "resume_attempts",
                ),
                counters,
            ):
                record[name] = counter
        resume_exhausted = RESUME_EXHAUSTED_RE.findall(detail)
        if resume_exhausted:
            record["resume_exhausted"] = int(resume_exhausted[-1])
        real_times = REAL_TIME_RE.findall(detail)
        user_times = USER_TIME_RE.findall(detail)
        sys_times = SYS_TIME_RE.findall(detail)
        if real_times:
            record["wall_seconds"] = float(real_times[-1])
        if user_times and sys_times:
            record["cpu_seconds"] = float(user_times[-1]) + float(sys_times[-1])
        records.append(record)
    return records


def stop_workers(workers):
    for process, _ in workers:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 5
    for process, _ in workers:
        if process.poll() is None:
            try:
                process.wait(timeout=max(0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
    for process, _ in workers:
        process.wait()


def verify(certificate):
    fields = certificate.split()
    if not fields or any(not field.isdigit() for field in fields):
        return False
    result = subprocess.run(
        [sys.executable, str(VERIFIER), *fields],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    return result.returncode == 0 and result.stdout.strip() == "True"


def has_seadata():
    result = subprocess.run(
        ["gp", "-q"],
        input='iferr(ellmodulareqn(211);print(1),e,print(0));\n',
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    return result.returncode == 0 and result.stdout.strip() == "1"


def is_probable_prime(p):
    result = subprocess.run(
        ["gp", "-q"],
        input=f"print(ispseudoprime({p}))\n",
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    return result.returncode == 0 and result.stdout.strip() == "1"


def main():
    args = parse_args()
    if args.cm_bound and args.cm_start > args.cm_bound:
        raise SystemExit("--cm-start must not exceed --cm-bound")
    if args.cm_only and not args.cm_bound:
        raise SystemExit("--cm-only requires --cm-bound")
    if args.cm_screen_only:
        if not args.cm_bound:
            raise SystemExit("--cm-screen-only requires --cm-bound")
        if args.resume_candidates is not None:
            raise SystemExit("--cm-screen-only cannot resume checkpoints")
        if not all(args.candidate_bits) or args.candidate_dir is None:
            raise SystemExit(
                "--cm-screen-only requires --candidate-bits and --candidate-dir"
            )
        args.cm_only = True
        args.cm_partition = "width"
        args.factor_seconds = [0]
        args.deep_factor_seconds = [0]
        args.prefactor_seconds = [0]
        args.pm1_bounds = [0]
        args.pp1_bounds = [0]
        args.ecm_bounds = [0]
        args.ecm_curves = [0]
        args.msieve_bits = [0]
        args.msieve_seconds = [0]
    if args.resume_only and args.resume_candidates is None:
        raise SystemExit("--resume-only requires --resume-candidates")
    if args.resume_only and args.cm_only:
        raise SystemExit("--resume-only and --cm-only are mutually exclusive")
    if not SHORT_GP.exists() or not VERIFIER.exists():
        raise SystemExit("short.gp and vsmallECPP.py must be beside this script")
    if not is_probable_prime(args.p):
        raise SystemExit("p must be a probable prime")
    seadata_available = has_seadata()
    if args.p.bit_length() > 50:
        if seadata_available:
            print(
                "PARI seadata detected; precomputed SEA tables are enabled",
                file=sys.stderr,
            )
        else:
            print(
                "warning: PARI seadata is not installed; large SEA calls will be "
                "significantly slower (see README.md)",
                file=sys.stderr,
            )

    run_id = str(uuid.uuid4())
    if any(args.candidate_bits):
        if args.candidate_dir is None:
            raise SystemExit("--candidate-dir is required when saving candidates")
        candidate_dir = args.candidate_dir.resolve()
        candidate_dir.mkdir(parents=True, exist_ok=True)
    else:
        candidate_dir = None

    seeds = args.seeds or list(range(args.seed, args.seed + args.workers))
    worker_count = len(seeds)
    resume_candidates = (
        load_resume_candidates(args.resume_candidates, args.p)
        if args.resume_candidates is not None
        else []
    )
    if args.resume_offset:
        if not resume_candidates:
            raise SystemExit("--resume-offset requires --resume-candidates")
        resume_candidates = resume_candidates[args.resume_offset :]
    if args.resume_top is not None:
        if not resume_candidates:
            raise SystemExit("--resume-top requires --resume-candidates")
        resume_candidates = resume_candidates[: args.resume_top]
    branch_curves = args.branch_curves_list or [args.branch_curves]
    configs = [
        {
            "worker": index + 1,
            "seed": seed,
            "factor_seconds": args.factor_seconds[index % len(args.factor_seconds)],
            "factor_flags": args.factor_flags[index % len(args.factor_flags)],
            "factor_rounds": args.factor_rounds[index % len(args.factor_rounds)],
            "deep_factor_bits": args.deep_factor_bits[
                index % len(args.deep_factor_bits)
            ],
            "deep_factor_seconds": args.deep_factor_seconds[
                index % len(args.deep_factor_seconds)
            ],
            "prefactor_seconds": args.prefactor_seconds[
                index % len(args.prefactor_seconds)
            ],
            "prefactor_flags": args.prefactor_flags[index % len(args.prefactor_flags)],
            "pm1_bound": args.pm1_bounds[index % len(args.pm1_bounds)],
            "pp1_bound": args.pp1_bounds[index % len(args.pp1_bounds)],
            "ecm_bound": args.ecm_bounds[index % len(args.ecm_bounds)],
            "ecm_curves": args.ecm_curves[index % len(args.ecm_curves)],
            "ecm_rounds": args.ecm_rounds[index % len(args.ecm_rounds)],
            "msieve_bits": args.msieve_bits[index % len(args.msieve_bits)],
            "msieve_seconds": args.msieve_seconds[
                index % len(args.msieve_seconds)
            ],
            "curve_family": args.curve_families[index % len(args.curve_families)],
            "root_smooth_bits": args.root_smooth_bits[
                index % len(args.root_smooth_bits)
            ],
            "cm_bound": args.cm_bound,
            "cm_start": args.cm_start,
            "cm_smooth_bits": args.cm_smooth_bits,
            "cm_slot": index,
            "cm_slots": worker_count,
            "cm_kind": {"both": 0, "odd": 1, "even": 2}[args.cm_kind],
            "cm_partition": {"sqrt": 0, "width": 1}[args.cm_partition],
            "cm_only": args.cm_only,
            "cm_screen_only": args.cm_screen_only,
            "resume_only": args.resume_only,
            "candidate_bits": args.candidate_bits[index % len(args.candidate_bits)],
            "candidate_file": (
                str(candidate_dir / f"{run_id}-w{index + 1}-seed{seed}.txt")
                if candidate_dir is not None
                and args.candidate_bits[index % len(args.candidate_bits)]
                else ""
            ),
            "resume_orders": (
                resume_candidates[index::worker_count][: args.resume_per_worker]
                if resume_candidates
                else []
            ),
            "sea_torsion": args.sea_torsions[index % len(args.sea_torsions)],
            "branch_curves": branch_curves[index % len(branch_curves)],
            "gp_threads": args.gp_threads,
        }
        for index, seed in enumerate(seeds)
    ]
    if any(config["curve_family"] not in (0, 1, 2, 3, 4, 5) for config in configs):
        raise SystemExit("curve families must be 0, 1, 2, 3, 4, or 5")
    if any(
        (config["ecm_bound"] and config["ecm_curves"])
        or config["pm1_bound"]
        or config["pp1_bound"]
        for config in configs
    ):
        if shutil.which("ecm") is None:
            raise SystemExit("external GMP factoring requested but ecm is absent")
    if any(config["msieve_bits"] and config["msieve_seconds"] for config in configs):
        if shutil.which("msieve") is None or shutil.which("timeout") is None:
            raise SystemExit("msieve strategy requested but msieve or timeout is absent")
    workers = []
    started = time.monotonic()
    started_at = utc_now()
    started_cpu = child_cpu_seconds()
    append_manifest(
        args.manifest,
        {
            "event": "start",
            "run_id": run_id,
            "started_at": started_at,
            "prime": str(args.p),
            "parallel_attempts": worker_count,
            "seadata_available": seadata_available,
            "workers": configs,
        },
    )
    try:
        for config in configs:
            log = tempfile.TemporaryFile(mode="w+")
            process = subprocess.Popen(
                ["/usr/bin/time", "-p", "gp", "-q", str(SHORT_GP)],
                cwd=ROOT,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=log,
                text=True,
                start_new_session=True,
            )
            process.stdin.write(
                worker_input(
                    args.p,
                    config["factor_seconds"],
                    config["factor_flags"],
                    config["factor_rounds"],
                    config["deep_factor_bits"],
                    config["deep_factor_seconds"],
                    config["prefactor_seconds"],
                    config["prefactor_flags"],
                    config["pm1_bound"],
                    config["pp1_bound"],
                    config["ecm_bound"],
                    config["ecm_curves"],
                    config["ecm_rounds"],
                    config["msieve_bits"],
                    config["msieve_seconds"],
                    config["curve_family"],
                    config["root_smooth_bits"],
                    config["candidate_bits"],
                    config["candidate_file"],
                    config["sea_torsion"],
                    config["branch_curves"],
                    config["seed"],
                    config["gp_threads"],
                    config["resume_orders"],
                    config["cm_start"],
                    config["cm_bound"],
                    config["cm_smooth_bits"],
                    config["cm_slot"],
                    config["cm_slots"],
                    config["cm_kind"],
                    config["cm_partition"],
                    config["cm_only"],
                    config["resume_only"],
                )
            )
            process.stdin.close()
            workers.append((process, log))
        print(
            f"started {worker_count} workers for p={args.p} "
            f"with seeds {seeds} "
            f"and factor limits "
            f"{[config['factor_seconds'] for config in configs]}; "
            f"factor flags "
            f"{[config['factor_flags'] for config in configs]}; "
            f"factor rounds "
            f"{[config['factor_rounds'] for config in configs]}; "
            f"deep-factor thresholds "
            f"{[config['deep_factor_bits'] for config in configs]}; "
            f"deep-factor limits "
            f"{[config['deep_factor_seconds'] for config in configs]}; "
            f"saved prefactor limits "
            f"{[config['prefactor_seconds'] for config in configs]}; "
            f"prefactor flags {[config['prefactor_flags'] for config in configs]}; "
            f"P-1 bounds {[config['pm1_bound'] for config in configs]}; "
            f"P+1 bounds {[config['pp1_bound'] for config in configs]}; "
            f"GMP-ECM bounds {[config['ecm_bound'] for config in configs]}; "
            f"ECM curves {[config['ecm_curves'] for config in configs]}; "
            f"ECM rounds {[config['ecm_rounds'] for config in configs]}; "
            f"msieve thresholds {[config['msieve_bits'] for config in configs]}; "
            f"msieve limits {[config['msieve_seconds'] for config in configs]}; "
            f"curve families {[config['curve_family'] for config in configs]}; "
            f"root smooth-bit thresholds "
            f"{[config['root_smooth_bits'] for config in configs]}; "
            f"CM discriminant range {args.cm_start}..{args.cm_bound}; "
            f"CM smooth-bit threshold {args.cm_smooth_bits}; "
            f"CM discriminant types {args.cm_kind}; "
            f"CM partition {args.cm_partition}; "
            f"CM-only mode {args.cm_only}; "
            f"CM screen-only mode {args.cm_screen_only}; "
            f"candidate save thresholds "
            f"{[config['candidate_bits'] for config in configs]}; "
            f"resume checkpoints "
            f"{len(resume_candidates)} ranked checkpoint(s); "
            f"resume queues {[len(config['resume_orders']) for config in configs]}; "
            f"resume-only mode {args.resume_only}; "
            f"SEA torsion filters {[config['sea_torsion'] for config in configs]}; "
            f"child-branch budgets {[config['branch_curves'] for config in configs]}; "
            f"{args.gp_threads} PARI thread(s) each",
            file=sys.stderr,
            flush=True,
        )

        last_report = started
        exhausted_workers = set()
        while True:
            for index, (process, log) in enumerate(workers):
                if index in exhausted_workers:
                    continue
                returncode = process.poll()
                if returncode is None:
                    continue
                output = process.stdout.read().strip()
                if returncode == 0 and verify(output):
                    elapsed = time.monotonic() - started
                    worker_detail = read_worker_log(log).splitlines()
                    stop_workers(workers)
                    cpu_seconds = child_cpu_seconds() - started_cpu
                    records = worker_records(workers, configs)
                    winning_cpu = records[index].get("cpu_seconds")
                    counter_detail = (
                        next(
                            (
                                line
                                for line in reversed(worker_detail)
                                if "curves=" in line
                            ),
                            worker_detail[-1],
                        )
                        if worker_detail
                        else ""
                    )
                    if args.output:
                        args.output.write_text(output + "\n")
                    print(output)
                    print(
                        f"worker {index + 1}/{worker_count} found and verified a certificate "
                        f"in {elapsed:.1f} seconds after launching "
                        f"{worker_count} parallel worker attempt(s); "
                        f"{cpu_seconds:.1f} aggregate child CPU seconds"
                        + (
                            f"; {winning_cpu:.1f} winning-worker CPU seconds"
                            if winning_cpu is not None
                            else ""
                        )
                        + (f" ({counter_detail})" if counter_detail else ""),
                        file=sys.stderr,
                    )
                    append_manifest(
                        args.manifest,
                        {
                            "event": "finish",
                            "run_id": run_id,
                            "finished_at": utc_now(),
                            "outcome": "success",
                            "parallel_attempts": worker_count,
                            "wall_seconds": round(elapsed, 3),
                            "child_cpu_seconds": round(cpu_seconds, 3),
                            "winning_worker": index + 1,
                            "winning_seed": configs[index]["seed"],
                            "workers": records,
                        },
                    )
                    return 0
                if (args.cm_only or args.resume_only) and returncode == 0 and not output:
                    exhausted_workers.add(index)
                    if len(exhausted_workers) == worker_count:
                        elapsed = time.monotonic() - started
                        cpu_seconds = child_cpu_seconds() - started_cpu
                        records = worker_records(workers, configs)
                        print(
                            f"completed the finite search without a certificate in "
                            f"{elapsed:.1f} seconds across {worker_count} parallel "
                            f"worker attempt(s); {cpu_seconds:.1f} aggregate child "
                            f"CPU seconds",
                            file=sys.stderr,
                        )
                        append_manifest(
                            args.manifest,
                            {
                                "event": "finish",
                                "run_id": run_id,
                                "finished_at": utc_now(),
                                "outcome": "exhausted",
                                "parallel_attempts": worker_count,
                                "wall_seconds": round(elapsed, 3),
                                "child_cpu_seconds": round(cpu_seconds, 3),
                                "workers": records,
                            },
                        )
                        return 0
                    continue
                detail = read_worker_log(log)
                raise RuntimeError(
                    f"worker {index} exited without a valid certificate "
                    f"(status {returncode}): {detail or output or 'no output'}"
                )

            now = time.monotonic()
            if now - last_report >= 60:
                progress = []
                for index, (_, log) in enumerate(workers):
                    counters = latest_counters(read_worker_log(log))
                    if counters:
                        (
                            curves,
                            sea_calls,
                            sea_aborts,
                            factor_attempts,
                            factor_aborts,
                            descents,
                            backtracks,
                            ecm_attempts,
                            ecm_factors,
                            q_candidates,
                            viable_q,
                            window_candidates,
                            residual_bits,
                            max_smooth_bits,
                            smooth_rejects,
                            deep_factor_attempts,
                            pm1_attempts,
                            pm1_factors,
                            pp1_attempts,
                            pp1_factors,
                            msieve_attempts,
                            msieve_factors,
                            exhausted_orders,
                            factor_recoveries,
                            cm_tests,
                            cm_last_D,
                            resume_attempts,
                        ) = counters
                        progress.append(
                            f"w{index + 1}"
                            f"{'/done' if index in exhausted_workers else ''}:"
                            f"{curves}c/{sea_calls}s/{sea_aborts}a "
                            f"{factor_attempts}f/{factor_aborts}x "
                            f"{descents}d/{backtracks}b "
                            f"{ecm_attempts}e/{ecm_factors}g "
                            f"{q_candidates}q/{viable_q}v/{window_candidates}w "
                            f"{residual_bits}r/{max_smooth_bits}m/{smooth_rejects}j "
                            f"{deep_factor_attempts}h "
                            f"{pm1_attempts}p-/{pm1_factors}g- "
                            f"{pp1_attempts}p+/{pp1_factors}g+ "
                            f"{msieve_attempts}ms/{msieve_factors}mg "
                            f"{exhausted_orders}z/{factor_recoveries}fr"
                            f" {cm_tests}ct/{cm_last_D}cd"
                            f" {resume_attempts}ra"
                        )
                    elif index in exhausted_workers:
                        progress.append(f"w{index + 1}:done")
                print(
                    f"searching: {now - started:.0f} seconds elapsed"
                    + (f"; {' '.join(progress)}" if progress else ""),
                    file=sys.stderr,
                )
                last_report = now
            time.sleep(0.25)
    except (KeyboardInterrupt, Exception) as error:
        stop_workers(workers)
        elapsed = time.monotonic() - started
        cpu_seconds = child_cpu_seconds() - started_cpu
        append_manifest(
            args.manifest,
            {
                "event": "finish",
                "run_id": run_id,
                "finished_at": utc_now(),
                "outcome": "interrupted" if isinstance(error, KeyboardInterrupt) else "error",
                "parallel_attempts": worker_count,
                "wall_seconds": round(elapsed, 3),
                "child_cpu_seconds": round(cpu_seconds, 3),
                "workers": worker_records(workers, configs),
                **({"error": str(error)} if not isinstance(error, KeyboardInterrupt) else {}),
            },
        )
        if isinstance(error, KeyboardInterrupt):
            print("search interrupted", file=sys.stderr)
            return 130
        raise
    finally:
        for _, log in workers:
            log.close()


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    raise SystemExit(main())
