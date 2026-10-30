#!/usr/bin/env python3
"""bppancestor is probed only when ancestors are asked for, and then once.

Usage:
    python3 tests/test_bppancestor_probe_only_when_needed.py PAGAN2

BppAncestors::test_executable() asks whether bppancestor exists by RUNNING
it through system(3): a shell and bppancestor itself, twice when it is not
beside the pagan2 binary. Writing an alignment called it before looking at
whether ancestors were wanted, so a placement run that asks for none paid it
anyway -- once per alignment written, so once per query with
--independent-queries, where on a host that reads programs from slow shared
storage it was most of each query's time. --no-bppancestors could not stop it
either. Now it is asked only when ancestors are wanted, once per process,
before any work, and processes forked later inherit the answer.

THE CHECK. A stand-in `bppancestor`, first on PATH, records every start and
exits 1 ("not usable"), so pagan2 falls back to its own ancestors and no run
depends on the real program. Counted:
  1. a placement run asking for no ancestors: 0 starts;
  2. the same with --independent-queries, when the binary has that option:
     0 starts for three queries;
  3. with --output-ancestors and --no-bppancestors: 0 starts;
  4. with --output-ancestors: exactly 1 start -- the probe still runs when it
     can matter (the PREMISE: the stand-in is what it finds), and once;
  5. --independent-queries --output-ancestors, when the binary has that
     option: exactly 1 start for three queries, not one per query.

Exit 0 and a PASS line, or exit 1 naming what failed.
"""
import os
import random
import shutil
import stat
import subprocess
import sys
import tempfile

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS]
TREE = "((A:0.05,B:0.05):0.02,(C:0.05,D:0.05):0.02);\n"


def _mutate(rng, cod, n):
    cod = list(cod)
    for i in rng.sample(range(1, len(cod)), n):
        cod[i] = rng.choice([c for c in SENSE if c != cod[i]])
    return cod


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(">%s\n%s\n" % (name, seq))


def make_data(workdir):
    rng = random.Random(20261030)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(99)]
    taxa = {"A": _mutate(rng, anc, 6), "B": _mutate(rng, anc, 6),
            "C": _mutate(rng, anc, 10), "D": _mutate(rng, anc, 10)}
    _write_fasta(os.path.join(workdir, "ref.aln.fas"),
                 [(t, "".join(c)) for t, c in sorted(taxa.items())])
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)
    _write_fasta(os.path.join(workdir, "queries.fas"),
                 [("q1", "".join(_mutate(rng, taxa["A"], 3))),
                  ("q2", "".join(_mutate(rng, taxa["C"], 3))),
                  ("q3", "".join(_mutate(rng, taxa["D"], 2)))])


def make_stand_in(bindir, log):
    """A `bppancestor` that records its start and reports itself unusable."""
    os.makedirs(bindir)
    path = os.path.join(bindir, "bppancestor")
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\necho started >> '%s'\nexit 1\n" % log)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def starts(log):
    if not os.path.exists(log):
        return 0
    with open(log) as fh:
        return sum(1 for _ in fh)


def run(binary, workdir, label, extra, env):
    out_dir = os.path.join(workdir, label)
    os.makedirs(os.path.join(out_dir, "tmp"))
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(out_dir, "out"),
           "--temp-folder", os.path.join(out_dir, "tmp")] + list(extra)
    with open(os.path.join(out_dir, "log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=out_dir,
                            env=env).returncode
    return rc, out_dir


def knows(binary, option):
    out = subprocess.run([binary, "--help"], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                         text=True, errors="replace").stdout
    return option in out


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    workdir = tempfile.mkdtemp(prefix="pagan2_bppancestor_probe.")
    make_data(workdir)
    failures = []
    cases = [("plain", [], 0),
             ("no_bppancestors", ["--output-ancestors", "--no-bppancestors"], 0),
             ("ancestors", ["--output-ancestors"], 1)]
    if knows(binary, "--independent-queries"):
        cases.insert(1, ("independent", ["--independent-queries"], 0))
        cases.append(("independent_ancestors",
                      ["--independent-queries", "--output-ancestors"], 1))
    seen = {}
    for label, extra, want in cases:
        log = os.path.join(workdir, label + ".starts")
        bindir = os.path.join(workdir, label + ".bin")
        make_stand_in(bindir, log)
        env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"])
        rc, out_dir = run(binary, workdir, label, extra, env)
        n = seen[label] = starts(log)
        if rc != 0:
            failures.append("%s: exit %d (see %s/log)" % (label, rc, out_dir))
        elif want is not None and n != want:
            failures.append("%s: bppancestor started %d time(s), expected %d"
                            % (label, n, want))
    if seen.get("ancestors", 0) < 1:
        failures.append("PREMISE: with --output-ancestors the stand-in was never "
                        "started, so a count of 0 elsewhere proves nothing")
    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and logs kept in %s)" % workdir)
        return 1
    print("PASS: bppancestor starts %s" % ", ".join(
        "%s=%d" % kv for kv in sorted(seen.items())))
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
