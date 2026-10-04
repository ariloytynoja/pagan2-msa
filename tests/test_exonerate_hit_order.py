#!/usr/bin/env python3
"""Read placement must not depend on the ORDER exonerate prints its hits in.

Usage:
    python3 tests/test_exonerate_hit_order.py PAGAN2

exonerate does not print the hits of a run in a reproducible order: on one
unchanged query and target file, exonerate 2.4.0 printed the same set of
sugar lines in four different orders in twelve runs. pagan2 read them
order-sensitively -- the hits on one target were folded in print order, so
which strand won depended on it; the targets were ranked in first-seen order
and std::sort is not stable -- so the same batch of reads could be placed,
and aligned, differently from one run to the next.

On the data below, v1.6 writes the same rows in a different ORDER for the
two hit orders, which is what this test catches. The order can also move
bases: with two queries carrying 30 nt insertions at the same site, v1.6
put one query's bases before or after the other's insertion columns
depending on the hit order. src/test_exonerate_fold_hits.cpp pins the fold
itself.

THE DATA (deterministic, standard library only, written under TMPDIR)
  Four codon sequences (ATG + 149 sense codons) on a four-taxon tree, where
  A and B are IDENTICAL, so every read derived from A scores the same against
  A as against B and the tie is decided by order alone. Six reads: A with a
  few substitutions each.

WHAT IS CHECKED
  PAGAN2 is run twice with `exonerate` on PATH replaced by a shim that runs
  the real exonerate and prints its sugar lines sorted ascending in one run
  and descending in the other -- the same hits, in two orders. The two
  out.fas files must be identical, and the shim must have been called in
  both runs (a pass with no exonerate call would prove nothing).

Exit 0 and a PASS line, or exit 1 naming what differed.
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
N_CODONS = 150
TREE = "((A:0.05,B:0.05):0.02,(C:0.05,D:0.05):0.02);\n"

SHIM = r"""#!/bin/bash
# Run the real exonerate, then print its output with the sugar lines in the
# order EXO_ORDER names. Its exit status is exonerate's.
out=$("$REAL_EXONERATE" "$@")
rc=$?
echo x >> "$EXO_CALLS"
printf '%s\n' "$out" | grep -v '^sugar:'
if [ "$EXO_ORDER" = reverse ]; then
    printf '%s\n' "$out" | grep '^sugar:' | LC_ALL=C sort -r
else
    printf '%s\n' "$out" | grep '^sugar:' | LC_ALL=C sort
fi
exit $rc
"""


def _mutate(rng, cod, n):
    cod = list(cod)
    for i in rng.sample(range(1, len(cod)), n):
        cod[i] = rng.choice([c for c in SENSE if c != cod[i]])
    return cod


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(">%s\n%s\n" % (name, seq))


def _read_fasta(path):
    rows, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                rows[name] = ""
            elif name is not None:
                rows[name] += line
    return rows


def make_data(workdir):
    rng = random.Random(20261004)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
    a = _mutate(rng, anc, 6)
    taxa = {"A": a, "B": list(a),
            "C": _mutate(rng, anc, 12), "D": _mutate(rng, anc, 12)}
    _write_fasta(os.path.join(workdir, "ref.aln.fas"),
                 [(t, "".join(c)) for t, c in sorted(taxa.items())])
    _write_fasta(os.path.join(workdir, "queries.fas"),
                 [("read_%d" % k, "".join(_mutate(rng, a, 3)))
                  for k in range(1, 7)])
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)


def run_pagan2(binary, workdir, shim_dir, order):
    out_dir = os.path.join(workdir, order)
    os.makedirs(os.path.join(out_dir, "tmp"))
    calls = os.path.join(out_dir, "exonerate_calls")
    env = dict(os.environ, EXO_ORDER=order, EXO_CALLS=calls,
               PATH=shim_dir + os.pathsep + os.environ.get("PATH", ""))
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(out_dir, "out"), "--guidetree", "--codons",
           "--temp-folder", os.path.join(out_dir, "tmp")]
    with open(os.path.join(out_dir, "log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=out_dir,
                            env=env).returncode
    out = os.path.join(out_dir, "out.fas")
    if rc != 0 or not os.path.isfile(out):
        raise SystemExit("FAIL: %s exited %d, no %s (see %s)"
                         % (binary, rc, out, os.path.join(out_dir, "log")))
    n_calls = 0
    if os.path.isfile(calls):
        with open(calls) as fh:
            n_calls = sum(1 for _ in fh)
    with open(out, "rb") as fh:
        return fh.read(), _read_fasta(out), n_calls


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    real = shutil.which("exonerate")
    if real is None:
        print("SKIP: no exonerate on PATH")
        return 0
    workdir = tempfile.mkdtemp(prefix="pagan2_exonerate_order.")
    shim_dir = os.path.join(workdir, "shim")
    os.makedirs(shim_dir)
    shim = os.path.join(shim_dir, "exonerate")
    with open(shim, "w") as fh:
        fh.write(SHIM)
    os.chmod(shim, stat.S_IRWXU)
    os.environ["REAL_EXONERATE"] = real
    make_data(workdir)

    fwd_bytes, fwd_rows, fwd_calls = run_pagan2(binary, workdir, shim_dir,
                                                "forward")
    rev_bytes, rev_rows, rev_calls = run_pagan2(binary, workdir, shim_dir,
                                                "reverse")
    failures = []
    if not fwd_calls or not rev_calls:
        failures.append("the exonerate shim was not called (forward %d, "
                        "reverse %d): nothing was tested"
                        % (fwd_calls, rev_calls))
    if fwd_rows != rev_rows:
        moved = sorted(set(fwd_rows) ^ set(rev_rows))
        changed = sorted(n for n in set(fwd_rows) & set(rev_rows)
                         if fwd_rows[n] != rev_rows[n])
        failures.append("the rows differ with exonerate's hits in the other "
                        "order: names only in one run %s, rows changed %s"
                        % (moved or "none", changed or "none"))
    elif fwd_bytes != rev_bytes:
        failures.append("the same rows were written in a different order")
    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and logs kept in %s)" % workdir)
        return 1
    print("PASS: %d rows identical with exonerate's hits in either order "
          "(%d + %d exonerate calls)" % (len(fwd_rows), fwd_calls, rev_calls))
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
