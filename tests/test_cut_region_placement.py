#!/usr/bin/env python3
"""Synthetic placement test: a region cut out of the query, and a region cut
out of the reference, placed with --guidetree --codons.

Usage:
    python3 tests/test_cut_region_placement.py PAGAN2 [BASELINE_PAGAN2]

Builds its own data (deterministic, standard library only) in a temporary
directory -- TMPDIR is honoured -- and runs PAGAN2 the way a codon-aware
read-placement pipeline does:

    pagan2 -a ref.aln.fas -t ref.tre -q queries.fas -o out --guidetree --codons

THE DATA
  * Four related codon sequences (ATG + 189 sense codons, no stop) on a
    four-taxon tree. Codons 120..129 (30 nt) are CUT OUT of all four, so a
    query that keeps them carries an INSERTION relative to the reference.
  * Queries, all from taxon A:
      cut_query   codons  60..69 cut out           -> a 30 nt DELETION
      cut_ref     keeps   120..129                 -> a 30 nt INSERTION
      cut_both    both of the above
      plain_1..3  A with a few substitutions, nothing cut
    so several reads are placed in ONE run and anything built once and
    reused across reads (the read-scoring model) is reused here.
  The codons flanking each cut are chosen so the gap cannot slide: the first
  codon of the cut differs from the codon after it, and the last from the
  codon before it.

WHAT IS CHECKED
  pagan2 writes a query once per tree node it is placed at (``q``, ``q.1``,
  ``q.2``); every copy is checked. Indels are judged PAIRWISE -- a query row
  against each reference row, dropping the columns where both have a gap --
  because a query is placed independently: another query's insertion adds
  columns that are gaps in both, and they say nothing about this one.
  1. No base is lost or invented: every query row (copies included) and
     every reference row, with its gaps removed, is exactly its input.
  2. cut_query / cut_both: against every reference, the deletion is ONE
     contiguous run of 30 reference bases facing query gaps, and it starts
     after exactly the 180 query bases that precede it.
  3. cut_ref / cut_both: against every reference, the insertion is ONE
     contiguous run of 30 query bases facing reference gaps, and they are
     the 30 bases that were cut out of the references.
  4. plain_*: no indel at all against any reference.
  5. Two runs of the same binary give the same rows (compared by name: the
     ORDER pagan2 writes the placed copies in is not stable run to run, in
     either build, and is reported but not failed).
  6. With BASELINE_PAGAN2: the same rows as PAGAN2's.

Exit 0 and a PASS line, or exit 1 naming each failed check.
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

CODON_TABLE_STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in CODON_TABLE_STOPS]

N_CODONS = 190                 # ATG + 189 sense codons
DEL = (60, 70)                 # codons cut out of the query
INS = (120, 130)               # codons cut out of every reference
TREE = "((A:0.05,B:0.05):0.02,(C:0.05,D:0.05):0.02);\n"


def _codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def _ancestor(rng):
    """ATG + sense codons, with codons that make both cuts unslidable."""
    while True:
        cod = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
        ok = True
        for lo, hi in (DEL, INS):
            if cod[lo] == cod[hi] or cod[hi - 1] == cod[lo - 1]:
                ok = False
        if ok:
            return cod


def _mutate(rng, cod, n):
    """*n* substitutions to other SENSE codons, never at a cut boundary."""
    cod = list(cod)
    protected = {DEL[0] - 1, DEL[0], DEL[1] - 1, DEL[1],
                 INS[0] - 1, INS[0], INS[1] - 1, INS[1], 0}
    sites = [i for i in range(len(cod)) if i not in protected]
    for i in rng.sample(sites, n):
        cod[i] = rng.choice([c for c in SENSE if c != cod[i]])
    return cod


def _cut(cod, *regions):
    keep = [c for i, c in enumerate(cod)
            if not any(lo <= i < hi for lo, hi in regions)]
    return "".join(keep)


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(">%s\n%s\n" % (name, seq))


def _read_fasta(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name is not None:
                out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def _pairwise_runs(qrow, rrow):
    """Indel runs of *qrow* against *rrow*, ignoring both-gap columns.

    ``(deleted, inserted)``: lists of ``(query_bases_before, length)`` for
    each INTERNAL run of reference bases facing query gaps (a deletion) and
    of query bases facing reference gaps (an insertion).
    """
    pairs = [(q, r) for q, r in zip(qrow, rrow) if q != "-" or r != "-"]
    deleted, inserted = [], []
    q_seen, i = 0, 0
    while i < len(pairs):
        q, r = pairs[i]
        if q == "-" or r == "-":
            kind = "del" if q == "-" else "ins"
            j, at = i, q_seen
            while j < len(pairs) and (
                    (kind == "del" and pairs[j][0] == "-")
                    or (kind == "ins" and pairs[j][1] == "-")):
                if pairs[j][0] != "-":
                    q_seen += 1
                j += 1
            internal = i > 0 and j < len(pairs)
            if internal:
                (deleted if kind == "del" else inserted).append((at, j - i))
            i = j
        else:
            q_seen += 1
            i += 1
    return deleted, inserted


def make_data(workdir):
    rng = random.Random(20261003)
    anc = _ancestor(rng)
    taxa = {"A": _mutate(rng, anc, 6), "B": _mutate(rng, anc, 6),
            "C": _mutate(rng, anc, 9), "D": _mutate(rng, anc, 9)}
    refs = [(t, _cut(c, INS)) for t, c in sorted(taxa.items())]
    a = taxa["A"]
    queries = [
        ("cut_query", _cut(a, DEL, INS)),
        ("cut_ref", _cut(a)),
        ("cut_both", _cut(a, DEL)),
    ]
    for k in range(1, 4):
        queries.append(("plain_%d" % k, _cut(_mutate(rng, a, 3), INS)))
    _write_fasta(os.path.join(workdir, "ref.aln.fas"), refs)
    _write_fasta(os.path.join(workdir, "queries.fas"), queries)
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)
    inserted = "".join(a[INS[0]:INS[1]])
    return dict(refs), dict(queries), inserted


def run_pagan2(binary, workdir, label):
    out_dir = os.path.join(workdir, label)
    os.makedirs(os.path.join(out_dir, "tmp"))
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(out_dir, "out"), "--guidetree", "--codons",
           "--temp-folder", os.path.join(out_dir, "tmp")]
    with open(os.path.join(out_dir, "log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=out_dir).returncode
    out = os.path.join(out_dir, "out.fas")
    if rc != 0 or not os.path.isfile(out):
        raise SystemExit("FAIL: %s exited %d, no %s (see %s)"
                         % (binary, rc, out, os.path.join(out_dir, "log")))
    with open(out, "rb") as fh:
        return fh.read(), _read_fasta(out)


def _base_name(name):
    stem, dot, tail = name.rpartition(".")
    return stem if dot and tail.isdigit() else name


def check(refs, queries, inserted, rows):
    failures = []
    for name in queries:
        if name not in rows:
            failures.append("%s: missing from the output" % name)
    for name, row in rows.items():
        seq = refs.get(name) or queries.get(_base_name(name))
        if seq is None:
            failures.append("%s: an output row nobody put in" % name)
        elif row.replace("-", "") != seq:
            failures.append("%s: bases differ from the input (lost or "
                            "invented)" % name)
    cut_len = 3 * (DEL[1] - DEL[0])
    for name, row in rows.items():
        base = _base_name(name)
        if base not in queries:
            continue
        want_del = [(3 * DEL[0], cut_len)] if base in (
            "cut_query", "cut_both") else []
        want_ins = [(3 * INS[0] - (cut_len if base == "cut_both" else 0),
                     cut_len)] if base in ("cut_ref", "cut_both") else []
        for ref in refs:
            deleted, inserted_runs = _pairwise_runs(row, rows[ref])
            if deleted != want_del:
                failures.append("%s vs %s: deletions %s, expected %s"
                                % (name, ref, deleted, want_del))
            if inserted_runs != want_ins:
                failures.append("%s vs %s: insertions %s, expected %s"
                                % (name, ref, inserted_runs, want_ins))
            elif want_ins:
                at = want_ins[0][0]
                if row.replace("-", "")[at:at + cut_len] != inserted:
                    failures.append("%s: the inserted bases are not the "
                                    "region cut from the references" % name)
    return failures


def _order(raw):
    return [l[1:].strip() for l in raw.decode().splitlines()
            if l.startswith(">")]


def main(argv):
    if len(argv) not in (2, 3):
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    baseline = os.path.abspath(argv[2]) if len(argv) == 3 else None
    workdir = tempfile.mkdtemp(prefix="pagan2_cut_region.")
    keep = False
    try:
        refs, queries, inserted = make_data(workdir)
        raw, rows = run_pagan2(binary, workdir, "run1")
        failures = check(refs, queries, inserted, rows)
        notes = []
        raw2, rows2 = run_pagan2(binary, workdir, "run2")
        if rows2 != rows:
            failures.append("two runs of %s give different rows" % binary)
        elif _order(raw2) != _order(raw):
            notes.append("two runs wrote the same rows in a different order")
        if baseline:
            raw_b, rows_b = run_pagan2(baseline, workdir, "baseline")
            if rows_b != rows:
                failures.append("the rows differ from the baseline %s"
                                % baseline)
            elif _order(raw_b) != _order(raw):
                notes.append("the baseline wrote the same rows in a "
                             "different order")
        for note in notes:
            print("note:", note)
        if failures:
            keep = True
            for f in failures:
                print("FAIL:", f)
            print("inputs and outputs kept in", workdir)
            return 1
        print("PASS: %d queries and %d references placed; no base lost or "
              "invented; the cut regions are contiguous single gap runs%s"
              % (len(queries), len(refs),
                 "; the same rows as the baseline" if baseline else ""))
        return 0
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
