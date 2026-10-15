#!/usr/bin/env python3
"""--independent-queries aligns each query exactly as a run given only it.

Usage:
    python3 tests/test_independent_queries.py PAGAN2

By default the queries of one --queryfile are placed into ONE growing tree
and alignment: every placed query, and the internal node made for it, is
part of what the next query is placed on and aligned against. So a query's
row depends on the other queries of the file. With --independent-queries the
reference is set up once and each query is aligned alone, in a child process
forked from that state, to "<outfile>.<n>"; "<outfile>.queries.tsv" lists
every query's outcome.

THE DATA (deterministic, standard library only, written under TMPDIR)
  Four codon sequences (ATG + 149 sense codons) on a four-taxon tree, and
  six queries: four carry a few substitutions against one taxon each, one
  of them also a two-codon insertion -- which adds alignment columns, so in
  a default run it changes the other queries' rows too -- and one is
  unrelated sequence that places nowhere.

WHAT IS CHECKED
  1. one --independent-queries run of all six: the manifest lists the six in
     file order, five tab-separated fields a row, and ends with the line
     "# complete: 6 queries, 5 placed, 1 unplaced, 0 failed" (no tab in it);
     each placed query's out.<n>.fas and out.<n>.nhx_tree are
     byte-identical to those of a run of that query alone, with the same
     options; the unrelated query is "unplaced" and, as alone, has no file;
  2. a second --independent-queries run writes byte-identical files;
  3. PREMISE: in a default run of all six, at least one query's row differs
     from its row alone (otherwise the data shows nothing about batches);
  4. --independent-queries without a reference alignment is refused, with a
     non-zero exit.

Exit 0 and a PASS line, or exit 1 naming what failed.
"""
import filecmp
import os
import random
import shutil
import subprocess
import sys
import tempfile

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS]
N_CODONS = 150
TREE = "((A:0.05,B:0.05):0.02,(C:0.05,D:0.05):0.02);\n"
OPTIONS = ["--guidetree", "--codons", "--one-placement-only"]


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
    rng = random.Random(20261014)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
    taxa = {"A": _mutate(rng, anc, 8), "B": _mutate(rng, anc, 8),
            "C": _mutate(rng, anc, 14), "D": _mutate(rng, anc, 14)}
    _write_fasta(os.path.join(workdir, "ref.aln.fas"),
                 [(t, "".join(c)) for t, c in sorted(taxa.items())])
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)
    q1 = _mutate(rng, taxa["A"], 3)
    q1 = q1[:70] + [rng.choice(SENSE), rng.choice(SENSE)] + q1[70:]
    queries = [("query_1", q1),
               ("query_2", _mutate(rng, taxa["A"], 4)),
               ("query_3", _mutate(rng, taxa["C"], 3)),
               ("unrelated", [rng.choice(SENSE) for _ in range(N_CODONS)]),
               ("query_5", _mutate(rng, taxa["B"], 5)),
               ("query_6", _mutate(rng, taxa["D"], 2))]
    queries = [(n, "".join(c)) for n, c in queries]
    _write_fasta(os.path.join(workdir, "queries.fas"), queries)
    for k, (n, s) in enumerate(queries, 1):
        _write_fasta(os.path.join(workdir, "query%d.fas" % k), [(n, s)])
    return [n for n, _ in queries]


def read_manifest(path):
    """(rows, trailer, problems): the rows of five fields, the final '#' line
    or None, and what is wrong with the file's shape."""
    problems = []
    with open(path) as fh:
        lines = fh.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if lines[:1] != ["query\tname\tstatus\tdetail\talignment"]:
        problems.append("manifest header is %r" % (lines[:1],))
    rows, trailer = [], None
    for k, line in enumerate(lines[1:], 2):
        if line.startswith("#"):
            if k != len(lines):
                problems.append("manifest line %d starts with '#' but is not the "
                                "last line" % k)
            trailer = line
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            problems.append("manifest line %d has %d fields, not 5: %r"
                            % (k, len(fields), line))
            continue
        rows.append(fields)
    return rows, trailer, problems


def run_pagan2(binary, workdir, label, queryfile, extra=(), reference=True):
    out_dir = os.path.join(workdir, label)
    os.makedirs(os.path.join(out_dir, "tmp"))
    cmd = [binary]
    if reference:
        cmd += ["-a", os.path.join(workdir, "ref.aln.fas"),
                "-t", os.path.join(workdir, "ref.tre")]
    cmd += ["-q", os.path.join(workdir, queryfile),
            "-o", os.path.join(out_dir, "out"),
            "--temp-folder", os.path.join(out_dir, "tmp")] + OPTIONS + list(extra)
    with open(os.path.join(out_dir, "log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=out_dir).returncode
    return rc, out_dir


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    workdir = tempfile.mkdtemp(prefix="pagan2_independent_queries.")
    names = make_data(workdir)
    failures = []

    alone = {}
    for k, name in enumerate(names, 1):
        rc, d = run_pagan2(binary, workdir, "alone%d" % k, "query%d.fas" % k)
        if rc != 0:
            failures.append("%s alone: exit %d (see %s/log)" % (name, rc, d))
        alone[k] = d

    runs = []
    for label in ("independent", "independent_again"):
        rc, d = run_pagan2(binary, workdir, label, "queries.fas",
                           ["--independent-queries"])
        if rc != 0:
            failures.append("%s: exit %d (see %s/log)" % (label, rc, d))
        runs.append(d)

    indep = runs[0]
    manifest = os.path.join(indep, "out.queries.tsv")
    if not os.path.isfile(manifest):
        failures.append("no %s" % manifest)
        rows = []
    else:
        rows, trailer, problems = read_manifest(manifest)
        failures.extend(problems)
        if [(r[0], r[1]) for r in rows] != [(str(k), n) for k, n in enumerate(names, 1)]:
            failures.append("manifest lists %r, expected the queries in file order"
                            % [(r[0], r[1]) for r in rows])
        expected = "# complete: 6 queries, 5 placed, 1 unplaced, 0 failed"
        if trailer != expected:
            failures.append("manifest ends with %r, expected %r" % (trailer, expected))

    n_identical = 0
    for k, name in enumerate(names, 1):
        row = rows[k - 1] if k <= len(rows) else None
        if row is None:
            continue
        a_fas = os.path.join(alone[k], "out.fas")
        i_fas = os.path.join(indep, "out.%d.fas" % k)
        if name == "unrelated":
            if row[2] != "unplaced" or row[4] != "":
                failures.append("unrelated: manifest says %r, expected unplaced, no file" % row)
            if os.path.exists(a_fas) or os.path.exists(i_fas):
                failures.append("unrelated: an alignment was written (alone %s, independent %s)"
                                % (os.path.exists(a_fas), os.path.exists(i_fas)))
            continue
        if row[2] != "placed" or row[4] != i_fas:
            failures.append("%s: manifest says %r, expected placed, %s" % (name, row, i_fas))
            continue
        if not os.path.isfile(a_fas):
            failures.append("%s: placed in the independent run but not alone" % name)
            continue
        for suffix in ("fas", "nhx_tree"):
            a = os.path.join(alone[k], "out." + suffix)
            i = os.path.join(indep, "out.%d.%s" % (k, suffix))
            if not filecmp.cmp(a, i, shallow=False):
                failures.append("%s: out.%d.%s differs from the run of it alone"
                                % (name, k, suffix))
                break
        else:
            n_identical += 1

    again = runs[1]
    produced = sorted(f for f in os.listdir(indep) if f.startswith("out."))
    if produced != sorted(f for f in os.listdir(again) if f.startswith("out.")):
        failures.append("the repeated independent run wrote other files")
    else:
        for f in produced:
            a, b = os.path.join(indep, f), os.path.join(again, f)
            if f == "out.queries.tsv":
                with open(a) as fa, open(b) as fb:
                    same = (fa.read().replace(indep, "") ==
                            fb.read().replace(again, ""))
            else:
                same = filecmp.cmp(a, b, shallow=False)
            if not same:
                failures.append("the repeated independent run differs in %s" % f)

    moved = []
    rc, batch = run_pagan2(binary, workdir, "batch", "queries.fas")
    batch_fas = os.path.join(batch, "out.fas")
    if rc != 0 or not os.path.isfile(batch_fas):
        failures.append("default run: exit %d, no out.fas" % rc)
    else:
        batch_rows = _read_fasta(batch_fas)
        moved = [n for k, n in enumerate(names, 1)
                 if os.path.isfile(os.path.join(alone[k], "out.fas"))
                 and batch_rows.get(n) != _read_fasta(os.path.join(alone[k], "out.fas")).get(n)]
        if not moved:
            failures.append("PREMISE: in the default run every row equals the row "
                            "alone, so the data shows nothing about batches")

    rc, d = run_pagan2(binary, workdir, "no_reference", "queries.fas",
                       ["--independent-queries"], reference=False)
    if rc == 0:
        failures.append("--independent-queries without a reference alignment "
                        "exited 0 (see %s/log)" % d)

    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and logs kept in %s)" % workdir)
        return 1
    print("PASS: %d of %d queries placed, each byte-identical to its run alone; "
          "the unrelated one unplaced; a repeat identical; the default run moved "
          "%d row(s)" % (n_identical, len(names), len(moved)))
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
