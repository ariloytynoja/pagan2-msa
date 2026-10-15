#!/usr/bin/env python3
"""--independent-queries: what the manifest says is what is on disk.

Usage:
    python3 tests/test_independent_queries_manifest.py PAGAN2

"<outfile>.queries.tsv" has a header, one row of five tab-separated fields
per query (query, name, status, detail, alignment) and, only once every
query has a row, the trailer "# complete: N queries, P placed, U unplaced,
F failed". A reader takes "placed" to mean the named alignment file is this
run's and complete, so the parent checks the file rather than trust the
child, and removes what a failed query left.

THE DATA (deterministic, standard library only, written under TMPDIR)
  Four codon sequences (ATG + 149 sense codons) on a four-taxon tree. The
  query file holds five records:
    query_1        a few substitutions against taxon A        -> n = 1
    query<TAB>2    the same against B, a tab in its name      -> n = 2
    empty          a header and no sequence: the reader drops it
    unrelated      random codons, placed nowhere              -> n = 3
    query_4        the same against C                         -> n = 4
  Before the run the output directory holds files of an "earlier run":
  out.3.fas and out.3_contigs.fas (n = 3 is unplaced now), out.7.fas and
  out.03.fas (names this run does not write), and out.4.fas as a symbolic
  link to /dev/full, so query 4's child writes its alignment nowhere and
  still reports it placed.

WHAT IS CHECKED
  1. every row has five fields; row 2's name is "query 2" (the tab written
     as a space); n counts the records READ: 1-4 for the four with a
     sequence, and the run warns once that 1 of the 5 records was dropped;
  2. query 4 is "failed", detail "no alignment file written" and naming the
     link it could not remove; its out.4.nhx_tree is gone;
  3. the earlier run's out.3.fas and out.3_contigs.fas are gone, out.7.fas
     and out.03.fas untouched;
  4. the manifest ends with "# complete: 4 queries, 2 placed, 1 unplaced,
     1 failed", which holds no tab; the run exits 1 (a query failed);
  5. FORK FAILURE: run under RLIMIT_NPROC = 1, so no fork() can succeed.
     Query 1 is "failed" with "fork: ..." and "no later query is run", every
     later query "failed" with "not run: ...", the trailer still closes the
     manifest, and the run exits 1. Where the limit cannot make fork() fail
     (running as root) this part prints a SKIP line.

Exit 0 and a PASS line, or exit 1 naming what failed.
"""
import os
import random
import resource
import shutil
import subprocess
import sys
import tempfile

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS]
N_CODONS = 150
TREE = "((A:0.05,B:0.05):0.02,(C:0.05,D:0.05):0.02);\n"
OPTIONS = ["--guidetree", "--codons", "--one-placement-only",
           "--independent-queries"]


def _mutate(rng, cod, n):
    cod = list(cod)
    for i in rng.sample(range(1, len(cod)), n):
        cod[i] = rng.choice([c for c in SENSE if c != cod[i]])
    return cod


def make_data(workdir):
    rng = random.Random(20261015)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
    taxa = {"A": _mutate(rng, anc, 8), "B": _mutate(rng, anc, 8),
            "C": _mutate(rng, anc, 14), "D": _mutate(rng, anc, 14)}
    with open(os.path.join(workdir, "ref.aln.fas"), "w") as fh:
        for t in sorted(taxa):
            fh.write(">%s\n%s\n" % (t, "".join(taxa[t])))
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)
    with open(os.path.join(workdir, "queries.fas"), "w") as fh:
        fh.write(">query_1\n%s\n" % "".join(_mutate(rng, taxa["A"], 3)))
        fh.write(">query\t2\n%s\n" % "".join(_mutate(rng, taxa["B"], 3)))
        fh.write(">empty\n")
        fh.write(">unrelated\n%s\n" % "".join(rng.choice(SENSE) for _ in range(N_CODONS)))
        fh.write(">query_4\n%s\n" % "".join(_mutate(rng, taxa["C"], 3)))


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


def run_pagan2(binary, workdir, out_dir, preexec_fn=None):
    os.makedirs(os.path.join(out_dir, "tmp"), exist_ok=True)
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(out_dir, "out"),
           "--temp-folder", os.path.join(out_dir, "tmp")] + OPTIONS
    log_path = os.path.join(out_dir, "log")
    with open(log_path, "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=out_dir,
                            preexec_fn=preexec_fn).returncode
    with open(log_path, errors="replace") as fh:
        return rc, fh.read()


def check_main_run(binary, workdir, failures):
    out_dir = os.path.join(workdir, "run")
    os.makedirs(out_dir)
    earlier = {"out.3.fas": ">stale\nACGT\n", "out.3_contigs.fas": ">stale\nACGT\n",
               "out.7.fas": ">other\nACGT\n", "out.03.fas": ">other\nACGT\n"}
    for name, text in earlier.items():
        with open(os.path.join(out_dir, name), "w") as fh:
            fh.write(text)
    has_dev_full = os.path.exists("/dev/full")
    if has_dev_full:
        os.symlink("/dev/full", os.path.join(out_dir, "out.4.fas"))

    rc, log = run_pagan2(binary, workdir, out_dir)
    manifest = os.path.join(out_dir, "out.queries.tsv")
    if not os.path.isfile(manifest):
        failures.append("no manifest (exit %d; see %s/log)" % (rc, out_dir))
        return None
    rows, trailer, problems = read_manifest(manifest)
    failures.extend(problems)
    got = [(r[0], r[1]) for r in rows]
    want = [("1", "query_1"), ("2", "query 2"), ("3", "unrelated"), ("4", "query_4")]
    if got != want:
        failures.append("manifest rows %r, expected %r" % (got, want))
    warned = [l for l in log.split("\n") if "records without a sequence" in l]
    if len(warned) != 1 or "query file: 1 of 5, dropped" not in warned[0]:
        failures.append("expected one warning that 1 of the 5 records was dropped, "
                        "got %r" % warned)

    by_n = dict((r[0], r) for r in rows)
    for n in ("1", "2"):
        r = by_n.get(n)
        if r and (r[2] != "placed" or r[4] != os.path.join(out_dir, "out.%s.fas" % n)
                  or not os.path.getsize(r[4])):
            failures.append("query %s: %r, expected placed with its alignment" % (n, r))
    r = by_n.get("3")
    if r and (r[2] != "unplaced" or r[4] != ""):
        failures.append("query 3: %r, expected unplaced without a file" % (r,))
    for name in ("out.3.fas", "out.3_contigs.fas"):
        if os.path.lexists(os.path.join(out_dir, name)):
            failures.append("%s of the earlier run is still there" % name)
    for name in ("out.7.fas", "out.03.fas"):
        if not os.path.isfile(os.path.join(out_dir, name)):
            failures.append("%s, which this run does not write, was removed" % name)

    n_failed = 0
    r = by_n.get("4")
    if has_dev_full and r:
        n_failed = 1
        link = os.path.join(out_dir, "out.4.fas")
        if r[2] != "failed" or not r[3].startswith("no alignment file written") \
                or link not in r[3] or r[4] != "":
            failures.append("query 4 (alignment written to /dev/full): %r, expected "
                            "failed, 'no alignment file written', naming %s" % (r, link))
        if os.path.lexists(os.path.join(out_dir, "out.4.nhx_tree")):
            failures.append("query 4 failed but its out.4.nhx_tree is still there")
    elif r and r[2] != "placed":
        failures.append("query 4: %r, expected placed" % (r,))

    expected = "# complete: 4 queries, %d placed, 1 unplaced, %d failed" % (3 - n_failed, n_failed)
    if trailer != expected or (trailer and "\t" in trailer):
        failures.append("manifest ends with %r, expected %r" % (trailer, expected))
    if rc != (1 if n_failed else 0):
        failures.append("exit status %d with %d failed queries" % (rc, n_failed))
    return has_dev_full


def check_fork_failure(binary, workdir, failures):
    if os.geteuid() == 0 or not hasattr(resource, "RLIMIT_NPROC"):
        return "SKIP (fork failure): RLIMIT_NPROC cannot make fork() fail here"
    out_dir = os.path.join(workdir, "nofork")
    os.makedirs(out_dir)

    def no_more_processes():
        resource.setrlimit(resource.RLIMIT_NPROC, (1, 1))

    rc, log = run_pagan2(binary, workdir, out_dir, preexec_fn=no_more_processes)
    manifest = os.path.join(out_dir, "out.queries.tsv")
    if not os.path.isfile(manifest):
        failures.append("fork failure: no manifest (exit %d; see %s/log)" % (rc, out_dir))
        return None
    rows, trailer, problems = read_manifest(manifest)
    failures.extend("fork failure: " + p for p in problems)
    if rows and rows[0][2] != "failed" and "fork" not in rows[0][3]:
        return "SKIP (fork failure): fork() still succeeds under RLIMIT_NPROC = 1"
    if len(rows) != 4:
        failures.append("fork failure: %d rows, expected 4" % len(rows))
        return None
    first = rows[0]
    if first[2] != "failed" or not first[3].startswith("fork: ") \
            or "no later query is run" not in first[3]:
        failures.append("fork failure: query 1 %r, expected failed, 'fork: ...', "
                        "'no later query is run'" % (first,))
    for r in rows[1:]:
        if r[2] != "failed" or not r[3].startswith("not run: fork: "):
            failures.append("fork failure: query %s %r, expected failed, 'not run: "
                            "fork: ...'" % (r[0], r))
    expected = "# complete: 4 queries, 0 placed, 0 unplaced, 4 failed"
    if trailer != expected:
        failures.append("fork failure: manifest ends with %r, expected %r"
                        % (trailer, expected))
    if rc != 1:
        failures.append("fork failure: exit status %d, expected 1" % rc)
    return None


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    workdir = tempfile.mkdtemp(prefix="pagan2_independent_manifest.")
    make_data(workdir)
    failures = []
    notes = []
    if not check_main_run(binary, workdir, failures):
        notes.append("SKIP (placed without a file): no /dev/full here")
    skip = check_fork_failure(binary, workdir, failures)
    if skip:
        notes.append(skip)
    for n in notes:
        print(n)
    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and logs kept in %s)" % workdir)
        return 1
    print("PASS: rows of five fields, the tab written as a space, n counting the "
          "records read; a placed query without its file failed and cleaned up; "
          "an earlier run's files for this run's n removed; the trailer; a fork "
          "failure stops the run")
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
