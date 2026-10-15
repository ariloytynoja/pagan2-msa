#!/usr/bin/env python3
"""A query placed onto a node made for an earlier query keeps every base.

Usage:
    python3 tests/test_placement_on_a_new_node_keeps_bases.py PAGAN2

When a query is placed, Reads_aligner::correct_sites_index() re-points the
parent of the target node to the new node, shifting each parent index by the
number of new columns before it. The node made for an earlier query of the
same run has insertion columns its parent does not point at yet ("corrected
later"). v1.6 summed the shift only at the columns the parent points at, so
new columns placed next to such an insertion were never counted: every later
parent index pointed too early, that clade was shifted against the rest of
the alignment, and its last columns fell off the 3' end. A batch of two
real 3.4-3.8 kb reads lost 168 nt of the 3' end of the reference row and
of a read in this way; in 40-read batches 41 of 55 rows lost 125-459 nt.

THE DATA (deterministic, standard library only, written under TMPDIR)
  Four codon sequences (ATG + 149 sense codons) on a four-taxon tree, and
  two queries in one file:
    q1  taxon A from codon 31 on (5' partial), with a 10-codon insertion
        after codon 60;
    q2  all of A, with q1's insertion and 3 more codons inside it.
  q1 is placed next to A. q2 matches the node made for q1 best -- q1's
  insertion from q1, A's 5' end from A -- so it is placed onto that node,
  and its 3 codons sit beside insertion columns the parent does not point
  at. v1.6 drops the last 3 codons of A, q1 and q2.

WHAT IS CHECKED
  PREMISE: q2 is placed onto the node holding A and q1 (in the output tree
  q2's sister clade is exactly {A, q1}); otherwise the data no longer
  reaches the code under test and the test says so.
  Every row of the output -- the four references and both queries -- holds
  exactly the bases it was given, gaps removed; and the four references,
  given aligned without gaps, still face each other base for base (no
  column where only one of two references has a base).

Exit 0 and a PASS line, or exit 1 naming what failed.
"""
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


def _parse_newick(text):
    """Nested lists of leaf names; branch lengths and NHX comments dropped."""
    pos = [0]

    def node():
        if text[pos[0]] == "(":
            pos[0] += 1
            kids = [node()]
            while text[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            assert text[pos[0]] == ")"
            pos[0] += 1
            _skip_label()
            return kids
        start = pos[0]
        while text[pos[0]] not in ",():;[":
            pos[0] += 1
        name = text[start:pos[0]]
        _skip_label()
        return name

    def _skip_label():
        while pos[0] < len(text) and text[pos[0]] not in ",);":
            if text[pos[0]] == "[":
                pos[0] = text.index("]", pos[0])
            pos[0] += 1

    return node()


def _leaves(tree):
    if isinstance(tree, str):
        return {tree}
    out = set()
    for kid in tree:
        out |= _leaves(kid)
    return out


def _sister_of(tree, name):
    if isinstance(tree, str):
        return None
    for i, kid in enumerate(tree):
        if kid == name:
            others = tree[:i] + tree[i + 1:]
            return set().union(*[_leaves(o) for o in others])
        found = _sister_of(kid, name)
        if found is not None:
            return found
    return None


def make_data(workdir):
    rng = random.Random(20261014)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
    taxa = {"A": _mutate(rng, anc, 8), "B": _mutate(rng, anc, 8),
            "C": _mutate(rng, anc, 14), "D": _mutate(rng, anc, 14)}
    ins1 = [rng.choice(SENSE) for _ in range(10)]
    ins2 = [rng.choice(SENSE) for _ in range(3)]
    a = taxa["A"]
    queries = {"q1": a[30:60] + ins1 + a[60:],
               "q2": a[:60] + ins1[:5] + ins2 + ins1[5:] + a[60:]}
    _write_fasta(os.path.join(workdir, "ref.aln.fas"),
                 [(t, "".join(c)) for t, c in sorted(taxa.items())])
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(TREE)
    _write_fasta(os.path.join(workdir, "queries.fas"),
                 [(n, "".join(queries[n])) for n in ("q1", "q2")])
    given = {n: "".join(c) for n, c in taxa.items()}
    given.update({n: "".join(c) for n, c in queries.items()})
    return given


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(argv[1])
    workdir = tempfile.mkdtemp(prefix="pagan2_new_node_bases.")
    given = make_data(workdir)
    os.makedirs(os.path.join(workdir, "tmp"))
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(workdir, "out"),
           "--temp-folder", os.path.join(workdir, "tmp")] + OPTIONS
    with open(os.path.join(workdir, "log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, cwd=workdir).returncode
    out = os.path.join(workdir, "out.fas")
    tree_file = os.path.join(workdir, "out.nhx_tree")
    if rc != 0 or not os.path.isfile(out) or not os.path.isfile(tree_file):
        print("FAIL: %s exited %d, or wrote no out.fas / out.nhx_tree "
              "(data and log kept in %s)" % (binary, rc, workdir))
        return 1

    failures = []
    with open(tree_file) as fh:
        sister = _sister_of(_parse_newick(fh.read().strip()), "q2")
    if sister != {"A", "q1"}:
        failures.append("PREMISE: q2 was not placed onto the node holding A "
                        "and q1 (its sister clade is %r), so the data no "
                        "longer reaches correct_sites_index()'s case"
                        % (sorted(sister) if sister else sister,))

    rows = _read_fasta(out)
    for name in sorted(given):
        if name not in rows:
            failures.append("%s: no row in out.fas" % name)
            continue
        bases = rows[name].replace("-", "").upper()
        if bases != given[name]:
            if given[name].startswith(bases):
                what = "its last %d bases are missing" % (len(given[name]) - len(bases))
            else:
                what = "%d bases written for %d given, not a prefix" % (
                    len(bases), len(given[name]))
            failures.append("%s: %s" % (name, what))

    # The reference alignment was given without gaps, so any two reference
    # rows must still face each other base for base: a column where only one
    # of them has a base means a clade was shifted against the other.
    refs = [t for t in ("A", "B", "C", "D") if t in rows]
    for i, x in enumerate(refs):
        for y in refs[i + 1:]:
            lone = sum(1 for p, q in zip(rows[x], rows[y]) if (p == "-") != (q == "-"))
            if lone:
                failures.append("%s and %s: %d columns where only one has a base "
                                "(they were given aligned without gaps)" % (x, y, lone))

    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and log kept in %s)" % workdir)
        return 1
    print("PASS: q2 placed onto the node made for q1; all %d rows keep every "
          "base" % len(given))
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
