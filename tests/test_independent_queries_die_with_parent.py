#!/usr/bin/env python3
"""--independent-queries: killing the parent pagan2 kills its query process.

Usage:
    python3 tests/test_independent_queries_die_with_parent.py PAGAN2

With --independent-queries each query is aligned in a child forked from the
parent pagan2. A caller that gives up on a run usually kills only the pid it
started -- a timeout that SIGKILLs the parent, a batch system's kill, an
operator -- and then removes the run's output and --temp-folder. A query
child that outlives its parent goes on computing at full speed for nobody,
into a directory that is being deleted.

THE DATA (deterministic, standard library only, written under TMPDIR)
  Codon sequences of 300 codons (ATG + 299 sense codons) on a balanced tree
  of 8 taxa, and four queries built alike (same length, same number of
  changes), run with --no-anchors --no-preselection --own-placement: every
  node of the tree is aligned against the query by full dynamic programming,
  so a query takes seconds at a peak memory of ~120 MB. These options also
  keep Exonerate out of the queries (only v1.6's argument-less presence probe
  runs it): with --no-anchors and --no-preselection alone it still ranks the
  candidate nodes, and its command line names files in the --temp-folder.
  The run uses a --temp-folder whose path is unique to this test; a forked
  child keeps its parent's command line and executable, so the pagan2
  processes of the run are those whose /proc/<pid>/cmdline holds that path
  and whose /proc/<pid>/exe is the binary under test (a helper program is
  not covered by the death signal, and is not counted).

WHAT IS CHECKED
  The parent is started and the first query's child (found by its parent
  pid in /proc/*/stat) is timed from its start to its manifest row (T). As
  soon as the second query's child exists, the PARENT pid alone is killed
  with SIGKILL. Within W = min(T/3, 0.5 s) no pagan2 process of the run may
  be left (a zombie counts as gone). The second query is like the first, so
  a child that ignored its parent's death would still have about 2/3 of T to
  run when W is over; one killed with its parent is gone in milliseconds. W
  follows the host's measured speed, so the verdict does not depend on it.
  The manifest the killed run left must not end with the "# complete" line.
  PREMISE: W >= 0.2 s (T >= 0.6 s). A host so fast that a query is over
  sooner gets a larger tree -- 16, then 32 taxa, each doubling T -- and only
  if that is still too fast does the test fail, saying nothing was tested.

Linux only (it reads /proc); elsewhere it prints SKIP and exits 0.
Exit 0 and a PASS line, or exit 1 naming what failed.
"""
import os
import random
import shutil
import signal
import subprocess
import sys
import tempfile
import time

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS]
N_CODONS = 300
TAXA_PER_ATTEMPT = (8, 16, 32)
N_QUERIES = 4
OPTIONS = ["--codons", "--no-anchors", "--no-preselection", "--own-placement",
           "--independent-queries"]
# The window in which the killed parent's child must be gone: a third of
# the measured first query, at most MAX_WINDOW, at least MIN_WINDOW.
WINDOW_FRACTION = 1.0 / 3.0
MAX_WINDOW_SECONDS = 0.5
MIN_WINDOW_SECONDS = 0.2
STARTUP_LIMIT_SECONDS = 300.0


def _mutate(rng, cod, n):
    cod = list(cod)
    for i in rng.sample(range(1, len(cod)), n):
        cod[i] = rng.choice([c for c in SENSE if c != cod[i]])
    return cod


def _tree(names):
    if len(names) == 1:
        return names[0] + ":0.05"
    half = len(names) // 2
    return "(%s,%s):0.02" % (_tree(names[:half]), _tree(names[half:]))


def make_data(workdir, n_taxa):
    rng = random.Random(20261014)
    anc = ["ATG"] + [rng.choice(SENSE) for _ in range(N_CODONS - 1)]
    names = ["T%d" % i for i in range(1, n_taxa + 1)]
    taxa = dict((t, _mutate(rng, anc, N_CODONS // 12)) for t in names)
    with open(os.path.join(workdir, "ref.aln.fas"), "w") as fh:
        for t in names:
            fh.write(">%s\n%s\n" % (t, "".join(taxa[t])))
    with open(os.path.join(workdir, "ref.tre"), "w") as fh:
        fh.write(_tree(names)[:-len(":0.02")] + ";\n")
    with open(os.path.join(workdir, "queries.fas"), "w") as fh:
        for k in range(1, N_QUERIES + 1):
            q = _mutate(rng, taxa[names[k % n_taxa]], N_CODONS // 20)
            fh.write(">query_%d\n%s\n" % (k, "".join(q)))


def _stat(pid):
    """(state, ppid) of *pid*, or None once it is gone."""
    try:
        with open("/proc/%d/stat" % pid) as fh:
            stat = fh.read()
    except (FileNotFoundError, ProcessLookupError):
        return None
    rest = stat[stat.rindex(")") + 2:].split()
    return rest[0], int(rest[1])


def _alive(pid):
    """True for a process that exists and is not a zombie."""
    st = _stat(pid)
    return st is not None and st[0] != "Z"


def _children(pid):
    """Live children of *pid*, by the parent pid in /proc/*/stat."""
    kids = []
    for entry in os.listdir("/proc"):
        if entry.isdigit():
            st = _stat(int(entry))
            if st is not None and st[1] == pid and st[0] != "Z":
                kids.append(int(entry))
    return kids


def _pagan2_processes_of_run(marker, binary):
    """Live pagan2 processes whose command line holds *marker*."""
    found = []
    needle = marker.encode()
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open("/proc/%s/cmdline" % entry, "rb") as fh:
                if needle not in fh.read():
                    continue
            if os.path.realpath("/proc/%s/exe" % entry) != binary:
                continue
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        if _alive(int(entry)):
            found.append(int(entry))
    return found


def _all_processes_of_run(marker):
    """Live pids whose command line holds *marker*: pagan2 and its helpers."""
    found = []
    needle = marker.encode()
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open("/proc/%s/cmdline" % entry, "rb") as fh:
                if needle in fh.read() and _alive(int(entry)):
                    found.append(int(entry))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return found


def _manifest_lines(path):
    try:
        with open(path) as fh:
            return [l for l in fh.read().split("\n")[1:] if l]
    except FileNotFoundError:
        return []


def attempt(binary, workdir, n_taxa):
    """One run: ("short", T) when a query is too fast to tell, else
    ("done", T, window, failures)."""
    make_data(workdir, n_taxa)
    marker = os.path.join(workdir, "tmp_%d_%d_%d" % (os.getpid(), n_taxa, int(time.time())))
    os.makedirs(marker)
    manifest = os.path.join(workdir, "out.queries.tsv")
    if os.path.exists(manifest):
        os.unlink(manifest)
    cmd = [binary, "-a", os.path.join(workdir, "ref.aln.fas"),
           "-t", os.path.join(workdir, "ref.tre"),
           "-q", os.path.join(workdir, "queries.fas"),
           "-o", os.path.join(workdir, "out"),
           "--temp-folder", marker] + OPTIONS
    failures = []
    with open(os.path.join(workdir, "log"), "w") as log:
        parent = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                  stdin=subprocess.DEVNULL, cwd=workdir)
    try:
        deadline = time.monotonic() + STARTUP_LIMIT_SECONDS
        first_child = None
        while first_child is None and time.monotonic() < deadline:
            if parent.poll() is not None:
                break
            kids = _children(parent.pid)
            if kids:
                first_child = min(kids)
                first_started = time.monotonic()
            else:
                time.sleep(0.01)
        query_seconds = None
        while first_child is not None and time.monotonic() < deadline:
            if parent.poll() is not None:
                break
            if _manifest_lines(manifest):
                query_seconds = time.monotonic() - first_started
                break
            time.sleep(0.01)
        if query_seconds is None:
            failures.append("never saw the first query finish (exit status %r; "
                            "see %s/log)" % (parent.poll(), workdir))
            return ("done", None, None, failures)
        window = min(query_seconds * WINDOW_FRACTION, MAX_WINDOW_SECONDS)
        if window < MIN_WINDOW_SECONDS:
            return ("short", query_seconds)
        second_child = None
        while time.monotonic() < deadline:
            if parent.poll() is not None:
                break
            kids = [k for k in _children(parent.pid) if k != first_child]
            if kids:
                second_child = kids[0]
                break
            time.sleep(0.005)
        if second_child is None:
            failures.append("never saw the second query's process start (exit "
                            "status %r; see %s/log)" % (parent.poll(), workdir))
            return ("done", query_seconds, window, failures)
        os.kill(parent.pid, signal.SIGKILL)
        parent.wait()
        killed_at = time.monotonic()
        left = _pagan2_processes_of_run(marker, binary)
        while left and time.monotonic() - killed_at < window:
            time.sleep(0.01)
            left = _pagan2_processes_of_run(marker, binary)
        if left:
            failures.append("%d pagan2 process(es) of the run still running "
                            "%.2f s after the parent was killed (a query takes "
                            "%.2f s here): %s"
                            % (len(left), window, query_seconds, left))
        tail = _manifest_lines(manifest)
        if tail and tail[-1].startswith("# complete"):
            failures.append("the killed run's manifest ends with %r" % tail[-1])
        return ("done", query_seconds, window, failures)
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
        # Leave nothing burning CPU: every process still carrying this run's
        # unique temp folder in its command line is one of ours.
        for pid in _all_processes_of_run(marker):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    if not sys.platform.startswith("linux"):
        print("SKIP: reads /proc, Linux only")
        return 0
    binary = os.path.realpath(argv[1])
    workdir = tempfile.mkdtemp(prefix="pagan2_die_with_parent.")
    failures = []
    tried = []
    result = None
    for n_taxa in TAXA_PER_ATTEMPT:
        result = attempt(binary, workdir, n_taxa)
        tried.append((n_taxa, result[1]))
        if result[0] == "done":
            failures = result[3]
            break
    else:
        failures.append("PREMISE: even with %d taxa a query took only %.2f s, so a "
                        "child killed with its parent cannot be told from one that "
                        "finished its query; nothing was tested (tried %s)"
                        % (tried[-1][0], tried[-1][1], tried))

    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        print("(data and log kept in %s)" % workdir)
        return 1
    print("PASS: the parent killed while query 2 ran (%d taxa; query 1 took "
          "%.2f s); no pagan2 process of the run left within %.2f s"
          % (tried[-1][0], result[1], result[2]))
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
