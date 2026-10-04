// Synthetic test for the RNG seed that pagan2's temp-file names come from.
//
// THE DEFECT.  Every helper wrapper names its scratch files from rand():
// Exonerate_queries writes "q<r>.fas"/"t<r>.fas", Mafft_alignment "m<r>.fas",
// and BppAncestors, Bppdist_tree, Raxml_tree and FastTree_tree do the same --
// all in the ONE directory get_temp_dir() returns (/tmp by default, or
// --temp-folder).  Each site already guards against reusing a name with a
// check-then-create retry loop, e.g. Exonerate_queries::write_exonerate_input():
//
//     while(true) { ... if(!q_file && !t_file) { open both; break; } *r = rand(); }
//
// That loop assumes two processes will walk DIFFERENT candidate names.  Under
// main.cpp's srand(time(0)) they do not: the seed is the clock in whole
// SECONDS, so two pagan2 processes started in the same second get the same
// seed, hence the same rand() sequence, hence the same first candidate and the
// same retries.  Both see the name free, both create it, both run their helper
// on whatever the other last wrote, and each delete_files() removes the
// other's files.  Nothing crashes -- the alignment is simply of the wrong
// sequences, or comes back empty and looks like an unlucky batch.
//
// It is not only a multi-process problem: pagan2 also runs threaded, and
// rand() is not thread-safe, so two threads can be handed the same value.
// Seeding does not fix that half; it fixes the process half, which is the one
// that bites whenever several pagan2 runs share a --temp-folder.
//
// WHY THIS TEST IS DETERMINISTIC.  It does not race two real processes and
// hope they land in the same second.  It applies each seeding RULE to a fixed
// timestamp and a set of pids, which is exactly the property in question and
// cannot flake:
//
//     old_seed(t, pid) = t                 -> identical for every pid
//     xor_seed(t, pid) = t ^ pid           -> distinct within one second, but
//                                             (t,pid) and (t+1,pid+1) collide
//                                             whenever t and pid are even
//     process_seed(t, pid)                 -> utils/process_seed.h, the rule
//                                             main() really uses
//
// and then shows the CONSEQUENCE the retry loop cares about: the sequence of
// candidate file names two processes would try.
//
// Build (from src/):
//   g++ -std=c++11 -w -I. -Iutils -Imain -o test_temp_file_seed_is_per_process \
//       test_temp_file_seed_is_per_process.cpp
// Run: ./test_temp_file_seed_is_per_process  (exit 0 and PASS lines, or FAIL + 1)
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "utils/process_seed.h"

using namespace std;
using namespace ppa;

static int failures = 0;

static void check(bool ok, const string &what)
{
    cout << (ok ? "PASS  " : "FAIL  ") << what << endl;
    if (!ok)
        ++failures;
}

// The two seeding rules, named so the test reads as the comparison it is.
static unsigned old_seed(time_t t, pid_t /*pid*/)
{
    return static_cast<unsigned>(t);
}

static unsigned xor_seed(time_t t, pid_t pid)
{
    return static_cast<unsigned>(t) ^ static_cast<unsigned>(pid);
}

static unsigned new_seed(time_t t, pid_t pid)
{
    return process_seed(static_cast<uint64_t>(t), static_cast<uint64_t>(pid));
}

// The first `n` candidate names a process would try, exactly as
// write_exonerate_input() generates them: seed once, then rand() per retry.
static vector<string> candidate_names(unsigned seed, int n)
{
    vector<string> names;
    srand(seed);
    int r = rand();
    for (int i = 0; i < n; ++i)
    {
        ostringstream os;
        os << "q" << r << ".fas";
        names.push_back(os.str());
        r = rand();
    }
    return names;
}

int main()
{
    // A fixed instant and several plausible pids. Same second, different
    // processes -- the situation the defect needs and the only one it needs.
    const time_t t = 1757000000;
    const pid_t pids[] = {4242, 4243, 9001, 31337};
    const int n_pids = sizeof(pids) / sizeof(pids[0]);

    // ---- PREMISE: the old rule really does collapse to one seed ----------
    {
        set<unsigned> seeds;
        for (int i = 0; i < n_pids; ++i)
            seeds.insert(old_seed(t, pids[i]));
        check(seeds.size() == 1,
              "PREMISE: srand(time(0)) gives every process in one second the "
              "SAME seed");
    }

    // ---- THE FIX: one seed per process ----------------------------------
    {
        set<unsigned> seeds;
        for (int i = 0; i < n_pids; ++i)
            seeds.insert(new_seed(t, pids[i]));
        check(static_cast<int>(seeds.size()) == n_pids,
              "process_seed(time, pid) gives each process its OWN seed");
    }

    // ---- THE CONSEQUENCE the retry loop depends on ----------------------
    // Two processes must not walk the same candidate names, or the
    // check-then-create loop can never diverge and both take name #1.
    {
        vector<string> a = candidate_names(old_seed(t, pids[0]), 5);
        vector<string> b = candidate_names(old_seed(t, pids[1]), 5);
        check(a == b,
              "PREMISE: under the old seed two processes try the IDENTICAL "
              "sequence of temp names, so the retry loop cannot help");
        cout << "        both would try: " << a[0] << ", " << a[1]
             << ", " << a[2] << " ..." << endl;
    }
    {
        vector<string> a = candidate_names(new_seed(t, pids[0]), 5);
        vector<string> b = candidate_names(new_seed(t, pids[1]), 5);
        check(a[0] != b[0],
              "under the new seed the two processes' FIRST candidate differs");
        set<string> all(a.begin(), a.end());
        all.insert(b.begin(), b.end());
        check(all.size() == a.size() + b.size(),
              "and none of their first five candidates coincide");
        cout << "        pid " << pids[0] << " tries " << a[0]
             << ", pid " << pids[1] << " tries " << b[0] << endl;
    }

    // ---- The seed is still clock-derived, so a rerun differs ------------
    // Guards against "fixing" this with a constant, which would make every
    // run of every process reuse one name.
    {
        check(new_seed(t, pids[0]) != new_seed(t + 1, pids[0]),
              "the seed still moves with the clock, so one process rerunning "
              "does not reuse its own names");
    }

    // ---- Nearby (time, pid) pairs: a launcher starting one a second -----
    // XOR collides along the diagonal (t+k, pid+k); the mixed seed must not
    // collide anywhere in a block of 64 seconds x 256 consecutive pids.
    {
        const pid_t p0 = 4242;
        check(xor_seed(t, p0) == xor_seed(t + 1, p0 + 1),
              "PREMISE: t ^ pid repeats for (t, pid) and (t+1, pid+1) when both are even");
        set<unsigned> seeds;
        int pairs = 0;
        for (int dt = 0; dt < 64; ++dt)
            for (int dp = 0; dp < 256; ++dp, ++pairs)
                seeds.insert(new_seed(t + dt, p0 + dp));
        ostringstream os;
        os << "no two of " << pairs << " nearby (time, pid) pairs share a seed ("
           << seeds.size() << " distinct)";
        check(static_cast<int>(seeds.size()) == pairs, os.str());
    }

    cout << (failures ? "\nFAILURES: " : "\nAll checks passed (")
         << (failures ? failures : 0) << (failures ? "\n" : ")\n");
    return failures ? 1 : 0;
}
