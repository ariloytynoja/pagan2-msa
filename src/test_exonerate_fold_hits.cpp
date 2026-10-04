// Synthetic test for Exonerate_queries::fold_hits(): one exonerate run's hits
// folded into one hit per (query, target), independently of the order
// exonerate printed them in.
//
// exonerate prints the same set of hits in a different order from run to run.
// The fold used to walk them in that order: same-strand hits were summed into
// the stored one, and a single opposite-strand hit replaced it when it scored
// more than the stored (possibly partial) sum. So +5, +6, -9 kept the plus
// strand at 11 but +5, -9, +6 kept the minus strand at 9.
//
// What is pinned:
// 1. Every one of the six orders of {+5, +6, -9} on one target gives the same
//    folded hit: the plus strand pair, score 11, spanning both plus hits.
// 2. Hits on other targets and queries are folded separately.
// 3. A tie between two strand pairs is decided by canonical_order(), i.e. by
//    the hits, not by which was printed first.
// 4. The result is ranked best first.
//
// Build (from src/, after a normal build has produced the objects; every
// object but main.o):
//   g++ -std=c++11 -w -fopenmp -I. -Iutils -Imain -o test_exonerate_fold_hits \
//     test_exonerate_fold_hits.cpp $(ls *.o | grep -v '^main.o$') \
//     -lboost_program_options -lboost_regex -lboost_thread -lboost_system \
//     -lgomp -lm -lz -lpthread -ldl
// Run: ./test_exonerate_fold_hits  (exit 0 and PASS lines, or FAIL + 1)
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "utils/exonerate_queries.h"

using namespace ppa;
using namespace std;

static int failures = 0;

static void check(bool ok, const string &what, const string &detail = "")
{
    cout << (ok ? "PASS " : "FAIL ") << what << endl;
    if (!detail.empty())
        cout << "      " << detail << endl;
    if (!ok)
        failures++;
}

static hit mk(const string &query, const string &node, int score,
              int q_start, int q_end, char q_strand,
              int t_start, int t_end, char t_strand)
{
    hit h;
    h.query = query; h.node = node; h.score = score;
    h.q_start = q_start; h.q_end = q_end; h.q_strand = q_strand;
    h.t_start = t_start; h.t_end = t_end; h.t_strand = t_strand;
    return h;
}

static string show(const hit &h)
{
    stringstream s;
    s << h.query << "/" << h.node << " " << h.q_strand << h.t_strand
      << " score " << h.score << " q " << h.q_start << "-" << h.q_end
      << " t " << h.t_start << "-" << h.t_end;
    return s.str();
}

static bool same(const hit &a, const hit &b)
{
    return a.query == b.query && a.node == b.node && a.score == b.score
        && a.q_start == b.q_start && a.q_end == b.q_end && a.q_strand == b.q_strand
        && a.t_start == b.t_start && a.t_end == b.t_end && a.t_strand == b.t_strand;
}

int main()
{
    // 1.
    vector<hit> three;
    three.push_back(mk("r", "A", 5, 0, 40, '+', 100, 140, '+'));
    three.push_back(mk("r", "A", 6, 60, 110, '+', 160, 210, '+'));
    three.push_back(mk("r", "A", 9, 0, 70, '+', 300, 230, '-'));
    hit want = mk("r", "A", 11, 0, 110, '+', 100, 210, '+');

    vector<int> idx;
    for (int i = 0; i < 3; i++)
        idx.push_back(i);
    int orders = 0, agree = 0;
    string first_wrong;
    do
    {
        vector<hit> in;
        for (int i = 0; i < 3; i++)
            in.push_back(three[idx[i]]);
        vector<hit> out = Exonerate_queries::fold_hits(in);
        orders++;
        if (out.size() == 1 && same(out[0], want))
            agree++;
        else if (first_wrong.empty())
            first_wrong = out.empty() ? string("nothing") : show(out[0]);
    }
    while (next_permutation(idx.begin(), idx.end()));
    check(orders == 6 && agree == 6,
          "all 6 orders of +5, +6, -9 fold to the plus strand pair at 11",
          agree == 6 ? show(want) : "one order gave " + first_wrong);

    // 2. and 4.
    vector<hit> mixed(three);
    mixed.push_back(mk("r", "B", 20, 5, 95, '+', 5, 95, '+'));
    mixed.push_back(mk("s", "A", 3, 0, 30, '+', 0, 30, '+'));
    vector<hit> out = Exonerate_queries::fold_hits(mixed);
    check(out.size() == 3 && out[0].node == "B" && out[0].score == 20
              && same(out[1], want) && out[2].query == "s" && out[2].score == 3,
          "targets and queries fold separately, best first");

    // 3. A tie: 7 on (+,+) against 3+4 on (+,-), in both print orders.
    vector<hit> tie;
    tie.push_back(mk("r", "C", 7, 10, 50, '+', 10, 50, '+'));
    tie.push_back(mk("r", "C", 3, 0, 20, '+', 90, 70, '-'));
    tie.push_back(mk("r", "C", 4, 30, 50, '+', 60, 40, '-'));
    vector<hit> tie_rev(tie.rbegin(), tie.rend());
    vector<hit> a = Exonerate_queries::fold_hits(tie);
    vector<hit> b = Exonerate_queries::fold_hits(tie_rev);
    check(a.size() == 1 && b.size() == 1 && same(a[0], b[0]) && a[0].score == 7,
          "a tie between strand pairs does not depend on the print order",
          a.empty() ? "" : show(a[0]));

    return failures ? 1 : 0;
}
