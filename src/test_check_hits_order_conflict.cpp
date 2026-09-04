// Synthetic, deterministic reproducer for the check_hits_order_conflict()
// std::out_of_range abort.
//
// The function walks [start_site, start_site+length) over hit_site1/
// hit_site2, two vector<bool> sized to each sequence's length, using .at().
// An Exonerate hit whose start_site lands at (or past) the end of the
// sequence therefore throws, and nothing catches it: the whole batch of
// queries being aligned dies, not just the one query whose anchor was out
// of range.
//
// This calls the function directly with a hand-built Substring_hit rather
// than going through Exonerate, because the trigger is a specific hit
// coordinate rather than any general property of the input sequences --
// driving it from a FASTA would depend on Exonerate's own behaviour and
// would not be a stable regression test.  The boundary chosen here,
// start_site_1 == len1, is exactly one past the last valid 0-based index.
//
// Build (from src/, after a normal build has produced the listed .o files):
//   g++ -std=c++11 -w -Iutils -Imain -I. -o test_check_hits_order_conflict \
//     test_check_hits_order_conflict.cpp find_anchors.o settings.o \
//     settings_handle.o log_output.o text_utils.o check_version.o \
//     -lboost_program_options -lboost_regex -lboost_thread -lboost_system \
//     -lgomp -lm -lz -lpthread -ldl
// Run: ./test_check_hits_order_conflict
//   Exit 0 + two PASS lines on the fixed code; SIGABRT (std::out_of_range)
//   on the pre-fix code, from the [boundary] scenario alone.
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include "utils/find_anchors.h"
#include "utils/settings_handle.h"

using namespace ppa;
using namespace std;

static void die(const string &msg)
{
    cerr << "FAIL: " << msg << endl;
    exit(1);
}

int main(int argc, char *argv[])
{
    // Populate Settings_handle::st with defaults (exonerate-hit-trim=5 etc.)
    // the same way main() does, without requiring any of the normally
    // mandatory alignment-input flags.
    const char *fake_argv[] = {"test_check_hits_order_conflict"};
    Settings_handle::st.read_command_line_arguments(1, const_cast<char**>(fake_argv));
    int trim = Settings_handle::st.get("exonerate-hit-trim").as<int>();

    string seq1(40, 'A');
    string seq2(40, 'A');
    int len1 = (int)seq1.length();
    int len2 = (int)seq2.length();

    // ---- Scenario 1: the out-of-range boundary hit ----
    // start_site_1 == len1 (one past the last valid 0-based index).  length
    // must survive "length -= trim*2" and still be positive, or the inner
    // loop's own bound (start_site_1+length) never exceeds start_site_1 and
    // the loop body -- where the .at() call that actually throws lives --
    // never runs even with the bug present.  A length of, say, 1 makes this
    // test pass against the buggy code too, i.e. proves nothing; keep it
    // comfortably above 2*exonerate-hit-trim.
    {
        vector<Substring_hit> hits;
        Substring_hit h;
        h.start_site_1 = len1;
        h.start_site_2 = 5;
        h.length = 20;
        h.score = 100;
        hits.push_back(h);

        cout << "[boundary] start_site_1=" << h.start_site_1 << " len1=" << len1 << endl;
        Find_anchors fa;
        fa.check_hits_order_conflict(&seq1, &seq2, &hits);
        cout << "PASS: boundary hit handled without throwing." << endl;
    }

    // ---- Scenario 2: a normal, comfortably in-bounds hit is preserved ----
    // Confirms the fix does not over-clamp: a hit that already fits inside
    // [0, len) after trimming must survive with its trimmed coordinates
    // intact, not be zeroed out by the new length-clamp.
    {
        vector<Substring_hit> hits;
        Substring_hit h;
        h.start_site_1 = 10;
        h.start_site_2 = 10;
        h.length = 30;
        h.score = 100;
        hits.push_back(h);

        Find_anchors fa;
        fa.check_hits_order_conflict(&seq1, &seq2, &hits);

        if (hits.size() != 1)
            die("normal in-bounds hit was dropped, expected 1 survivor");
        const Substring_hit &out = hits.front();
        int expect_start = 10 + trim;
        int expect_length = 30 - trim*2;
        if (out.start_site_1 != expect_start || out.start_site_2 != expect_start)
            die("normal hit's start was altered beyond the documented trim");
        if (out.length != expect_length)
            die("normal hit's length was altered beyond the documented trim");
        cout << "PASS: normal in-bounds hit preserved (start=" << out.start_site_1
             << " length=" << out.length << ")." << endl;
    }

    return 0;
}
