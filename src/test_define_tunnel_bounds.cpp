// Synthetic reproducer for the Find_anchors::define_tunnel() out_of_range
// abort.
//
// index1/index2 hold the GAPPED position of each NON-GAP character (so
// index1.size() == count of non-gap chars in str1), but start_site_1 is a
// count in that same ungapped space -- an Exonerate hit whose reported
// [start_site_1, start_site_1+length) range runs past the end of a
// sequence that carries internal gaps can have start_site_1+length exceed
// index1.size() even though it already passed check_hits_order_conflict's
// clamp (that clamp bounds against str1's GAPPED length, which is >=
// index1.size() whenever str1 contains '-'). This constructs exactly that:
// a short string with a gap near the end, and a hit that legitimately fits
// the gapped length but overruns the ungapped index.
//
// Build (from src/, after a normal build has produced the listed .o files):
//   g++ -std=c++11 -w -Iutils -Imain -I. -o test_define_tunnel_bounds \
//     test_define_tunnel_bounds.cpp find_anchors.o settings.o \
//     settings_handle.o log_output.o text_utils.o check_version.o \
//     -lboost_program_options -lboost_regex -lboost_thread -lboost_system \
//     -lgomp -lm -lz -lpthread -ldl
// Run: ./test_define_tunnel_bounds
//   Exit 0 + a PASS line on the fixed code; SIGABRT (std::out_of_range) on
//   the pre-fix code.
#include <iostream>
#include <vector>
#include <string>
#include "utils/find_anchors.h"
#include "utils/settings_handle.h"

using namespace ppa;
using namespace std;

int main()
{
    const char *fake_argv[] = {"test_define_tunnel_bounds"};
    Settings_handle::st.read_command_line_arguments(1, const_cast<char**>(fake_argv));

    // 10 real characters + 2 trailing gaps = gapped length 12, but only 10
    // non-gap positions -- index1/index2 will each have size 10.
    string str1 = "AAAAAAAAAA--";
    string str2 = "AAAAAAAAAA--";

    vector<Substring_hit> hits;
    Substring_hit h;
    // Fits inside the GAPPED length (12) but not the ungapped one (10):
    // start_site_1+length == 5+7 == 12 <= gapped length, yet walks
    // index1.at(5) .. index1.at(11), and index1 only has 10 entries.
    h.start_site_1 = 5;
    h.start_site_2 = 5;
    h.length = 7;
    h.score = 100;
    h.plus_strand_1 = true;
    h.plus_strand_2 = true;
    hits.push_back(h);

    vector<int> upper_bound, lower_bound;
    Find_anchors fa;
    fa.define_tunnel(&hits, &upper_bound, &lower_bound, &str1, &str2);

    cout << "PASS: define_tunnel returned without throwing." << endl;
    return 0;
}
