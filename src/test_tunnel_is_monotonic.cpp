// Synthetic, deterministic reproducer for a std::out_of_range abort in
// Tunnel_matrix::initTunnelEntries(), reached through the exonerate-anchored
// tunnel (Viterbi_alignment::define_tunnel(), non-NCBI path):
//
//     fa.check_hits_order_conflict(&s1,&s2,&hits);
//     fa.define_tunnel(&hits,&upper_bound,&lower_bound,&s1,&s2);
//     align_array M(left_length-1, right_length-1, upper_bound, lower_bound, &empty);
//
// Tunnel_matrix requires both bounds to be monotonically increasing (its own
// doc comment says so). define_tunnel() produces monotonic bounds as long as
// the hits it is given are co-linear -- start_site_2 increases with
// start_site_1 -- and check_hits_order_conflict() is what is meant to ensure
// that: it sorts the hits by start_site_1 and, for every adjacent pair whose
// start_site_2 go the other way, erases the lower-scoring one.
//
// THE DEFECT. When the erased hit is the EARLIER one of the pair, the loop
// carries on from the survivor without comparing it with the hit before the
// one it erased. Three hits are enough:
//
//     A  start_site_1 0    start_site_2 100   score 30
//     B  start_site_1 40   start_site_2 150   score 30
//     C  start_site_1 80   start_site_2 40    score 40
//
// A,B is in order; B,C is not and B scores lower, so B is erased -- and A,C,
// which also cross, are never compared. define_tunnel() then walks a diagonal
// that jumps from ~130 back down to ~45, the bounds go non-monotonic, and the
// Tunnel_matrix constructor throws. Nothing catches it, so the whole run
// aborts, not just the one query.
//
// WHAT IS CHECKED
//   [order]   the hits check_hits_order_conflict() keeps are co-linear;
//   [tunnel]  for those hits define_tunnel() returns monotonic bounds and a
//             Tunnel_matrix can be built from them (this is the abort);
//   [anchors] every cell an anchor lies on is inside the tunnel;
//   [guard]   define_tunnel() handed crossing hits directly (no caller in the
//             tree does this, but nothing in its signature prevents it) still
//             returns monotonic bounds, and they still contain every anchor
//             cell -- the repair may only WIDEN the tunnel, never cut an
//             anchor out of it.
//
// Build (from src/, after a normal build has produced the listed .o files):
//   g++ -std=c++11 -w -Iutils -Imain -I. -o test_tunnel_is_monotonic \
//     test_tunnel_is_monotonic.cpp find_anchors.o settings.o \
//     settings_handle.o log_output.o text_utils.o check_version.o \
//     -lboost_program_options -lboost_regex -lboost_thread -lboost_system \
//     -lgomp -lm -lz -lpthread -ldl
// Run: ./test_tunnel_is_monotonic
//   Exit 0 and four PASS lines on the fixed code; FAIL lines and exit 1 on
//   v1.6, where [tunnel] reports the std::out_of_range the run would die of.
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "utils/find_anchors.h"
#include "utils/settings_handle.h"
#include "utils/tunnel_matrix.h"

using namespace ppa;
using namespace std;

static int failures = 0;

static void report(bool ok, const string &what, const string &detail = "")
{
    cout << (ok ? "PASS  " : "FAIL  ") << what << endl;
    if(!detail.empty())
        cout << "        " << detail << endl;
    if(!ok)
        failures++;
}

static Substring_hit hit(int s1, int s2, int length, int score)
{
    Substring_hit h;
    h.start_site_1 = s1;
    h.start_site_2 = s2;
    h.length = length;
    h.score = score;
    return h;
}

// Index of the first column whose bound is lower than the previous one's, or
// -1 if both bounds are monotonically increasing.
static int first_descent(const vector<int> &upper, const vector<int> &lower)
{
    for(int i = 1; i < (int)upper.size(); i++)
        if(upper[i] < upper[i-1] || lower[i] < lower[i-1])
            return i;
    return -1;
}

// The cells define_tunnel() anchors on: hit position k of start_site_1+k /
// start_site_2+k sits at column index1[..] and row index2[..], i.e. at the
// 1-based gapped positions. Neither test sequence has gaps, so that is +1.
static string anchor_outside(const vector<Substring_hit> &hits,
                             const vector<int> &upper, const vector<int> &lower)
{
    for(size_t h = 0; h < hits.size(); h++)
        for(int k = 0; k < hits[h].length; k++)
        {
            int x = hits[h].start_site_1 + k + 1;
            int y = hits[h].start_site_2 + k + 1;
            if(y < upper.at(x) || y > lower.at(x))
            {
                stringstream s;
                s << "anchor cell (" << x << "," << y << ") lies outside the tunnel ["
                  << upper.at(x) << "," << lower.at(x) << "] at that column";
                return s.str();
            }
        }
    return "";
}

static string build_matrix(const vector<int> &upper, const vector<int> &lower, int length2)
{
    vector<int> starts(upper), ends(lower);
    try
    {
        Tunnel_matrix<int> m((int)starts.size(), length2 + 1, starts, ends);
    }
    catch(const std::exception &e)
    {
        return string("Tunnel_matrix threw ") + e.what();
    }
    return "";
}

int main()
{
    // Defaults (anchors-offset 15, exonerate-hit-trim 5) as main() sets them.
    const char *fake_argv[] = {"test_tunnel_is_monotonic"};
    Settings_handle::st.read_command_line_arguments(1, const_cast<char**>(fake_argv));

    string s1(150, 'A');
    string s2(250, 'A');

    vector<Substring_hit> crossing;
    crossing.push_back(hit(0, 100, 30, 30));    // A
    crossing.push_back(hit(40, 150, 30, 30));   // B
    crossing.push_back(hit(80, 40, 40, 40));    // C

    // [order]
    vector<Substring_hit> kept(crossing);
    Find_anchors fa;
    fa.check_hits_order_conflict(&s1, &s2, &kept);

    bool colinear = true;
    stringstream kept_s;
    for(size_t i = 0; i < kept.size(); i++)
    {
        kept_s << (i ? ", " : "") << kept[i].start_site_1 << "->" << kept[i].start_site_2;
        if(i > 0 && kept[i].start_site_2 < kept[i-1].start_site_2)
            colinear = false;
    }
    report(colinear, "[order]   check_hits_order_conflict() keeps only co-linear hits",
           "kept (start_site_1->start_site_2): " + kept_s.str());

    // [tunnel] and [anchors]
    vector<int> upper, lower;
    fa.define_tunnel(&kept, &upper, &lower, &s1, &s2);
    int d = first_descent(upper, lower);
    string threw = build_matrix(upper, lower, (int)s2.length());
    stringstream tunnel_s;
    if(d >= 0)
        tunnel_s << "bounds descend at column " << d << " (start " << upper[d-1] << "->" << upper[d]
                 << ", end " << lower[d-1] << "->" << lower[d] << "); ";
    tunnel_s << (threw.empty() ? "Tunnel_matrix built" : threw);
    report(d < 0 && threw.empty(), "[tunnel]  the tunnel for those hits is monotonic and buildable",
           tunnel_s.str());

    string outside = anchor_outside(kept, upper, lower);
    report(outside.empty(), "[anchors] every anchor cell is inside that tunnel", outside);

    // [guard] crossing hits straight into define_tunnel(). The trim
    // check_hits_order_conflict() would have applied is irrelevant here.
    vector<Substring_hit> raw(crossing);
    vector<int> gupper, glower;
    fa.define_tunnel(&raw, &gupper, &glower, &s1, &s2);
    d = first_descent(gupper, glower);
    threw = build_matrix(gupper, glower, (int)s2.length());
    // C overwrites nothing of A or B on the diagonal, so all three hits'
    // cells are anchors.
    outside = anchor_outside(raw, gupper, glower);
    stringstream guard_s;
    if(d >= 0)
        guard_s << "bounds descend at column " << d << "; ";
    guard_s << (threw.empty() ? "Tunnel_matrix built" : threw);
    if(!outside.empty())
        guard_s << "; " << outside;
    report(d < 0 && threw.empty() && outside.empty(),
           "[guard]   define_tunnel() given crossing hits widens to a monotonic tunnel holding every anchor",
           guard_s.str());

    return failures ? 1 : 0;
}
