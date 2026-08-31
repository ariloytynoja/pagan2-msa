// Synthetic test for the object lifetime that made
// Reads_aligner::query_placement_all() read and write freed memory.
//
// In the reverse-strand branch the code kept node_rc and discarded node:
//
//     current_root = node_rc;
//     node->has_left_child(false);
//     delete node;
//     ...
//     n<<node->right_child->get_name()<<"."<<it->second;   // freed
//     node->right_child->set_name(n.str());                // freed, and a WRITE
//
// has_left_child(false) spares only the LEFT child -- the shared reference
// subtree that must survive. ~Node() still deletes right_child, which is the
// read node itself, so both lines above touch memory that `delete node` just
// released. The forward branch a few lines earlier is correct because there
// the surviving node IS node; only the reverse branch got it backwards. The
// same situation in query_placement_one() was already written correctly, using
// node_rc->right_child.
//
// This test asserts that asymmetry directly, without invoking any undefined
// behaviour: it watches the two child addresses through a global operator
// delete and checks which of them ~Node() releases.
//
// Build (from src/, after a normal build has produced the .o files), as one
// command line -- the object list is everything node.o pulls in:
//
//   g++ -std=c++11 -w -I. -Iutils -Imain
//       -o test_node_delete_frees_right_child
//       test_node_delete_frees_right_child.cpp
//       node.o sequence.o viterbi_alignment.o basic_alignment.o
//       reference_alignment.o model_factory.o evol_model.o db_matrix.o
//       int_matrix.o eigen.o settings.o settings_handle.o log_output.o
//       text_utils.o codon_translation.o find_anchors.o exonerate_queries.o
//       check_version.o tunnel_matrix.o fasta_reader.o xml_writer.o
//       -lboost_program_options -lboost_regex -lboost_thread -lboost_system
//       -lgomp -lm -lz -lcurl -lpthread -ldl
//
// Run: ./test_node_delete_frees_right_child   (exit 0 and PASS lines, or FAIL)
//
// End-to-end reproducer for the bug itself, for a build with
// -fsanitize=address (ref.fas holds two ~120 nt sequences differing in the
// last base, q.fas holds the reverse complement of refA TWICE under one name
// -- the duplicate name is what reaches the renaming branch):
//
//   pagan2 -a ref.fas -t ref.tre -q q.fas -o out --fragments --both-strands
//
// Before the fix AddressSanitizer reports
//   heap-use-after-free ... READ of size 8
//     #0 ppa::Reads_aligner::query_placement_all(...) main/reads_aligner.cpp:573
//   freed by thread T0 here:
//     #1 ppa::Reads_aligner::query_placement_all(...) main/reads_aligner.cpp:566
// After the fix the same command completes cleanly.

#include <iostream>
#include <string>
#include <cstdlib>
#include <new>
#include "main/node.h"
#include "utils/settings_handle.h"

using namespace ppa;
using namespace std;

// Two addresses to watch. Recording is done with plain flags and malloc/free
// so the hooks never allocate and cannot recurse.
static void *g_left  = 0;
static void *g_right = 0;
static bool  g_left_freed  = false;
static bool  g_right_freed = false;

// noexcept rather than the C++03 dynamic exception specification: the latter
// is gone in C++17, which this project has already had to work around once.
void *operator new(std::size_t size)
{
    void *p = std::malloc(size ? size : 1);
    if(!p) throw std::bad_alloc();
    return p;
}

void operator delete(void *p) noexcept
{
    if(p != 0)
    {
        if(p == g_left)  g_left_freed  = true;
        if(p == g_right) g_right_freed = true;
        std::free(p);
    }
}

void *operator new[](std::size_t size) { return operator new(size); }
void operator delete[](void *p) noexcept { operator delete(p); }

#ifdef __cpp_sized_deallocation
// Only emitted from C++14 on; the project builds as C++11, where these are
// never called, but defining them keeps the test correct if it is built with
// a newer -std.
void operator delete(void *p, std::size_t) noexcept { operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { operator delete(p); }
#endif

static int failures = 0;

static void check(bool ok, const string &what)
{
    cout << (ok ? "PASS: " : "FAIL: ") << what << endl;
    if(!ok) failures++;
}

int main()
{
    const char *fake_argv[] = {"test_node_delete_frees_right_child"};
    Settings_handle::st.read_command_line_arguments(1, const_cast<char**>(fake_argv));

    // The shape query_placement_all builds: a temporary parent whose left
    // child is the shared reference subtree and whose right child is the read.
    Node *parent = new Node();
    Node *reference_subtree = new Node();
    Node *read_node = new Node();

    reference_subtree->set_name("reference_subtree");
    read_node->set_name("read");

    parent->add_left_child(reference_subtree);
    parent->add_right_child(read_node);

    g_left  = static_cast<void*>(reference_subtree);
    g_right = static_cast<void*>(read_node);

    // Exactly what the discard path does before deleting the losing node.
    parent->has_left_child(false);
    delete parent;

    check(!g_left_freed,
          "has_left_child(false) spares the left child (shared reference subtree)");
    check(g_right_freed,
          "~Node() still deletes the right child, so node->right_child is "
          "dangling after 'delete node'");

    // The left child was deliberately not freed by ~Node(); release it here so
    // the test itself leaks nothing.
    g_left = 0;
    delete reference_subtree;

    if(failures == 0)
        cout << "All checks passed." << endl;
    else
        cout << failures << " check(s) failed." << endl;

    return failures == 0 ? 0 : 1;
}
