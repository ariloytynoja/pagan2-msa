// Synthetic test for Reads_aligner::find_target_node().
//
// Reads_aligner looked up a target node name with nodes.find(name) at five
// separate sites and used the result as nit->second without ever comparing
// it to nodes.end().  When the name is absent that is undefined behaviour:
// it does not reliably fault at the lookup, it yields a garbage Node* that
// is then passed to read_match_score()/query_match_score() and crashes
// somewhere else entirely -- observed as a SIGSEGV inside
// Node::get_distance_to_parent() with this=0x13, several frames away from
// the map lookup that actually caused it.
//
// find_target_node() replaces that pattern with a lookup that returns NULL
// for an absent name.  This test builds a small synthetic node map and
// checks both directions.
//
// Build (from src/):
//   g++ -std=c++11 -w -I. -Iutils -Imain -o test_find_target_node \
//       test_find_target_node.cpp
// Run: ./test_find_target_node   (exit 0 and PASS lines, or FAIL + exit 1)
#include <iostream>
#include <map>
#include <string>
#include <cstdlib>
#include "main/reads_aligner.h"

using namespace ppa;
using namespace std;

static int failures = 0;

static void check(bool ok, const string &what)
{
    if(ok)
    {
        cout << "PASS: " << what << endl;
    }
    else
    {
        cout << "FAIL: " << what << endl;
        failures++;
    }
}

int main()
{
    // A synthetic node map standing in for one built by
    // Node::get_all_nodes(): two named nodes, nothing else.
    //
    // Deliberately heap-allocated and never deleted.  Node's destructor
    // lives in node.cpp and tears down a whole subtree, which would drag
    // most of PAGAN2's object graph into this test's link line; the
    // constructor alone is header-inline.  Leaking two nodes from a test
    // process that is about to exit costs nothing and keeps the test
    // buildable from a single translation unit.
    Node *a = new Node();
    a->set_name("#1#");
    Node *b = new Node();
    b->set_name("querySeqA");

    map<string,Node*> nodes;
    nodes.insert(make_pair(a->get_name(),a));
    nodes.insert(make_pair(b->get_name(),b));

    // Present names must resolve to exactly the node that was inserted.
    check(Reads_aligner::find_target_node(nodes,"#1#") == a,
          "an internal node name present in the tree resolves to its node");
    check(Reads_aligner::find_target_node(nodes,"querySeqA") == b,
          "a leaf name present in the tree resolves to its node");

    // The regression: a name that is NOT in the tree must come back NULL
    // rather than a dereferenced end() iterator.  Any of these is what an
    // aligner hit can name when its target index and the tree disagree.
    check(Reads_aligner::find_target_node(nodes,"#2#") == 0,
          "an absent internal node name returns NULL, not garbage");
    check(Reads_aligner::find_target_node(nodes,"querySeqB") == 0,
          "an absent leaf name returns NULL, not garbage");
    check(Reads_aligner::find_target_node(nodes,"") == 0,
          "an empty name returns NULL, not garbage");

    // An empty map is the degenerate case of the same thing.
    map<string,Node*> empty_nodes;
    check(Reads_aligner::find_target_node(empty_nodes,"#1#") == 0,
          "any lookup in an empty node map returns NULL");

    if(failures > 0)
    {
        cout << failures << " check(s) failed." << endl;
        return 1;
    }
    cout << "All checks passed." << endl;
    return 0;
}
