// Synthetic test for Exonerate_queries::read_full_line().
//
// Exonerate's output was read with
//
//     char line[256];
//     while ( fgets( line, sizeof line, fpipe))
//         ... split_sugar_string(string(line), &h) ...
//
// and split_sugar_string()'s regex is anchored on a trailing "\n". fgets()
// stops after 255 characters, so any output line longer than that arrives in
// two pieces: the first has no newline and the second does not start with
// "sugar:", so neither matches and the hit is dropped silently -- no warning,
// no error, just missing anchors.
//
// A sugar line is
//
//     sugar: <query> <qs> <qe> <strand> <target> <ts> <te> <strand> <score>
//
// so the two sequence names sit in a budget of roughly 225 characters. That is
// easily exceeded by ordinary identifiers. Measured on a real run, with a
// 193-character query name against PAGAN's own short internal target names
// (#1#, #2#, ...), Exonerate emitted 89 hits of which 47 landed on lines of
// 255 characters or more: 42 of 89 parsed, 47 silently lost. The same data
// with a short query name parsed 89 of 89.
//
// read_full_line() reads a line of any length. This test feeds it a file whose
// first line is far past the old limit and checks it comes back whole.
//
// Build (from src/, after a normal build has produced the .o files), as one
// command line:
//
//   g++ -std=c++11 -w -I. -Iutils -Imain
//       -o test_exonerate_long_line_parsing
//       test_exonerate_long_line_parsing.cpp exonerate_queries.o
//       settings.o settings_handle.o log_output.o text_utils.o
//       check_version.o node.o sequence.o viterbi_alignment.o
//       basic_alignment.o reference_alignment.o model_factory.o evol_model.o
//       db_matrix.o int_matrix.o eigen.o codon_translation.o find_anchors.o
//       tunnel_matrix.o fasta_reader.o xml_writer.o
//       -lboost_program_options -lboost_regex -lboost_thread -lboost_system
//       -lgomp -lm -lz -lcurl -lpthread -ldl
//
// Run: ./test_exonerate_long_line_parsing   (exit 0 and PASS lines, or FAIL)

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include "utils/exonerate_queries.h"

using namespace ppa;
using namespace std;

static int failures = 0;

static void check(bool ok, const string &what)
{
    cout << (ok ? "PASS: " : "FAIL: ") << what << endl;
    if(!ok) failures++;
}

int main()
{
    const char *path = "test_exonerate_long_line_parsing.input";

    // A realistic over-long sugar line. On the run this was measured from, a
    // 193-character query name against a 148-character target name produced
    // lines of 371 to 379 characters; reproduce that shape here.
    string query_name = "hCoV-19_England_MILK-9E05B3_2020_EPI_ISL_601443_2020-09-20";
    while(query_name.size() < 193)
        query_name += "x";

    string target_name = "SARS-CoV-1_Urbani_2003_consensus_sample_replicate";
    while(target_name.size() < 148)
        target_name += "y";

    string long_line = "sugar: " + query_name + " 0 3822 + " + target_name
                     + " 0 3822 + 19110\n";
    string short_line = "sugar: q 0 100 + #2# 0 100 + 500\n";
    string last_line  = "sugar: r 0 50 + #3# 0 50 + 250";   // no trailing newline

    {
        ofstream out(path);
        out << long_line << short_line << last_line;
    }

    check(long_line.size() > 256,
          "the synthetic sugar line is longer than the old 256-byte buffer");

    FILE *f = fopen(path, "r");
    if(f == NULL)
    {
        cout << "FAIL: could not open " << path << endl;
        return 1;
    }

    string got;

    check(Exonerate_queries::read_full_line(f, &got), "first line is read");
    check(got == long_line,
          "the over-long line comes back whole, newline included");

    check(Exonerate_queries::read_full_line(f, &got), "second line is read");
    check(got == short_line, "a short line is unchanged");

    check(Exonerate_queries::read_full_line(f, &got), "final line is read");
    check(got == last_line,
          "a final line without a trailing newline is still returned");

    check(!Exonerate_queries::read_full_line(f, &got),
          "end of input reports false");

    fclose(f);
    remove(path);

    if(failures == 0)
        cout << "All checks passed." << endl;
    else
        cout << failures << " check(s) failed." << endl;

    return failures == 0 ? 0 : 1;
}
