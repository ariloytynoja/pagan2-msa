// Synthetic test for the helper-presence probe used by FastTree_tree.
//
// Two separate defects are pinned here, both of which made
// FastTree_tree::test_executable() report "FastTree is not available" on
// machines where FastTree is installed and works.
//
// 1. WHICH EXIT STATUS MEANS "PRESENT".  The probe used to require exit 0
//    from `fasttree` invoked with no arguments.  FastTree given no input file
//    reads the alignment from standard input, so it never exits 0 there: it
//    reads EOF, has nothing to build a tree from, and exits 1.  The old code
//    compensated with a hard-coded `gethostname() == "wasabi2"` case that
//    accepted exit 1 on one host.  helper_was_found() instead keys on the
//    SHELL's convention, which does not depend on the helper at all:
//    /bin/sh returns 127 for "not found" and 126 for "not executable".
//
// 2. WHETHER THE PROBE CAN HANG.  Because FastTree reads standard input, a
//    probe that does not redirect it BLOCKS FOREVER whenever the parent's
//    own stdin is open -- an interactive shell, a pipeline, a job launcher.
//    The `</dev/null` in the probe command is load-bearing.
//
// The second half runs the REAL FastTree_tree::test_executable() and
// infer_phylogeny(), with --fasttree-exec naming a stand-in script in a
// directory whose name holds a space. The stand-in records what its standard
// input is and, when given arguments, prints a tree. This test's own stdin is
// made a pipe that stays open, as a launcher's would be. Checked: the probe
// finds the program by that path (it is ONE shell word), its stdin is
// /dev/null rather than the inherited pipe, and the tree comes back.
//
// Build (from src/, after a normal build has produced the objects; every
// object but main.o):
//   g++ -std=c++11 -w -fopenmp -I. -Iutils -Imain \
//       -o test_fasttree_probe_exit_codes test_fasttree_probe_exit_codes.cpp \
//       $(ls *.o | grep -v '^main.o$') -lboost_program_options -lboost_regex \
//       -lboost_thread -lboost_system -lgomp -lm -lz -lpthread -ldl
// Run: ./test_fasttree_probe_exit_codes  (exit 0 and PASS lines, or FAIL + 1)
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include "utils/helper_probe.h"
#include "utils/fasttree_tree.h"
#include "utils/settings_handle.h"

using namespace ppa;
using namespace std;

static int failures = 0;

static void check(bool ok, const string &what)
{
    if(ok)
        cout << "PASS: " << what << endl;
    else
    {
        cout << "FAIL: " << what << endl;
        failures++;
    }
}

int main()
{
    // ---- 1. which status means "the program ran" ------------------------

    check(helper_was_found(system("true")) == true,
          "exit 0 means the program ran");

    // The case the old probe got wrong: a program that is present, runs, and
    // exits non-zero because it did not like its arguments. FastTree with no
    // input file is exactly this.
    check(helper_was_found(system("false")) == true,
          "exit 1 means the program ran (it just failed)");

    check(helper_was_found(
              system("no_such_program_qzx7 >/dev/null 2>/dev/null")) == false,
          "shell 127 (command not found) means absent");

    // 126: present on disk but not executable. Build one and try to run it.
    {
        const char *path = "./test_probe_not_executable.tmp";
        FILE *f = fopen(path, "w");
        if(f)
        {
            fputs("#!/bin/sh\nexit 0\n", f);
            fclose(f);
            chmod(path, 0644);
            string cmd = string(path) + " >/dev/null 2>/dev/null";
            check(helper_was_found(system(cmd.c_str())) == false,
                  "shell 126 (found but not executable) means absent");
            remove(path);
        }
        else
            cout << "SKIP: could not create the non-executable fixture" << endl;
    }

    // Death by signal is still evidence the program ran.
    check(helper_was_found(system("kill -TERM $$")) == true,
          "killed by a signal still means the program ran");

    // system() reports -1 when the fork or wait failed: nothing ran.
    check(helper_was_found(-1) == false,
          "system() == -1 (fork/wait failed) means absent");

    // ---- 2. the real probe and run, through --fasttree-exec ------------
    {
        const char *tmp = getenv("TMPDIR");
        string base = string(tmp && *tmp ? tmp : ".") + "/test_fasttree_probe.XXXXXX";
        vector<char> tmpl(base.begin(), base.end());
        tmpl.push_back('\0');
        string root = mkdtemp(&tmpl[0]) ? string(&tmpl[0]) : string();
        string dir = root + "/dir with space";
        string exe = dir + "/fake FastTree";
        string record = root + "/stdin.txt";
        if(root.empty() || mkdir(dir.c_str(), 0700) != 0)
        {
            cout << "FAIL: could not create the scratch directory under " << base << endl;
            return 1;
        }
        {
            ofstream f(exe.c_str());
            f << "#!/bin/sh\n"
                 "if [ $# -eq 0 ]; then\n"
                 "  readlink /proc/self/fd/0 > \"$FT_STDIN_RECORD\" 2>/dev/null"
                 " || echo unknown > \"$FT_STDIN_RECORD\"\n"
                 "  exit 1\n"
                 "fi\n"
                 "echo '(a:0.1,b:0.1);'\n";
        }
        chmod(exe.c_str(), 0700);
        setenv("FT_STDIN_RECORD", record.c_str(), 1);

        // Our stdin becomes a pipe whose write end stays open.
        int fds[2];
        if(pipe(fds) == 0)
        {
            dup2(fds[0], 0);
            close(fds[0]);
        }

        const char *fake_argv[] = {"test_fasttree_probe_exit_codes", "--fasttree-exec", exe.c_str()};
        Settings_handle::st.read_command_line_arguments(3, const_cast<char**>(fake_argv));

        FastTree_tree ft;
        check(ft.test_executable(),
              "--fasttree-exec naming a path with a space is found");

        string seen;
        ifstream r(record.c_str());
        getline(r, seen);
        check(seen == "/dev/null",
              "the probe's stdin is /dev/null, not the inherited open pipe (got '" + seen + "')");

        vector<Fasta_entry> seqs(2);
        seqs[0].name = "a"; seqs[0].sequence = "ACGTACGTAC";
        seqs[1].name = "b"; seqs[1].sequence = "ACGTACGTAA";
        string tree = ft.infer_phylogeny(&seqs, false, 1);
        check(tree.find("(a:0.1,b:0.1);") != string::npos,
              "infer_phylogeny() runs the program by that path and reads its tree");

        remove(exe.c_str());
        remove(record.c_str());
        rmdir(dir.c_str());
        rmdir(root.c_str());
    }

    if(failures)
    {
        cout << failures << " failure(s)" << endl;
        return 1;
    }
    cout << "all checks passed" << endl;
    return 0;
}
