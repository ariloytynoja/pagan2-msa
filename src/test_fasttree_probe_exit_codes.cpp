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
//    The `</dev/null` in the probe command is load-bearing.  The second half
//    of this test demonstrates both directions with a bounded `timeout` so
//    that the failing case cannot hang the test run itself.
//
// Build (from src/):
//   g++ -std=c++11 -w -I. -Iutils -Imain -o test_fasttree_probe_exit_codes \
//       test_fasttree_probe_exit_codes.cpp
// Run: ./test_fasttree_probe_exit_codes  (exit 0 and PASS lines, or FAIL + 1)
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include "utils/helper_probe.h"

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

// Is `timeout` (coreutils) available? The hang half of the test needs it to
// bound the failing case; without it we skip rather than risk hanging.
static bool have_timeout()
{
    return helper_was_found(system("timeout --help >/dev/null 2>/dev/null"))
           && WEXITSTATUS(system("timeout --help >/dev/null 2>/dev/null")) == 0;
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

    // ---- 2. the redirect is what stops the probe hanging ----------------

    if(!have_timeout())
        cout << "SKIP: `timeout` unavailable; not testing the stdin hang"
             << endl;
    else
    {
        // `cat` stands in for FastTree here: both read standard input when
        // given no file. The parent side of the pipe stays open for longer
        // than the timeout, which is what an interactive shell or a job
        // launcher looks like to the child.
        //
        // Without the redirect the probe hangs and `timeout` kills it (124).
        int hung = system("sleep 10 | timeout 2 cat >/dev/null 2>/dev/null");
        check(WIFEXITED(hung) && WEXITSTATUS(hung) == 124,
              "a stdin-reading helper with NO redirect hangs (timeout kills it)");

        // With the redirect it returns immediately, whatever the parent's
        // stdin is doing.
        int fine = system(
            "sleep 10 | timeout 2 cat </dev/null >/dev/null 2>/dev/null");
        check(WIFEXITED(fine) && WEXITSTATUS(fine) == 0,
              "the same helper with `</dev/null` returns at once");
    }

    if(failures)
    {
        cout << failures << " failure(s)" << endl;
        return 1;
    }
    cout << "all checks passed" << endl;
    return 0;
}
