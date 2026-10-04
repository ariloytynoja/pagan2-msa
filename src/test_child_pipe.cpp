// Synthetic test for utils/child_pipe.h: helpers are run WITHOUT A SHELL.
//
// What is pinned, each a way the replacement for popen(3) could look right
// and not be:
//
// 1. NO SHELL RE-PARSES THE ARGUMENTS. An argument holding a space, a `;`
//    and a `$VAR` reaches the program as ONE literal argument; under popen it
//    would have been split, and the `;` would have run a second command.
// 2. STDERR IS IN THE STREAM, as the old `2>&1` put it there -- exonerate's
//    error text must still reach pagan2's reader.
// 3. THE EXIT STATUS SURVIVES, readable with WEXITSTATUS as pclose's was.
// 4. A MISSING PROGRAM IS REPORTED BY OPEN FAILING, not by a shell printing
//    "not found" into the output stream and exiting 127 later.
// 5. A bare name is found on PATH, like the shell found it; a name with a
//    '/' is run as given.
// 6. The reader sees EOF: the parent must not keep the pipe's write end.
// 7. ...nor may a helper spawned by ANOTHER THREAD while this pipe exists.
//    pagan2 runs exonerate from inside OpenMP-parallel alignment; a pipe that
//    is not close-on-exec is inherited by every helper spawned meanwhile, and
//    this reader then waits for EOF until that unrelated helper exits.
// 8. A parent running with stdin and stdout CLOSED still gets the child's
//    stdout. pipe() then returns descriptors 0 and 1, and the child's
//    "dup2(w,1); close(w)" would close the very stdout it had just set up.
//
// Build (from src/):
//   g++ -std=c++11 -w -I. -Iutils -o test_child_pipe test_child_pipe.cpp
// Run: ./test_child_pipe  (exit 0 and PASS lines, or FAIL + 1)
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <vector>

#include "utils/child_pipe.h"

static int failures = 0;

static void check(bool ok, const std::string &what)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
    if (!ok)
        failures++;
}

static std::string slurp(const std::vector<std::string> &argv, int *status,
                         bool *opened)
{
    Child_pipe cp;
    *opened = child_pipe_open(argv, &cp);
    std::string out;
    if (!*opened)
        return out;
    char buf[256];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, cp.fp)) > 0)
        out.append(buf, n);
    *status = child_pipe_close(&cp);
    return out;
}

int main()
{
    int st = -1;
    bool opened = false;

    std::vector<std::string> a;
    a.push_back("printf");
    a.push_back("[%s]\\n");
    a.push_back("a b; echo INJECTED $HOME");
    std::string out = slurp(a, &st, &opened);
    check(opened && out == "[a b; echo INJECTED $HOME]\n",
          "an argument with spaces, ';' and '$' arrives as one literal argument");
    check(out.find("INJECTED\n") == std::string::npos,
          "no second command was run");

    std::vector<std::string> b;
    b.push_back("/bin/sh");
    b.push_back("-c");
    b.push_back("echo to-stdout; echo to-stderr >&2; exit 3");
    out = slurp(b, &st, &opened);
    check(opened && out.find("to-stdout") != std::string::npos
                 && out.find("to-stderr") != std::string::npos,
          "stderr is merged into the stream, as 2>&1 did");
    check(WIFEXITED(st) && WEXITSTATUS(st) == 3,
          "the exit status is returned like pclose's");

    std::vector<std::string> c;
    c.push_back("no_such_program_qzx7");
    out = slurp(c, &st, &opened);
    check(!opened, "a missing program fails at open, not as output text");

    std::vector<std::string> d;
    d.push_back("/bin/echo");
    d.push_back("absolute");
    out = slurp(d, &st, &opened);
    check(opened && out == "absolute\n" && WEXITSTATUS(st) == 0,
          "a path with '/' is executed as given");

    std::vector<std::string> e;
    child_pipe_append_flags(&e, "  -T dna\t-Q  dna ");
    check(e.size() == 4 && e[0] == "-T" && e[1] == "dna" && e[2] == "-Q"
                        && e[3] == "dna",
          "literal flag text splits into separate arguments");

    // EOF must arrive: `true` prints nothing and exits; if the parent kept
    // the write end the read below would block for ever.
    std::vector<std::string> f;
    f.push_back("true");
    out = slurp(f, &st, &opened);
    check(opened && out.empty() && WEXITSTATUS(st) == 0,
          "the reader sees EOF when the child exits");

    // 7. Deterministic stand-in for the thread race: hold one pipe open, as
    // a thread between child_pipe_make() and its own spawn does, start a
    // 3-second helper meanwhile, then close our write end. EOF must come at
    // once; a write end the helper inherited delays it until the helper exits.
    int held[2];
    bool made = child_pipe_make(held);
    std::vector<std::string> g;
    g.push_back("sleep");
    g.push_back("3");
    Child_pipe slow;
    bool slow_opened = child_pipe_open(g, &slow);
    struct timeval t0, t1;
    gettimeofday(&t0, NULL);
    double waited = -1;
    if (made)
    {
        close(held[1]);
        char ch;
        while (read(held[0], &ch, 1) > 0) {}
        gettimeofday(&t1, NULL);
        waited = (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1e6;
        close(held[0]);
    }
    if (slow_opened)
        child_pipe_close(&slow);
    check(made && slow_opened && waited >= 0 && waited < 1.0,
          "a helper spawned meanwhile does not inherit another pipe's write end");

    // 8. In a forked copy of this test with stdin and stdout closed.
    pid_t pid = fork();
    if (pid == 0)
    {
        close(0);
        close(1);
        std::vector<std::string> h;
        h.push_back("/bin/sh");
        h.push_back("-c");
        h.push_back("echo to-stdout; echo to-stderr >&2");
        int hst = -1;
        bool hopened = false;
        std::string hout = slurp(h, &hst, &hopened);
        _exit(hopened && hout.find("to-stdout") != std::string::npos
                      && hout.find("to-stderr") != std::string::npos ? 0 : 1);
    }
    int fst = -1;
    while (pid > 0 && waitpid(pid, &fst, 0) < 0 && errno == EINTR) {}
    check(pid > 0 && WIFEXITED(fst) && WEXITSTATUS(fst) == 0,
          "with stdin and stdout closed, the child's stdout still reaches the reader");

    return failures ? 1 : 0;
}
