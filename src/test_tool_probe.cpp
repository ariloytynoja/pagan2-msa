// Synthetic test for utils/tool_probe.h: find a helper the way the shell's
// command search would, without running it.
//
// What is pinned:
// 1. A DIRECTORY with the tool's name is not the tool. access(X_OK) alone
//    accepts it (a directory's x bit means "searchable"), and a directory
//    called `exonerate` beside pagan2 or on $PATH would then be reported as
//    exonerate and fail only when the alignment tries to run it.
// 2. A regular executable file is found, beside the binary (the directory is
//    returned, with its '/') and on $PATH (returned as "", i.e. by bare name).
// 3. With $PATH UNSET the system default path is searched, as the shell
//    does -- `sh` is found -- instead of nothing.
// 4. Concurrent first lookups are safe: the probe is called from
//    OpenMP-parallel alignment and memoises into a shared map. Build the
//    test with -fsanitize=thread as well to have the race checked rather
//    than merely survived.
//
// Build (from src/):
//   g++ -std=c++11 -w -I. -Iutils -pthread -o test_tool_probe test_tool_probe.cpp
//   (race check: add -fsanitize=thread -g)
// Run: ./test_tool_probe  (exit 0 and PASS lines, or FAIL + 1)
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

#include "utils/tool_probe.h"

static int failures = 0;

static void check(bool ok, const std::string &what)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
    if (!ok)
        failures++;
}

int main()
{
    const char *tmp = getenv("TMPDIR");
    std::string base = std::string(tmp && *tmp ? tmp : ".") + "/test_tool_probe.XXXXXX";
    std::vector<char> tmpl(base.begin(), base.end());
    tmpl.push_back('\0');
    if (mkdtemp(&tmpl[0]) == NULL)
    {
        std::cout << "FAIL cannot create a scratch directory under " << base << std::endl;
        return 1;
    }
    std::string root(&tmpl[0]);
    std::string with_dir = root + "/with_dir/";
    std::string with_file = root + "/with_file/";
    mkdir(with_dir.c_str(), 0700);
    mkdir(with_file.c_str(), 0700);
    mkdir((with_dir + "fake_tool_qzx7").c_str(), 0700);           // a directory
    {
        std::ofstream f((with_file + "fake_tool_qzx7").c_str());
        f << "#!/bin/sh\nexit 0\n";
    }
    chmod((with_file + "fake_tool_qzx7").c_str(), 0700);           // the tool

    std::string where = "unset";
    check(!ppa::find_tool(with_dir, "fake_tool_qzx7", &where),
          "a directory beside the binary is not the tool");

    where = "unset";
    check(ppa::find_tool(with_file, "fake_tool_qzx7", &where) && where == with_file,
          "an executable file beside the binary is found, with its directory");

    std::string saved_path = getenv("PATH") ? getenv("PATH") : "";

    setenv("PATH", (root + "/with_dir").c_str(), 1);
    check(!ppa::find_tool("", "fake_tool_qzx7", &where),
          "a directory on $PATH is not the tool");

    setenv("PATH", (root + "/with_dir:" + root + "/with_file").c_str(), 1);
    where = "unset";
    check(ppa::find_tool("/nonexistent_dir_qzx7/", "fake_tool_qzx7", &where) && where.empty(),
          "the executable file later on $PATH is found, to be run by bare name");

    unsetenv("PATH");
    check(ppa::find_tool("", "sh", &where) && where.empty(),
          "with $PATH unset the system default path is searched");

    if (!saved_path.empty())
        setenv("PATH", saved_path.c_str(), 1);

    // Concurrent first lookups of many distinct tools: every insert is a
    // cache miss, which is the case that writes to the shared map.
    const int n_threads = 8, per_thread = 200;
    std::vector<int> wrong(n_threads, 0);
    std::vector<std::thread> threads;
    for (int t = 0; t < n_threads; t++)
        threads.push_back(std::thread([t, &wrong, &with_file, per_thread]() {
            for (int i = 0; i < per_thread; i++)
            {
                char name[64];
                snprintf(name, sizeof name, "no_such_tool_%d_%d", t, i);
                std::string w;
                if (ppa::find_tool(with_file, name, &w))
                    wrong[t]++;
                if (!ppa::find_tool(with_file, "fake_tool_qzx7", &w) || w != with_file)
                    wrong[t]++;
            }
        }));
    for (size_t t = 0; t < threads.size(); t++)
        threads[t].join();
    int total_wrong = 0;
    for (int t = 0; t < n_threads; t++)
        total_wrong += wrong[t];
    check(total_wrong == 0, "concurrent lookups from 8 threads all answer correctly");

    rmdir((with_dir + "fake_tool_qzx7").c_str());
    remove((with_file + "fake_tool_qzx7").c_str());
    rmdir(with_dir.c_str());
    rmdir(with_file.c_str());
    rmdir(root.c_str());
    return failures ? 1 : 0;
}
