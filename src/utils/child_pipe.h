/***************************************************************************
 *   Copyright (C) 2026 by Martin Mokrejs                                  *
 *   GPL-3.0-or-later, as the rest of PAGAN2                               *
 ***************************************************************************/

#ifndef CHILD_PIPE_H
#define CHILD_PIPE_H

/*
 * Run a helper program and read what it prints, WITHOUT A SHELL.
 *
 * popen(3) runs `/bin/sh -c "<command>"`, so every helper invocation costs a
 * shell process before the helper even starts -- and the exonerate
 * preselection runs once PER QUERY. On a busy host with ~150 concurrent
 * pagan2 processes that shell was measured as the step blocked in
 * `rpc_wait_bit_killable`, i.e. waiting on the network filesystem, while
 * pagan2 itself sat in read() on the pipe; the whole alignment rate was
 * then gated by the shell's start-up, not by exonerate.
 *
 * child_pipe_open() instead posix_spawn()s the program with an explicit
 * argument vector: one process per call, no quoting, and a path containing a
 * space or a shell metacharacter is passed as one argument, not re-parsed.
 * The child's stdout AND stderr go to the pipe, which is what the old
 * command's trailing `2>&1` asked for, so the reader sees the same stream.
 *
 * A program name without a '/' is looked up on PATH (posix_spawnp), exactly
 * as the shell did; a name with one is executed as given.
 *
 * The pipe is created close-on-exec, as popen(3)'s is: the helpers are run
 * from inside OpenMP-parallel alignment, and a write end inherited by a
 * helper another thread spawns at the same moment keeps this reader from
 * seeing EOF until that unrelated helper exits.
 *
 * child_pipe_close() returns the child's wait status like pclose(3) does
 * (use WIFEXITED/WEXITSTATUS on it), or -1 if it could not be waited for.
 */

#include <cerrno>
#include <cstdio>
#include <fcntl.h>
#include <spawn.h>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

extern char **environ;

struct Child_pipe
{
    FILE *fp;
    pid_t pid;
    Child_pipe() : fp(NULL), pid(-1) {}
};

// Split a whitespace-separated list of literal flags ("-T dna -Q dna") into
// arguments. Only for constant flag text written in this program, never for
// file names or anything user supplied.
inline void child_pipe_append_flags(std::vector<std::string> *argv,
                                    const std::string &flags)
{
    std::string tok;
    for (size_t i = 0; i <= flags.size(); i++)
    {
        if (i == flags.size() || flags[i] == ' ' || flags[i] == '\t')
        {
            if (!tok.empty())
                argv->push_back(tok);
            tok.clear();
        }
        else
            tok += flags[i];
    }
}

// A close-on-exec pipe whose two ends are both above stderr. If this process
// runs with stdin or stdout closed, pipe() hands out 0 or 1, and the write
// end landing on 1 would make the child's "dup2(w,1); close(w)" close its own
// stdout. Where pipe2() is missing, FD_CLOEXEC is set right after pipe(),
// which leaves a short window for a concurrent spawn.
inline bool child_pipe_make(int fds[2])
{
#if defined(__linux__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
    if (pipe2(fds, O_CLOEXEC) != 0)
        return false;
#else
    if (pipe(fds) != 0)
        return false;
    fcntl(fds[0], F_SETFD, FD_CLOEXEC);
    fcntl(fds[1], F_SETFD, FD_CLOEXEC);
#endif
    for (int i = 0; i < 2; i++)
    {
        if (fds[i] > 2)
            continue;
        int moved = fcntl(fds[i], F_DUPFD_CLOEXEC, 3);
        if (moved < 0)
        {
            int e = errno;
            close(fds[0]);
            close(fds[1]);
            errno = e;
            return false;
        }
        close(fds[i]);
        fds[i] = moved;
    }
    return true;
}

inline bool child_pipe_open(const std::vector<std::string> &argv,
                            Child_pipe *out)
{
    if (argv.empty() || argv[0].empty())
        return false;

    int fds[2];
    if (!child_pipe_make(fds))
        return false;

    posix_spawn_file_actions_t fa;
    if (posix_spawn_file_actions_init(&fa) != 0)
    {
        close(fds[0]);
        close(fds[1]);
        return false;
    }
    // Child: the pipe's write end becomes stdout and stderr ("2>&1"); both
    // original pipe descriptors are closed so the reader sees EOF when the
    // helper exits.
    posix_spawn_file_actions_addclose(&fa, fds[0]);
    posix_spawn_file_actions_adddup2(&fa, fds[1], 1);
    posix_spawn_file_actions_adddup2(&fa, fds[1], 2);
    posix_spawn_file_actions_addclose(&fa, fds[1]);

    std::vector<char *> args;
    for (size_t i = 0; i < argv.size(); i++)
        args.push_back(const_cast<char *>(argv[i].c_str()));
    args.push_back(NULL);

    pid_t pid = -1;
    int rc;
    if (argv[0].find('/') != std::string::npos)
        rc = posix_spawn(&pid, args[0], &fa, NULL, &args[0], environ);
    else
        rc = posix_spawnp(&pid, args[0], &fa, NULL, &args[0], environ);
    posix_spawn_file_actions_destroy(&fa);
    close(fds[1]);

    if (rc != 0)
    {
        close(fds[0]);
        errno = rc;
        return false;
    }

    FILE *fp = fdopen(fds[0], "r");
    if (fp == NULL)
    {
        close(fds[0]);
        int st;
        while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {}
        return false;
    }
    out->fp = fp;
    out->pid = pid;
    return true;
}

inline int child_pipe_close(Child_pipe *cp)
{
    if (cp->fp != NULL)
        fclose(cp->fp);
    cp->fp = NULL;
    if (cp->pid <= 0)
        return -1;
    int st = 0;
    pid_t r;
    while ((r = waitpid(cp->pid, &st, 0)) < 0 && errno == EINTR) {}
    cp->pid = -1;
    return r < 0 ? -1 : st;
}

#endif // CHILD_PIPE_H
