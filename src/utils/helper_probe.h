#ifndef HELPER_PROBE_H
#define HELPER_PROBE_H

#include <sys/wait.h>

namespace ppa {

// Did the helper program named in a `system()` probe actually run?
//
// Every test_executable() in pagan2 asks this question, and each one answers
// it by comparing against the exit status the helper happens to produce when
// invoked with no useful arguments -- 0 for one program, 1 for the next. That
// couples a presence check to a detail of the helper's own CLI, and it is
// wrong whenever the helper changes that detail or was never surveyed on the
// platform in hand.
//
// The status that does NOT depend on the helper is the shell's. `system()`
// runs the command through /bin/sh, and POSIX fixes two codes for it:
//
//     127  command not found
//     126  found, but could not be executed (not executable, bad interpreter)
//
// Anything else -- including a non-zero exit the helper chose itself, and
// including death by signal -- means the program was there and ran. That is
// all a presence probe needs to establish; whether it liked its arguments is
// a different question, and not one the probe asked.
//
// `system()` returns -1 when the fork or wait failed, in which case nothing
// ran at all.
inline bool helper_was_found(int status)
{
    if(status == -1)
        return false;
    if(!WIFEXITED(status))
        return true;
    int code = WEXITSTATUS(status);
    return code != 127 && code != 126;
}

}

#endif // HELPER_PROBE_H
