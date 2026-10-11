#ifndef TOOL_PROBE_H
#define TOOL_PROBE_H

/*
 * Is a helper program available? Answer by LOOKING FOR IT, not by running it.
 *
 * Exonerate_queries::test_executable() builds
 * `<dir>/exonerate >/dev/null 2>/dev/null`, hands it to system(3), and on
 * failure retries with a bare name. Two shells and up to two real exonerate
 * processes -- and viterbi_alignment.cpp constructs an Exonerate_queries and
 * calls it INSIDE the alignment path, so it is paid PER PAIRWISE ALIGNMENT.
 *
 * Measured with `strace -f -e trace=execve` on one run, 900 nt reference:
 *
 *     queries   execve   exonerate   sh
 *           1       27          11   13
 *           4      147          71   73
 *
 * Of the 4-query run's 71 exonerate processes, 48 were the probe (24 by bare
 * name, 24 by full path) against 23 real alignments -- so two thirds of every
 * exonerate this program starts exist only to ask whether exonerate exists.
 *
 * HOW IT LOOKS. Like the shell's own command search: the first REGULAR file
 * with execute permission, first beside pagan2's binary, then along $PATH
 * (an empty element is the current directory; an unset $PATH is the system
 * default path, as the shell uses). A directory of the same name does not
 * count. access(2) checks the real rather than the effective user, which only
 * differs for a set-uid pagan2.
 *
 * WHAT IS GIVEN UP, STATED PLAINLY. What is no longer established is
 * "...and it exits 1". An exonerate that is present and executable but
 * broken now reports AVAILABLE and fails at the real invocation, where the
 * failure is visible, instead of reporting absent and having pagan2 silently
 * align without anchors. That direction is deliberate: a silent fallback is
 * the harder failure to diagnose.
 *
 * The answer is memoised per (dir, tool): the filesystem is not expected to
 * change under a running alignment. The cache is guarded by a mutex because
 * the probe is called from OpenMP-parallel alignment.
 */

#include <string>
#include <map>
#include <mutex>
#include <utility>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>

namespace ppa {

inline bool tool_is_executable_file(const std::string &path)
{
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode)
           && access(path.c_str(), X_OK) == 0;
}

/*
 * Is <tool> executable in <dir>, or on $PATH?
 *
 * On success *found_dir is set the way the caller's own variable expects: the
 * directory (with its trailing '/') when the tool sits beside pagan2's binary,
 * and "" when it was found on $PATH and must therefore be invoked by bare
 * name. That is exactly what the system(3) probe assigned, so the caller's
 * later command building is unchanged.
 */
inline bool find_tool(const std::string &dir, const std::string &tool,
                      std::string *found_dir)
{
    static std::mutex lock;
    static std::map<std::string, std::pair<bool, std::string> > cache;
    std::lock_guard<std::mutex> guard(lock);

    std::string key = dir;
    key += '\x1f';                       // not a path character
    key += tool;

    std::map<std::string, std::pair<bool, std::string> >::const_iterator it
        = cache.find(key);
    if(it != cache.end())
    {
        if(found_dir != 0)
            *found_dir = it->second.second;
        return it->second.first;
    }

    bool ok = false;
    std::string where = "";

    if(!dir.empty() && tool_is_executable_file(dir + tool))
    {
        ok = true;
        where = dir;
    }

    if(!ok)
    {
        std::string path;
        const char *path_env = getenv("PATH");
        if(path_env != 0)
            path = path_env;
        else
        {
            size_t n = confstr(_CS_PATH, NULL, 0);
            if(n > 0)
            {
                std::string def(n, '\0');
                confstr(_CS_PATH, &def[0], n);
                path = def.c_str();
            }
            else
                path = "/bin:/usr/bin";
        }

        std::string::size_type start = 0;
        while(start <= path.size())
        {
            std::string::size_type end = path.find(':', start);
            if(end == std::string::npos)
                end = path.size();

            std::string element = path.substr(start, end - start);
            // POSIX: an empty $PATH element means the current directory.
            if(element.empty())
                element = ".";
            if(element[element.size() - 1] != '/')
                element += '/';

            if(tool_is_executable_file(element + tool))
            {
                ok = true;
                where = "";          // invoked by bare name, as before
                break;
            }

            if(end == path.size())
                break;
            start = end + 1;
        }
    }

    cache[key] = std::make_pair(ok, where);
    if(found_dir != 0)
        *found_dir = where;
    return ok;
}

}

#endif // TOOL_PROBE_H
