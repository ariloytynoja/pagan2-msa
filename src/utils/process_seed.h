#ifndef PROCESS_SEED_H
#define PROCESS_SEED_H

#include <stdint.h>

namespace ppa {

// The seed main() gives srand(), from the clock and the process id.
//
// rand() names the scratch files every helper wrapper writes ("q<r>.fas",
// "t<r>.fas", "m<r>.fas", ...) in the one directory get_temp_dir() returns,
// and each of those sites relies on a check-then-create retry loop that
// assumes two processes walk DIFFERENT candidate names. Two processes that
// share a seed walk the same names and use each other's files.
//
// time(0) alone is the same for every process started in one second. Mixing
// in the pid with XOR is not enough either: t ^ pid == (t+1) ^ (pid+1)
// whenever t and pid are both even, which is exactly a launcher starting one
// pagan2 a second with consecutive pids. Both values are therefore put
// through a 64-bit mixer (the splitmix64 finaliser), so that seeds collide
// only by chance, not along any pattern of nearby (time, pid) pairs.
inline unsigned process_seed(uint64_t t, uint64_t pid)
{
    uint64_t z = t * 0x9E3779B97F4A7C15ULL ^ pid;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return static_cast<unsigned>(z ^ (z >> 32));
}

}

#endif // PROCESS_SEED_H
