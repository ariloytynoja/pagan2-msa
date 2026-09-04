#ifndef FASTTREE_TREE_H
#define FASTTREE_TREE_H

#include "utils/settings_handle.h"
#include "utils/temp_dir.h"
#include "utils/fasta_entry.h"
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace std;

namespace ppa {

class FastTree_tree
{
    string progpath;

    std::string get_temp_dir()
    {
        return resolve_temp_dir();
    }

    void delete_files(int r);

public:
    FastTree_tree();
    bool test_executable();
    string infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein, int n_threads);
};
}

#endif // FASTTREE_TREE_H
