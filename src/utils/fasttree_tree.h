#ifndef FASTTREE_TREE_H
#define FASTTREE_TREE_H

#include "utils/settings_handle.h"
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
        std::string tmp_dir = "/tmp/";

        if(Settings_handle::st.is("temp-folder"))
            tmp_dir = Settings_handle::st.get("temp-folder").as<string>()+"/";

        struct stat st;
        if(stat(tmp_dir.c_str(),&st) != 0)
            tmp_dir = "";

        return tmp_dir;
    }

    void delete_files(int r);

public:
    FastTree_tree();
    bool test_executable();
    string infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein, int n_threads);
};
}

#endif // FASTTREE_TREE_H
