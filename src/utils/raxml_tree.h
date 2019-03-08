/***************************************************************************
 *   Copyright (C) 2010-2019 by Ari Loytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef RAXML_TREE_H
#define RAXML_TREE_H

#include "utils/settings_handle.h"
#include "utils/fasta_entry.h"
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <stdlib.h>

using namespace std;

namespace ppa {

class RAxML_tree
{
    string raxmlpath;

    std::string get_temp_dir()
    {
        std::string tmp_dir = "/tmp/";

        if(Settings_handle::st.is("temp-folder"))
        {
            tmp_dir = Settings_handle::st.get("temp-folder").as<string>()+"/";

            char resolved_path[200];
            char* rp = realpath(tmp_dir.c_str(), resolved_path);
            tmp_dir = string(resolved_path)+"/";
        }
        struct stat st;
        if(stat(tmp_dir.c_str(),&st) != 0)
            tmp_dir = "";

        return tmp_dir;
    }

    void delete_files(int r);

public:
    RAxML_tree();
    bool test_executable();
    string infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein, int n_threads);
};
}

#endif // RAXML_TREE_H
