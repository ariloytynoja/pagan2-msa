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

#ifndef BPPDIST_TREE_H
#define BPPDIST_TREE_H

#include "utils/settings_handle.h"
#include "utils/fasta_entry.h"
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace std;

namespace ppa {

class BppDist_tree
{
    string bppdistpath;

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
    BppDist_tree();
    bool test_executable();
    string infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein, int n_threads);
    string infer_phylogeny_fifo(std::vector<Fasta_entry> *sequences,bool is_protein, int n_threads);
};
}

#endif // BPPDIST_TREE_H
