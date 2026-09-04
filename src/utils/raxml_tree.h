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
#include "utils/temp_dir.h"
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
        std::string tmp_dir = resolve_temp_dir();

        // RAxML is invoked with "-w <dir>", which it requires to be an
        // absolute path, so this one caller canonicalises what the others
        // use verbatim.
        //
        // The buffer is allocated by realpath() rather than being a fixed
        // char[200] on the stack: a resolved path may be up to PATH_MAX
        // (4096) bytes, so the old buffer could be overflowed by a
        // sufficiently deep --temp-folder.  A NULL return is also handled
        // now; it was previously ignored and the uninitialised buffer used.
        if(!tmp_dir.empty())
        {
            char *resolved = realpath(tmp_dir.c_str(), NULL);
            if(resolved != NULL)
            {
                tmp_dir = std::string(resolved)+"/";
                free(resolved);
            }
        }

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
