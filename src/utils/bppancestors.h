/***************************************************************************
 *   Copyright (C) 2013 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
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

#ifndef BPPANCESTORS_H
#define BPPANCESTORS_H

#include "utils/fasta_entry.h"
#include "main/node.h"
#include "utils/settings_handle.h"

#include <sys/stat.h>

namespace ppa
{

class BppAncestors
{
    std::string bppdistpath;

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
    BppAncestors();
    bool test_executable();
    bool infer_ancestors(Node *root,vector<Fasta_entry> *aligned_sequences,bool isCodon=false);
    void count_events(Node *root,vector<Fasta_entry> *aligned_sequences,string outfile,bool isCodon=false);
};
}

#endif // BPPANCESTORS_H
