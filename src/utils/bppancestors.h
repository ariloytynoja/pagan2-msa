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
#include "utils/temp_dir.h"

#include <sys/stat.h>

namespace ppa
{

class BppAncestors
{
    std::string bppdistpath;

    std::string get_temp_dir()
    {
        return resolve_temp_dir();
    }

    void delete_files(int r);
    bool probe_executable();

public:
    BppAncestors();
    // Are ancestors asked for (--events, --output-ancestors, --ancestors,
    // --xml, --xml-nhx)? Whether BppAncestors may produce them is
    // test_executable()'s question, and --no-bppancestors the caller's.
    static bool ancestors_requested();
    // Can bppancestor be run? Asked ONCE per process; a process forked after
    // the first call inherits the answer.
    bool test_executable();
    bool infer_ancestors(Node *root,vector<Fasta_entry> *aligned_sequences,bool isCodon=false);
    void count_events(Node *root,vector<Fasta_entry> *aligned_sequences,string outfile,bool isCodon=false);
};
}

#endif // BPPANCESTORS_H
