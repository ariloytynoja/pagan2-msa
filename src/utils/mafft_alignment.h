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

#ifndef MAFFT_ALIGNMENT_H
#define MAFFT_ALIGNMENT_H

#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include "utils/settings_handle.h"
#include "utils/fasta_entry.h"

namespace ppa{

using namespace std;

class Mafft_alignment
{
    string mafftpath;

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

    std::string remove_last_whitespaces(const std::string & s)
    {
        // Copy sequence
        std::string st (s);

        while (st.size() > 0 && this->is_whitespace_character(st[st.size() - 1]))
        {
            st.erase(st.end() - 1);
        }

        // Send result
        return st;
    }

    std::string remove_whitespaces(const std::string & s)
    {
        std::string st="";

        for (unsigned int i = 0; i < s.size(); i++)
        {
            if (!this->is_whitespace_character(s[i]))
            {
                st+=s[i];
            }
        }
        return st;
    }

    bool is_whitespace_character(char c)
    {
        return (c == ' ')
               || (c == '\t')
               || (c == '\n')
               || (c == '\r')
               || (c == '\f');
    }

    void delete_files(int r);

public:
    Mafft_alignment();
    bool test_executable();
    void align_sequences(vector<Fasta_entry> *sequences);
    void align_sequences_fifo(vector<Fasta_entry> *sequences);
};
}

#endif // MAFFT_ALIGNMENT_H
