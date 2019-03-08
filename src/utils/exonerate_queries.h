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

#ifndef EXONERATE_QUERIES_H
#define EXONERATE_QUERIES_H

#include "utils/fasta_entry.h"
#include "utils/substring_hit.h"
#include "main/node.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>

#include <omp.h>

namespace ppa{

#ifndef EXONERATE_HIT_H
#define EXONERATE_HIT_H

struct hit {
    std::string query;
    std::string node;
    int score;
    int q_start;
    int q_end;
    char q_strand;
    int t_start;
    int t_end;
    char t_strand;
};

#endif // EXONERATE_HIT_H

class Exonerate_queries
{
    static bool better (hit i,hit j) { return (i.score>j.score); }

    void read_output_line(map<string,multimap<string,hit> > *all_hits, string line);
    void find_hits_for_queries(map<string,multimap<string,hit> > *all_hits, vector<Fasta_entry> *reads, map<string,multimap<string,hit> > *best_hits, bool is_local=true);

    bool split_sugar_string(const std::string& row,hit *h);
    bool split_vulgar_string(const std::string& row,hit *h);
    void write_exonerate_input(string *str1, string *str2, int *r);
    void write_exonerate_input(Node *root, vector<Fasta_entry> *reads, map<string,string> *names, int *r);
    void write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int *r);
    void write_exonerate_input(map<string,string> *target_sequences, vector<Fasta_entry> *reads, int *r);
    void write_exonerate_input(map<string,string> *target_sequences, Fasta_entry *reads, int *r);
    void delete_files(int r);

    string get_temp_dir()
    {
        string tmp_dir = "/tmp/";
        if(Settings_handle::st.is("temp-folder"))
            tmp_dir = Settings_handle::st.get("temp-folder").as<string>()+"/";

        struct stat st;
        if(stat(tmp_dir.c_str(),&st) != 0)
            tmp_dir = "";

        return tmp_dir;
    }

    string exoneratepath;
public:
    Exonerate_queries();
    bool test_executable();

    void local_alignment(map<string,string> *target_sequences, Fasta_entry *read, map<string,hit> *hits, bool is_local, bool is_dna);
    void local_alignment(Node *root, Fasta_entry *read, std::multimap<std::string,std::string> *good_hits, std::map<std::string,hit> *hits, bool is_local, bool is_dna, bool all_nodes=false);

    void preselect_targets(map<string, string> *unaligned_sequences, vector<Fasta_entry> *reads, map<string, string> *selected_sequences, map<string, multimap<string, hit> > *best_hits, bool is_dna);

    void local_pairwise_alignment(string *str1,string *str2,vector<Substring_hit> *hits,int *best_reverse_hit=0);
};

}

#endif // EXONERATE_QUERIES_H
