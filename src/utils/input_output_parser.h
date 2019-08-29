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

#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <vector>
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "main/node.h"
#include "utils/model_factory.h"
#include "utils/codon_translation.h"

namespace ppa{

class Input_output_parser
{
    bool warn_tree_change = false;
public:
    Input_output_parser();

    void parse_input_sequences(ppa::Fasta_reader *fr,vector<Fasta_entry> *sequences, bool *reference_alignment);
    Node * parse_input_tree(ppa::Fasta_reader *fr,vector<Fasta_entry> *sequences, bool reference_alignment, int n_threads);
    void match_sequences_and_tree(ppa::Fasta_reader *fr,vector<Fasta_entry> *sequences, Node *root, bool reference_alignment,int *data_type);
    void define_alignment_model(ppa::Fasta_reader *fr,Model_factory *mf, int data_type);
    void output_aligned_sequences(ppa::Fasta_reader *fr,vector<Fasta_entry> *sequences, Node *root);
    void prune_extended_alignment(Fasta_reader *fr,Node *root,vector<Fasta_entry> *aligned_sequences);
    void output_pruned_alignment(Fasta_reader *fr,Node *root,Node *tmp_root,vector<Fasta_entry> *aligned_sequences,string desc, string prefix);


    void translate_codons(vector<Fasta_entry> *sequences, vector<Fasta_entry> *translated)
    {
        Codon_translation ct;
        ct.define_translation_tables();
        for(int i=0;i<(int)sequences->size();i++)
        {
            string s = sequences->at(i).sequence;
            s = ct.gapped_DNA_to_protein(&s);

            Fasta_entry e;
            e.name = sequences->at(i).name;
            e.sequence = s;

            translated->push_back(e);
        }
    }
};
}

#endif // INPUT_PARSER_H
