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

#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include <string>
#include <vector>
#include "utils/db_matrix.h"
#include "utils/int_matrix.h"
#include "utils/evol_model.h"
#include "utils/settings.h"

namespace ppa{

struct Char_symbol
{
    int index;
    char symbol;
    int n_units;
    int first_residue;
    int second_residue;
    std::vector<char> residues;
};

struct Codon_symbol
{
    int index;
    std::string symbol;
    int n_units;
    int first_codon;
    int second_codon;
    std::vector<std::string> codons;
};

class Model_factory
{

    int char_as;  // alphabet size
    int char_fas; // full_alphabet size
    int sequence_data_type;

    std::string char_alphabet;
    std::string full_char_alphabet;

    static std::string dna_char_alphabet;
    static std::string dna_full_char_alphabet;
    static std::string protein_char_alphabet;
    static std::string protein_full_char_alphabet;

    std::vector<std::string> *character_alphabet;
    std::vector<std::string>  *full_character_alphabet;

    static std::vector<std::string> dna_character_alphabet;
    static std::vector<std::string> dna_full_character_alphabet;
    static std::vector<std::string> protein_character_alphabet;
    static std::vector<std::string> protein_full_character_alphabet;
    static std::vector<std::string> codon_character_alphabet;
    static std::vector<std::string> codon_full_character_alphabet;

    static std::vector<std::string> ancestral_character_alphabet;

    static Db_matrix characterPi;

    Db_matrix *charPi;
    float char_ins_rate;
    float char_del_rate;
    float char_ext_prob;
    float char_end_ext_prob;
    float char_break_ext_prob;

    Int_matrix *parsimony_table;
    Int_matrix *child_parsimony_table;
    Int_matrix *mostcommon_table;

    Db_matrix * charU;
    Db_matrix * charV;
    Db_matrix * charRoot;

    Db_matrix *char_ambiguity;

    std::vector<Char_symbol> char_symbols;
    std::vector<Codon_symbol> codon_symbols;

    void define_dna_alphabet();
    void define_protein_alphabet();
    void define_protein_alphabet_groups();
    void define_codon_alphabet();

    void print_char_alphabet();
    void print_codon_alphabet();

    void build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot);
    void print_char_q_matrices(Db_matrix *charQ);

public:
    Model_factory(int sequence_data_type);
    ~Model_factory();

    enum Data_type {dna,protein,codon};

    int get_sequence_data_type() { return sequence_data_type; }

    void dna_model(float *char_pi,Settings *st);
    void dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob, float end_ext_prob, float break_ext_prob);

    void protein_model(Settings *st);
    void protein_model(float ins_rate,float del_rate, float ext_prob, float end_ext_prob);

    void codon_model(Settings *st);
    void codon_model(float ins_rate,float del_rate, float ext_prob, float end_ext_prob);

    Evol_model alignment_model(double distance);

    void print_int_matrix(Int_matrix *m);
    void print_char_p_matrices(Evol_model &model);

    /********************************/

    std::string get_char_alphabet() { return char_alphabet; }
    std::string get_full_char_alphabet() { return full_char_alphabet; }

    static std::string get_dna_char_alphabet() { return dna_char_alphabet; }
    static std::string get_dna_full_char_alphabet() { return dna_full_char_alphabet; }

    static std::string get_protein_char_alphabet() { return protein_char_alphabet; }
    static std::string get_protein_full_char_alphabet()
    {
        if(protein_full_char_alphabet.empty())
        {
            protein_full_char_alphabet = protein_char_alphabet+"X";

            for(int i=0;i<(int)protein_char_alphabet.length()-1;i++)
                for(int j=i+1;j<(int)protein_char_alphabet.length();j++)
                    protein_full_char_alphabet += tolower(protein_char_alphabet.at(i));
        }
        return protein_full_char_alphabet;
    }

    /********************************/

    std::vector<std::string> *get_character_alphabet() { return character_alphabet; }
    std::vector<std::string> *get_full_character_alphabet() { return full_character_alphabet; }

    static std::vector<std::string> *get_dna_character_alphabet()
    {
        if(dna_character_alphabet.empty())
        {
            for(int i=0;i<(int)dna_char_alphabet.length();i++)
                dna_character_alphabet.push_back(std::string(1,dna_char_alphabet.at(i)));
        }
        return &dna_character_alphabet;
    }

    static std::vector<std::string> *get_dna_full_character_alphabet()
    {
        if(dna_full_character_alphabet.empty())
        {
            for(int i=0;i<(int)dna_full_char_alphabet.length();i++)
                dna_full_character_alphabet.push_back(std::string(1,dna_full_char_alphabet.at(i)));
        }
        return &dna_full_character_alphabet;
    }

    static std::vector<std::string> *get_protein_character_alphabet()
    {
        if(protein_character_alphabet.empty())
        {
            for(int i=0;i<(int)protein_char_alphabet.length();i++)
                protein_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(i)));
        }
        return &protein_character_alphabet;
    }

    static std::vector<std::string> *get_protein_full_character_alphabet()
    {
        if(protein_full_character_alphabet.empty())
        {
            for(int i=0;i<(int)protein_char_alphabet.length();i++)
                protein_full_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(i)));

            protein_full_character_alphabet.push_back("X");

            for(int i=0;i<(int)protein_char_alphabet.length()-1;i++)
                for(int j=i+1;j<(int)protein_char_alphabet.length();j++)
                    protein_full_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(i))+std::string(1,protein_char_alphabet.at(j)));

        }
        return &protein_full_character_alphabet;
    }

    static std::vector<std::string> *get_codon_character_alphabet()
    {
        if(codon_character_alphabet.empty())
        {
            std::string full_alpha = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTNNN";

            for(int i=0;i<62;i++)
            {
                codon_character_alphabet.push_back(full_alpha.substr(i*3,3));
            }
        }
        return &codon_character_alphabet;
    }

    static std::vector<std::string> *get_codon_full_character_alphabet()
    {

        if(codon_full_character_alphabet.empty())
        {
            std::string full_alpha = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTNNN";

            for(int i=0;i<62;i++)
                codon_full_character_alphabet.push_back(full_alpha.substr(i*3,3));

            for(int i=0;i<60;i++)
                for(int j=i+1;j<61;j++)
                    codon_full_character_alphabet.push_back(full_alpha.substr(i*3,3)+full_alpha.substr(j*3,3));

        }
        return &codon_full_character_alphabet;
    }

    /********************************/

    static std::string get_ancestral_character_alphabet_at(int i) { return ancestral_character_alphabet.at(i); }
    /********************************/

    std::string get_character_in_full_alphabet_at(int i) { return full_character_alphabet->at(i); }

    int parsimony_state(int left_state,int right_state) { return parsimony_table->g(left_state,right_state); }
    int get_child_parsimony_state(int parent_state,int child_state) { return child_parsimony_table->g(parent_state,child_state);}
    int mostcommon_state(int left_state,int right_state) { return mostcommon_table->g(left_state,right_state); }
};

}
#endif // MODEL_FACTORY_H
