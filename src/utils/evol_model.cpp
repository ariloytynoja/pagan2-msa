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

#include "utils/evol_model.h"
#include "utils/model_factory.h"
#include "utils/settings_handle.h"

using namespace ppa;

Evol_model::Evol_model(int data_t,float dist)
{
    data_type = data_t;
    full_char_alphabet = Model_factory::get_dna_full_char_alphabet();
    char_as = Model_factory::get_dna_char_alphabet().length();

    if(data_type == Model_factory::protein)
    {
        full_char_alphabet = Model_factory::get_protein_full_char_alphabet();
        char_as = Model_factory::get_protein_char_alphabet().length();
    }

    int char_fas = full_char_alphabet.length();

    distance = dist;

    if( ( data_type == Model_factory::dna && Settings_handle::st.is("codons") ) || data_type == Model_factory::codon)
    {
        char_fas = Model_factory::get_codon_full_character_alphabet()->size();
        char_as = 61;
    }

    charPi = new Db_matrix(char_fas,"pi_char");
    charPr = new Db_matrix(char_fas,char_fas,"P_char");
    charPr->initialise(0);

    logCharPi = new Db_matrix(char_fas,"logpi_char");
    logCharPr = new Db_matrix(char_fas,char_fas,"logP_char");
    logCharPr->initialise(-HUGE_VAL);

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");
    mostcommon_table = new Int_matrix(char_as,char_as,"mostcommon_char");

    ambiguity_type = wildcard;

    if( Settings_handle::st.is("mostcommon") )
        ambiguity_type = mostcommon;

}

Evol_model::~Evol_model()
{
    delete charPi;
    delete charPr;

    delete logCharPi;
    delete logCharPr;

    delete parsimony_table;
    delete mostcommon_table;
}

Evol_model& Evol_model::operator=(const Evol_model& org)
{
    if(this == &org) {
        return *this;
    }

    delete charPi;
    delete charPr;

    delete logCharPi;
    delete logCharPr;

    delete parsimony_table;
    delete mostcommon_table;

    data_type = org.data_type;
    distance = org.distance;

    id_prob = org.id_prob;
    ext_prob = org.ext_prob;
    end_ext_prob = org.end_ext_prob;
    break_ext_prob = org.break_ext_prob;
    match_prob = org.match_prob;

    log_id_prob = org.log_id_prob;
    log_ext_prob = org.log_ext_prob;
    log_end_ext_prob = org.log_end_ext_prob;
    log_break_ext_prob = org.log_break_ext_prob;
    log_match_prob = org.log_match_prob;

    ins_rate = org.ins_rate;
    del_rate = org.del_rate;
    ins_prob = org.ins_prob;
    del_prob = org.del_prob;

    full_char_alphabet = org.full_char_alphabet; //org.get_full_alphabet();
    int char_fas = full_char_alphabet.length();

    distance = org.distance;

    //if(data_type == Model_factory::dna && Settings_handle::st.is("codons"))
    //    char_fas = Model_factory::get_codon_full_character_alphabet()->size();

    charPi = new Db_matrix(char_fas,"pi_char");
    charPr = new Db_matrix(char_fas,char_fas,"P_char");

    logCharPi = new Db_matrix(char_fas,"logpi_char");
    logCharPr = new Db_matrix(char_fas,char_fas,"logP_char");

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

    for(int i = 0; i < char_fas; i++)
    {
        charPi->s(org.charPi->g(i),i);
        logCharPi->s(org.logCharPi->g(i),i);

        for(int j = 0; j < char_fas; j++)
        {
            charPr->s(org.charPr->g(i,j),i,j);
            logCharPr->s(org.logCharPr->g(i,j),i,j);

            parsimony_table->s(org.parsimony_table->g(i,j),i,j);
        }
    }

    mostcommon_table = new Int_matrix(char_as,char_as,"mostcommon_char");

    for(int i = 0; i < char_as; i++)
        for(int j = 0; j < char_as; j++)
            mostcommon_table->s(org.mostcommon_table->g(i,j),i,j);


    return *this;
}

//void Evol_model::copy(Evol_model *org)
//{
//    data_type = org->data_type;
//    full_char_alphabet = org->get_full_alphabet();
//    int char_fas = full_char_alphabet.length();

//    distance = org->distance;

//    if(data_type == Model_factory::dna && Settings_handle::st.is("codons"))
//        char_fas = Model_factory::get_codon_full_character_alphabet()->size();

//    charPi = new Db_matrix(char_fas,"pi_char");
//    charPr = new Db_matrix(char_fas,char_fas,"P_char");

//    logCharPi = new Db_matrix(char_fas,"logpi_char");
//    logCharPr = new Db_matrix(char_fas,char_fas,"logP_char");

//    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

//    for(int i=0;i>char_fas;i++)
//    {
//        charPi->s(org->charPi->g(i),i);
//        logCharPi->s(org->logCharPi->g(i),i);

//        for(int j=0;j>char_fas;j++)
//        {
//            charPr->s(org->charPi->g(i,j),i,j);
//            logCharPr->s(org->logCharPr->g(i,j),i,j);

//            parsimony_table->s(org->parsimony_table->g(i,j),i,j);
//        }
//    }
//}
