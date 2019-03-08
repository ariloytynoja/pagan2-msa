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

#ifndef MATRIX_POINTER_H
#define MATRIX_POINTER_H

#include <cmath>
#include "utils/model_factory.h"
#include "utils/log_output.h"

namespace ppa {

using namespace std;
using namespace ppa;

struct Matrix_pointer
{
    double score;
    double full_score;
    double fwd_score;
    double bwd_score;
    int x_ind;
    int y_ind;
    int x_edge_ind;
    int y_edge_ind;
    int matrix;
    int path_index;

    Matrix_pointer() : score(-HUGE_VAL), full_score(0), fwd_score(0), bwd_score(0), x_ind(-1), y_ind(-1),
                        x_edge_ind(-1), y_edge_ind(-1), matrix(-1), path_index(-1) {}
    Matrix_pointer(double s,int x, int y, int m) : score(s), full_score(0), fwd_score(0), bwd_score(0),
                        x_ind(x), y_ind(y), x_edge_ind(-1), y_edge_ind(-1), matrix(m), path_index(-1) {}
};

struct Path_pointer
{

    Matrix_pointer mp;
    bool real_site;

    float branch_length_increase;
    int branch_count_increase;

    Path_pointer(Matrix_pointer m, bool real, float l, int i) : mp(m), real_site(real), branch_length_increase(l), branch_count_increase(i) {}
    Path_pointer(Matrix_pointer m, bool real) : mp(m), real_site(real), branch_length_increase(0), branch_count_increase(0) {}
    Path_pointer(bool real) : real_site(real), branch_length_increase(0), branch_count_increase(0) {}

};

struct Edge_history
{
    int left_real_site_index;
    int right_real_site_index;
    int left_skip_site_index;
    int right_skip_site_index;

    int real_site_index;
    int match_site_index;
    int path_state;

    Edge_history(int lr,int rr,int ls=-1,int rs=-1) : left_real_site_index(lr), right_real_site_index(rr),
                                                      left_skip_site_index(ls), right_skip_site_index(rs),
                                                      real_site_index(-1), match_site_index(0), path_state(-1) {}
};


}
#endif // MATRIX_POINTER_H

#ifndef BASIC_ALIGNMENT_H
#define BASIC_ALIGNMENT_H

#include "utils/evol_model.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
//#include "main/graph_reconstruction.h"
#include "boost/multi_array.hpp"
#include <string>

namespace ppa {

using namespace std;
using namespace ppa;


class Basic_alignment
{
protected:
    enum Matrix_pt {x_mat,y_mat,m_mat};
    enum Gap_type {normal_gap,end_gap,pair_break_gap};

    Sequence *left;
    Sequence *right;
    Sequence *ancestral_sequence;

    float left_branch_length;
    float right_branch_length;

    Evol_model *model;
//    string full_char_alphabet;

    vector<Path_pointer> path;

    // parameters that may be set by program arguments
    //
    float ins_del_ratio;   // (not used so far)
    float del_ins_ratio;   // (not used so far)

    float max_allowed_skip_distance;
    int   max_allowed_skip_branches;
    int   max_allowed_match_skip_branches;
    float branch_skip_weight;
    float branch_skip_probability;
    bool  weighted_branch_skip_penalty;

    bool reduced_terminal_gap_penalties;
    bool pair_end_reads;
    bool edges_for_skipped_flanked_by_gaps;
    int x_length;
    int y_length;
    int x_read1_length;
    int y_read1_length;

    bool compute_full_score;
    bool weight_edges;
    double bwd_full_probability; // for control

    void build_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path,bool is_reads_sequence=false);

    void create_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path,bool is_reads_sequence);
    void create_ancestral_edges(Sequence *sequence);
    void check_skipped_boundaries(Sequence *sequence);

    void delete_edge_range(Sequence *sequence,int edge_ind,int skip_start_site);

    void transfer_child_edge(Sequence *sequence, Edge *child, vector<int> *child_index, float branch_length,
                             bool adjust_posterior_weight = true, float branch_weight = 1.0);
    void transfer_child_edge(Sequence *sequence, Edge edge, Edge *child, float branch_length,
                             bool adjust_posterior_weight = true, float branch_weight = 1.0);

    /*********************************/

    void compute_site_consensus(Site *site,Sequence *left,int l_pos,Sequence *right,int r_pos,bool is_dna)
    {
        int ndl = 1;
        int ndr = 1;

        if(Settings_handle::st.is("use-duplicate-weigths"))
        {
            if(left->is_terminal_sequence())
                ndl = left->get_num_duplicates();
            if(right->is_terminal_sequence())
                ndr = right->get_num_duplicates();
        }
//        cout<<"\n"<<l_pos<<" "<<ndl<<"; "<<r_pos<<" "<<ndr<<" "<<left->is_read_sequence()<<" "<<right->is_read_sequence()<<endl;
        if(!is_dna)
        {
/*
            // old->
            if(l_pos>=0)
            {
                if(!left->is_terminal_sequence())
                    site->add_sumAmino(left->get_site_at(l_pos)->get_sumAmino());
                else if(left->is_read_sequence())
                    site->add_sumAmino(ndl);
            }
            if(r_pos>=0)
            {
                if(!right->is_terminal_sequence())
                    site->add_sumAmino(right->get_site_at(r_pos)->get_sumAmino());
                else if(right->is_read_sequence())
                    site->add_sumAmino(ndr);
            }
            // <-old

            for(int i=0;i<20;i++)
                site->set_sumAA(i,0);

            if(l_pos>=0)
            {
                if(left->is_read_descendants())
                {
                    for(int i=0;i<20;i++)
                    {
                        site->add_sumAA(i,left->get_site_at(l_pos)->get_sumAA(i));
                    }
                }
                else if(left->is_read_sequence())
                {
                    site->add_sumAA(left->get_site_at(l_pos)->get_state(),ndl);
                }
            }
            if(r_pos>=0)
            {
                if(right->is_read_descendants())
                {
                    for(int i=0;i<20;i++)
                    {
                        site->add_sumAA(i,right->get_site_at(r_pos)->get_sumAA(i));
                    }
                }
                else if(right->is_read_sequence())
                {
                    site->add_sumAA(right->get_site_at(r_pos)->get_state(),ndr);
                }
            }

            int max_c=1;
            int max_i=site->get_state();
            for(int i=0;i<20;i++)
            {
                int c = site->get_sumAA(i);

                if(c>max_c)
                {
                    max_i = i;
                    max_c = c;
                }
                else if(max_c == c)
                {
                    max_i = model->parsimony_state(max_i,i);
                }
            }
            site->set_state(max_i);

            return;
*/
            return;
        }

        int lsA = 0; int lsC = 0; int lsG = 0; int lsT = 0;
        int rsA = 0; int rsC = 0; int rsG = 0; int rsT = 0;

        if(l_pos>=0)
        {
            if(!left->is_terminal_sequence())
            {
                lsA += left->get_site_at(l_pos)->get_sumA();
                lsC += left->get_site_at(l_pos)->get_sumC();
                lsG += left->get_site_at(l_pos)->get_sumG();
                lsT += left->get_site_at(l_pos)->get_sumT();
            }
            else if(left->is_read_sequence())
            {
                int s = left->get_site_at(l_pos)->get_state();
                if(s == 0)
                    lsA += ndl;
                else if(s == 1)
                    lsC += ndl;
                else if(s == 2)
                    lsG += ndl;
                else if(s == 3)
                    lsT += ndl;
                else if(s == 4) {
                    lsA += ndl;
                    lsG += ndl;    }
                else if(s == 5) {
                    lsC += ndl;
                    lsT += ndl;    }
                else if(s == 6) {
                    lsA += ndl;
                    lsC += ndl;    }
                else if(s == 7) {
                    lsG += ndl;
                    lsT += ndl;    }
                else if(s == 8) {
                    lsA += ndl;
                    lsT += ndl;    }
                else if(s == 9) {
                    lsC += ndl;
                    lsG += ndl;    }
                else if(s == 10) {
                    lsC += ndl;
                    lsG += ndl;
                    lsT += ndl;    }
                else if(s == 11) {
                    lsA += ndl;
                    lsG += ndl;
                    lsT += ndl;    }
                else if(s == 12) {
                    lsA += ndl;
                    lsC += ndl;
                    lsT += ndl;    }
                else if(s == 13) {
                    lsA += ndl;
                    lsC += ndl;
                    lsG += ndl;    }
                else if(s == 14) {
                    lsA += ndl;
                    lsC += ndl;
                    lsG += ndl;
                    lsT += ndl;    }
                else
                    Log_output::write_out("Basic_alignment: compute_site_consensus: no such option (l)\n",1);
            }
        }

        if(r_pos>=0)
        {
            if(!right->is_terminal_sequence())
            {
                rsA += right->get_site_at(r_pos)->get_sumA();
                rsC += right->get_site_at(r_pos)->get_sumC();
                rsG += right->get_site_at(r_pos)->get_sumG();
                rsT += right->get_site_at(r_pos)->get_sumT();
            }
            else if(right->is_read_sequence())
            {
                int s = right->get_site_at(r_pos)->get_state();
                if(s == 0)
                    rsA += ndr;
                else if(s == 1)
                    rsC += ndr;
                else if(s == 2)
                    rsG += ndr;
                else if(s == 3)
                    rsT += ndr;
                else if(s == 4) {
                    rsA += ndr;
                    rsG += ndr;    }
                else if(s == 5) {
                    rsC += ndr;
                    rsT += ndr;    }
                else if(s == 6) {
                    rsA += ndr;
                    rsC += ndr;    }
                else if(s == 7) {
                    rsG += ndr;
                    rsT += ndr;    }
                else if(s == 8) {
                    rsA += ndr;
                    rsT += ndr;    }
                else if(s == 9) {
                    rsC += ndr;
                    rsG += ndr;    }
                else if(s == 10) {
                    rsC += ndr;
                    rsG += ndr;
                    rsT += ndr;    }
                else if(s == 11) {
                    rsA += ndr;
                    rsG += ndr;
                    rsT += ndr;    }
                else if(s == 12) {
                    rsA += ndr;
                    rsC += ndr;
                    rsT += ndr;    }
                else if(s == 13) {
                    rsA += ndr;
                    rsC += ndr;
                    rsG += ndr;    }
                else if(s == 14) {
                    rsA += ndr;
                    rsC += ndr;
                    rsG += ndr;
                    rsT += ndr;    }
                else
                    Log_output::write_out("Basic_alignment: compute_site_consensus: no such option (r)\n",1);
            }
        }

        if(lsA+lsC+lsG+lsT+rsA+rsC+rsG+rsT>0)
        {
            int sA = lsA+rsA;
            int sC = lsC+rsC;
            int sG = lsG+rsG;
            int sT = lsT+rsT;
//            cout<<sA<<","<<sC<<","<<sG<<","<<sT<<"\n";
            site->add_sumA(sA);
            site->add_sumC(sC);
            site->add_sumG(sG);
            site->add_sumT(sT);

            if(Settings_handle::st.is("use-consensus"))
            {
                if(sA>sC && sA>sG && sA>sT)
                    site->set_state(0);
                else if(sC>sA && sC>sG && sC>sT)
                    site->set_state(1);
                else if(sG>sA && sG>sC && sG>sT)
                    site->set_state(2);
                else if(sT>sA && sT>sC && sT>sG)
                    site->set_state(3);
                else if(sA>sC && sA==sG && sA>sT)
                    site->set_state(4);
                else if(sC>sA && sC>sG && sC==sT)
                    site->set_state(5);
                else if(sA==sC && sA>sG && sA>sT)
                    site->set_state(6);
                else if(sG>sA && sG>sC && sG==sT)
                    site->set_state(7);
                else if(sA>sC && sA>sG && sA==sT)
                    site->set_state(8);
                else if(sC>sA && sC==sG && sC>sT)
                    site->set_state(9);
                else if(sC>sA && sC==sG && sC==sT)
                    site->set_state(10);
                else if(sA>sC && sA==sG && sA==sT)
                    site->set_state(11);
                else if(sA==sC && sA>sG && sA==sT)
                    site->set_state(12);
                else if(sA==sC && sA==sG && sA>sT)
                    site->set_state(13);
                else if(sA==sC && sA==sG && sA==sT)
                    site->set_state(14);
                else
                    Log_output::write_out("Basic_alignment: compute_site_consensus: no such option (s) "+
                                Log_output::itos(sA)+" "+Log_output::itos(sC)+" "+Log_output::itos(sG)+" "+Log_output::itos(sT)+" "+"\n",1);
            }
        }
    }

    /*********************************/

    std::string itos(int i) // convert int to string
    {
        std::stringstream s;
        s << i;
        return s.str();
    }

    std::string ftos(float f) // convert float to string
    {
        std::stringstream s;
        s << f;
        return s.str();
    }

    /*********************************/

    bool first_is_bigger(double a,double b)
    {
        if (a==-HUGE_VAL && b==-HUGE_VAL)
            return false;
        else if (a>b)
            return true;
        else if (a<b)
            return false;
        else
            return false;
//            if ((double)rand()/(double)RAND_MAX>0.5)
//                return true;
        return false;
    }


    /********************************************/

    double get_log_edge_weight(Edge *edge) { return (*this.*log_edge_weight)(edge); }
    double (ppa::Basic_alignment::*log_edge_weight)(Edge *edge);

    double get_edge_weight(Edge *edge) { return (*this.*edge_weight)(edge); }
    double (ppa::Basic_alignment::*edge_weight)(Edge *edge);

    double edge_posterior_weight(Edge *edge) { return edge->get_posterior_weight(); }
    double edge_log_posterior_weight(Edge *edge) { return edge->get_log_posterior_weight(); }

    double edge_equal_weight(Edge *) const { return 1.0; }
    double edge_log_equal_weight(Edge *) const { return 0.0; }

    /********************************************/

    double get_transformed_edge_weight(double w) { return (*this.*transform_edge_weight)(w); }
    double (ppa::Basic_alignment::*transform_edge_weight)(double w);

    double square_root_edge_weight(double w) { return sqrt( w );}
    double cube_root_edge_weight(double w) { return exp((1.0/3.0)*(log(w)));}
    double plain_edge_weight(double w) { return w;}

    /********************************************/

    float get_log_gap_open_penalty(int prev_site, bool is_x_matrix)
    {
        if(reduced_terminal_gap_penalties)
        {
            if(prev_site==0)
            {
                return 0;
            }

            if(pair_end_reads)
            {
                if(is_x_matrix && prev_site == x_read1_length)
                {
                    return 0;
                }
                else if(!is_x_matrix && prev_site == y_read1_length)
                {
                    return 0;
                }
            }
        }

        return model->log_gap_open();
    }

    float get_log_gap_close_penalty(int this_site, bool is_x_matrix)
    {
        if(reduced_terminal_gap_penalties)
        {
            if(is_x_matrix && this_site==x_length)
            {
                return 0;
            }
            else if(!is_x_matrix && this_site==y_length)
            {
                return 0;
            }

            if(pair_end_reads)
            {
                if(is_x_matrix && this_site == x_read1_length+1)
                {
                    return 0;
                }
                else if(!is_x_matrix && this_site == y_read1_length+1)
                {
                    return 0;
                }
            }
        }

        return model->log_gap_close();
    }

    /********************************************/

    void set_basic_settings()
    {
        del_ins_ratio = 1;//model->del_prob/model->ins_prob;
        ins_del_ratio = 1;//model->ins_prob/model->del_prob;

        edges_for_skipped_flanked_by_gaps = false;

        weighted_branch_skip_penalty = false; // by default, use penalty *per node*

        max_allowed_skip_distance = 0.5;
        max_allowed_skip_branches = 10;
        max_allowed_match_skip_branches = 5;

        branch_skip_weight = 1;
        branch_skip_probability = 0.9; // 0.2

        weight_edges = false;
        compute_full_score = false;
        reduced_terminal_gap_penalties = false;
        pair_end_reads = false;

        x_length = -1;
        y_length = -1;

    }

    void set_reads_alignment_settings()
    {

        max_allowed_skip_distance = 5;
        max_allowed_skip_branches = 50000;
        max_allowed_match_skip_branches = 50000;

        branch_skip_weight = 1;
        branch_skip_probability = 1;

        //branch_skip_probability = 0.2;

//        if(Settings_handle::st.is("pair-end"))
//            pair_end_reads = true;
    }

    void set_reference_alignment_settings()
    {
        max_allowed_skip_distance = 5;
        max_allowed_skip_branches = 50000;
        max_allowed_match_skip_branches = 50000;
    }

    void set_additional_settings()
    {
        if(Settings_handle::st.is("branch-length-confirm-insertion"))
            max_allowed_skip_distance = Settings_handle::st.get("branch-length-confirm-insertion").as<float>();

        if(Settings_handle::st.is("any-skips-confirm-insertion"))
            max_allowed_skip_branches = Settings_handle::st.get("any-skips-confirm-insertion").as<int>();

        if(Settings_handle::st.is("match-skips-confirm-insertion"))
            max_allowed_match_skip_branches = Settings_handle::st.get("match-skips-confirm-insertion").as<int>();

        if(Settings_handle::st.is("branch-skip-weight-per-distance"))
        {
            branch_skip_weight = Settings_handle::st.get("branch-skip-weight-per-distance").as<float>();
            weighted_branch_skip_penalty = true;
        }

        if(Settings_handle::st.is("branch-skip-penalty-per-branch"))
        {
            branch_skip_probability = Settings_handle::st.get("branch-skip-penalty-per-branch").as<float>();
            weighted_branch_skip_penalty = false;
        }

        if( Settings_handle::st.is("weight-sampled-edges") && Settings_handle::st.get("sample-additional-paths").as<int>() > 0)
            weight_edges = true;

        if( Settings_handle::st.is("full-probability") ||
            Settings_handle::st.is("mpost-posterior-plot-file") ||
            Settings_handle::st.is("sample-path") ||
            Settings_handle::st.get("sample-additional-paths").as<int>() > 0 )
            compute_full_score = true;

        if( !Settings_handle::st.is("no-reduced-terminal-penalties") )
            reduced_terminal_gap_penalties = true;
    }

    /********************************************/

    void mark_no_gap_penalty_sites(Sequence *left, Sequence *right)
    {

        x_length = left->sites_length();
        y_length = right->sites_length();
        x_read1_length = -1;
        y_read1_length = -1;

        if(pair_end_reads)
        {
            // this is a clumsy way to transfer the information of the break point
            // in a combined pair-end read but has to for now
            for(int i=0;i<left->sites_length();i++)
            {
                if( left->get_site_at(i)->get_site_type() == Site::break_start_site )
                {
                    x_read1_length = i;
                    left->get_site_at(i)->set_site_type( Site::real_site );
                }

                if( left->get_site_at(i)->get_site_type() == Site::break_stop_site )
                {
                    left->get_site_at(i)->set_site_type( Site::real_site );
                    break;
                }
            }

            for(int i=0;i<right->sites_length();i++)
            {
                if( right->get_site_at(i)->get_site_type() == Site::break_start_site )
                {
                    y_read1_length = i;
                    right->get_site_at(i)->set_site_type( Site::real_site );
                }

                if( right->get_site_at(i)->get_site_type() == Site::break_stop_site )
                {
                    right->get_site_at(i)->set_site_type( Site::real_site );
                    break;
                }
            }
        }

    }

    /********************************************/

    string print_pairwise_alignment(vector<Site> *sites);

    void debug_print_input_sequences(int noise_level)
    {

        if(Settings::noise<noise_level)
            return;

        vector<string> *full_char_alphabet = Model_factory::get_dna_full_character_alphabet();

        if(model->get_data_type() == Model_factory::protein)
            full_char_alphabet = Model_factory::get_protein_full_character_alphabet();
        else if( ( model->get_data_type() == Model_factory::dna && Settings_handle::st.is("codons") ) || model->get_data_type() == Model_factory::codon)
            full_char_alphabet = Model_factory::get_codon_full_character_alphabet();

        stringstream ss;
        ss<<"sequences:"<<endl<<" ";
        for(int i=1;i<left->sites_length()-1;i++)
            ss<<Model_factory::get_ancestral_character_alphabet_at(left->get_site_at(i)->get_state());
        ss<<endl<<" ";
        for(int i=1;i<right->sites_length()-1;i++)
        {
            ss<<Model_factory::get_ancestral_character_alphabet_at(right->get_site_at(i)->get_state());
        }
        ss<<endl;

        if(Settings::noise>4)
        {
            ss<<"\nLEFT\n"<<left->print_sequence(left->get_sites());
            ss<<"\nRIGHT\n"<<right->print_sequence(right->get_sites());
            ss<<endl;
        }
        Log_output::write_out(ss.str(),noise_level);
    }

    string print_path(vector<Path_pointer> *path)
    {
        stringstream ss;
        ss<<endl;
        for(unsigned int i=0;i<path->size();i++)
        {
            Path_pointer *tsite =  &path->at(i);
            if(tsite->real_site)
                ss<<i<<" r: ";
            else
                ss<<i<<" s: ";

            ss<<tsite->mp.x_ind<<" "<<tsite->mp.y_ind<<" "<<tsite->mp.matrix<<": "<<tsite->mp.score<<" "<<
                   fixed<<log(tsite->mp.full_score)<<" "<<fixed<<log(tsite->mp.bwd_score)<<endl;
        }
        return ss.str();
    }

    /********************************************/

public:
    Basic_alignment();

    Sequence* get_simple_sequence() { return ancestral_sequence; }

};

}
#endif // BASIC_ALIGNMENT_H
