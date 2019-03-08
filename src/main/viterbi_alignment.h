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

#ifndef VITERBI_ALIGNMENT_H
#define VITERBI_ALIGNMENT_H

#include "utils/evol_model.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/tunnel_matrix.h"
#include "main/sequence.h"
#include "main/basic_alignment.h"
#include "boost/multi_array.hpp"
#include "utils/find_anchors.h"
#include <string>
#include <stack>

namespace ppa {

using namespace std;
using namespace ppa;


class Viterbi_alignment: public Basic_alignment
{
    vector<int> upper_bound;
    vector<int> lower_bound;
    bool tunnel_defined;

    vector<Tunnel_block> empty_tunnel_blocks;

    //typedef boost::multi_array<Matrix_pointer, 2> balign_array;
    typedef Tunnel_matrix<Matrix_pointer> align_array;
    //typedef balign_array::index index;

    align_array *match;
    align_array *xgap;
    align_array *ygap;

    //typedef boost::multi_array_types::index_range range;
    //balign_array::index_gen indices;
    //typedef align_array::array_view<1>::type align_slice;
    typedef Tunnel_slice<Matrix_pointer> align_slice;

    /*********************************/

    void merge_sampled_sequence(Sequence *ancestral_sequence, Sequence *sampled_sequence);

    void initialise_array_corner();
    void initialise_array_corner_bwd();

    /*********************************/

    void compute_fwd_scores(int i,int j);
    void compute_bwd_full_score(int i,int j);
    void compute_posterior_score(int i,int j,double full_score);

    void backtrack_new_path(vector<Path_pointer> *path,Path_pointer pp);
    void sample_new_path(vector<Path_pointer> *path,Path_pointer pp);

    /*********************************/

    void iterate_bwd_edges_for_gap(Site * site,align_slice *x_slice,align_slice *y_slice,align_slice *m_slice,
                                   Matrix_pointer *max,bool is_x_matrix, int gap_type = Viterbi_alignment::normal_gap);
    void iterate_bwd_edges_for_match(Site * left_site,Site * right_site,Matrix_pointer *max);
    void iterate_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max);    
    void iterate_fwd_edges_for_gap(Site * site,align_slice *g_slice,
                                   Matrix_pointer *max_s,Matrix_pointer *max_d,Matrix_pointer *max_m);
    void iterate_fwd_edges_for_match(Site * left_site,Site * right_site,
                                     Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);

    /*********************************/

    void iterate_bwd_edges_for_sampled_gap(int site_index1,int site_index2,Matrix_pointer *bwd_p,bool is_x_matrix);
    void iterate_bwd_edges_for_sampled_match(int left_index,int right_index,Matrix_pointer *bwd_p);

    void iterate_bwd_edges_for_sampled_end_corner(Site * left_site,Site * right_site,Matrix_pointer *sample);

    void add_sample_m_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double m_match = 0);
    void add_sample_x_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double x_match = 0);
    void add_sample_y_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double y_match = 0);

    void add_sample_gap_ext(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_double(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_open(Edge * edge,align_slice *m_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_close(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);

    /*********************************/

    void score_m_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_x_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_y_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);

    void score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type = Viterbi_alignment::normal_gap);
    void score_gap_double(Edge *edge,align_slice *w_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_open(Edge *edge,align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_close(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix);

    /*********************************/

    void score_gap_ext_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max);
    void score_gap_double_bwd(Edge *edge,align_slice *w_slice,Matrix_pointer *max);
    void score_gap_open_bwd(Edge *edge,align_slice *m_slice,Matrix_pointer *max);

    void score_match_bwd(Edge * left_edge,Edge * right_edge,
                                       Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);

    /*********************************/

    void insert_gap_path_pointer(stack<Path_pointer> *path, int i, int j, int matrix,float branch_length)
    {
        Matrix_pointer mp(-1,i,j,matrix);
        if(matrix == Viterbi_alignment::x_mat)
        {
            mp.fwd_score = (*xgap)[i][j].fwd_score;
            mp.bwd_score = (*xgap)[i][j].bwd_score;
            mp.full_score = (*xgap)[i][j].full_score;
        }
        else
        {
            mp.fwd_score = (*ygap)[i][j].fwd_score;
            mp.bwd_score = (*ygap)[i][j].bwd_score;
            mp.full_score = (*ygap)[i][j].full_score;
        }
        Path_pointer pp( mp, false, branch_length,1 );
        path->push(pp);
    }

    void insert_preexisting_gap(stack<Path_pointer> *path,int *i, int *j, int x_ind, int y_ind)
    {
        bool add_one_more = false;
        bool mark_used = false;

        while(x_ind<*i)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0 && mark_used)
                left->get_edges()->at(ind).is_used(true);

            this->insert_gap_path_pointer(path,*i-1,*j,Viterbi_alignment::x_mat,left_branch_length);
            --*i;

            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0 && mark_used)
                left->get_edges()->at(ind).is_used(true);
        }

        add_one_more = false;

        while(y_ind<*j)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0 && mark_used)
                right->get_edges()->at(ind).is_used(true);

            this->insert_gap_path_pointer(path,*i,*j-1,Viterbi_alignment::y_mat,right_branch_length);
            --*j;
            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0 && mark_used)
                right->get_edges()->at(ind).is_used(true);
        }
    }


    void insert_new_path_pointer(stack<Path_pointer> *path,int *i, int *j,Path_pointer pp)
    {
        if( *i>0 || *j>0 )
            path->push(pp);
    }


    /********************************************/


    static int plot_number;
    void plot_posterior_probabilities_up();
    void plot_posterior_probabilities_down();

    string print_matrices();

    /********************************************/

public:
    Viterbi_alignment();

    float define_tunnel(Sequence *left_sequence,Sequence *right_sequence, string left_sequence_id, string right_sequence_id, Evol_model *evol_model,bool compute_coverage=false);

    void align(Sequence *left_sequence,Sequence *right_sequence,Evol_model *model,
               float left_branch_length=0,float right_branch_length=0,bool is_reads_sequence=false);
    bool replace_largest_tunnel_block_with_gap_tunnel();
    long long int get_predicted_memory_consumption(int sequence_1_length, int sequence_2_length);


};

}
#endif // VITERBI_ALIGNMENT_H
