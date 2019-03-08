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

#ifndef REFERENCE_ALIGNMENT_H
#define REFERENCE_ALIGNMENT_H

#include "utils/evol_model.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
#include "main/basic_alignment.h"
#include <string>

namespace ppa {

using namespace std;
using namespace ppa;


class Reference_alignment: public Basic_alignment
{

    vector<Matrix_pointer> *mvectp;
    vector<Matrix_pointer> *xvectp;
    vector<Matrix_pointer> *yvectp;

    vector<int> *left_child_site_to_path_index_p;
    vector<int> *right_child_site_to_path_index_p;

    vector<int> *left_child_site_to_last_path_index_p;
    vector<int> *right_child_site_to_last_path_index_p;

    vector<int> *path_to_left_child_site_index_p;
    vector<int> *path_to_right_child_site_index_p;


    void make_alignment_path(vector<Matrix_pointer> *simple_path);

    void backtrack_new_vector_path(vector<Path_pointer> *path,Path_pointer fp,vector<Matrix_pointer> *simple_path);

    /*********************************/

    void iterate_bwd_edges_for_known_gap(Site * left_site,Site * right_site,vector<Matrix_pointer> *z_slice,vector<Matrix_pointer> *w_slice,
                                                     vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type,bool alignment_end=false);

    void iterate_bwd_edges_for_known_match(Site * left_site,Site * right_site,Matrix_pointer *max,int prev_mat);

    void iterate_bwd_edges_for_vector_end(Site * left_site,Site * right_site,Matrix_pointer *max,int last_matrix);

    /*********************************/

    void score_m_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max, bool allow_any=false);
    void score_x_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max);
    void score_y_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max);

    void score_gap_ext_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix,int gap_type,bool alignment_end = false);
    void score_gap_double_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *w_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_open_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix,bool alignment_end = false);
    void score_gap_close_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix);

    /*********************************/

public:
    Reference_alignment();

    void read_alignment(Sequence *left_sequence,Sequence *right_sequence,Evol_model *model,
                        float left_branch_length=0,float right_branch_length=0);

};

}
#endif // REFERENCE_ALIGNMENT_H
