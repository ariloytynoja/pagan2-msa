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

#include "main/reference_alignment.h"
#include "main/sequence.h"
#include "utils/exceptions.h"
#include "main/node.h"
#include <iomanip>
#include <fstream>
#include <stack>

using namespace std;
using namespace ppa;

Reference_alignment::Reference_alignment(){ reduced_terminal_gap_penalties = false; }


/********************************************/

void Reference_alignment::read_alignment(Sequence *left_sequence,Sequence *right_sequence,Evol_model *evol_model,float l_branch_length,float r_branch_length)
{

    left = left_sequence;
    right = right_sequence;
    model = evol_model;

    left_branch_length = l_branch_length;
    right_branch_length = r_branch_length;


    this->debug_print_input_sequences(3);


    // set the basic parameters (copy from Settings)
    //
    this->set_basic_settings();
    if(!Settings_handle::st.is("perfect-reference"))
        this->set_reference_alignment_settings();

    this->set_additional_settings();

    // Set the edge weighting scheme, define the dynamic-programming matrices
    //
    log_edge_weight = &ppa::Reference_alignment::edge_log_posterior_weight;
    edge_weight = &ppa::Reference_alignment::edge_posterior_weight;

    transform_edge_weight = &ppa::Reference_alignment::square_root_edge_weight;
    if( Settings_handle::st.is("no-weight-transform") )
        transform_edge_weight = &ppa::Reference_alignment::plain_edge_weight;
    if( Settings_handle::st.is("cuberoot-weight-transform") )
        transform_edge_weight = &ppa::Reference_alignment::cube_root_edge_weight;

    int left_length = left->sites_length();    int right_length = right->sites_length();

    Log_output::write_out("Reference_alignment: lengths: "+this->itos(left_length)+" "+this->itos(right_length),3);

    string gapped_anc = "";

    string *gapped_left = left->get_gapped_sequence();
    string *gapped_right = right->get_gapped_sequence();


    if(gapped_left->length() != gapped_right->length())
    {
        Log_output::write_out("Error: gapped sequences of different length. Exiting.\n\n",0);
        exit(1);
    }


    bool is_codons = Settings_handle::st.is("codons") || evol_model->get_data_type()==Model_factory::codon;

    string::iterator lgi = gapped_left->begin();
    string::iterator rgi = gapped_right->begin();
    int li = 0;
    int ri = 0;

    vector<Matrix_pointer> simple_path;

    for(;lgi!=gapped_left->end();lgi++,rgi++)
    {
        bool lgap = (*lgi == '-');
        bool rgap = (*rgi == '-');

        if(is_codons)
        {
            lgi++;rgi++;
            bool lgap2 = (*lgi == '-');
            bool rgap2 = (*rgi == '-');
            lgi++;rgi++;
            bool lgap3 = (*lgi == '-');
            bool rgap3 = (*rgi == '-');

            if(! ( ( (lgap && lgap2 && lgap3) || (!lgap && !lgap2 && !lgap3) ) &&
                   ( (rgap && rgap2 && rgap3) || (!rgap && !rgap2 && !rgap3) ) )
              )
            {
               Log_output::write_out("Reading frame error in a codon refrence alignment. Exiting.\n",0);
               exit(1);
            }
        }

        if(!lgap && rgap)
        {
            gapped_anc.append("A");
            if(is_codons)
                gapped_anc.append("CC");

            Matrix_pointer bwd_p;
            bwd_p.matrix = Reference_alignment::x_mat;
            bwd_p.x_ind = -1;
            bwd_p.y_ind = li;

            simple_path.push_back(bwd_p);

            li++;
        }
        else if(lgap && !rgap)
        {
            gapped_anc.append("A");
            if(is_codons)
                gapped_anc.append("GG");

            Matrix_pointer bwd_p;
            bwd_p.matrix = Reference_alignment::y_mat;
            bwd_p.x_ind = ri;
            bwd_p.y_ind = -1;

            simple_path.push_back(bwd_p);

            ri++;
        }
        else if(!lgap && !rgap)
        {
            gapped_anc.append("A");
            if(is_codons)
                gapped_anc.append("TT");

            Matrix_pointer bwd_p;
            bwd_p.matrix = Reference_alignment::m_mat;
            bwd_p.x_ind = ri;
            bwd_p.y_ind = li;

            simple_path.push_back(bwd_p);

            li++;
            ri++;
        }
        else if(lgap && rgap)
        {
            gapped_anc.append("-");
            if(is_codons)
                gapped_anc.append("--");
        }
    }


    // Now build the sequence forward following the path saved in a vector;
    //
    ancestral_sequence = new Sequence(path.size(),model->get_data_type(),gapped_anc);

    this->make_alignment_path(&simple_path);


    Log_output::write_out("Reference_alignment: sequence built",3);

}

void Reference_alignment::make_alignment_path(vector<Matrix_pointer> *simple_path)
{

    int left_length = left->sites_length();    int right_length = right->sites_length();

    Log_output::write_out("Reference_alignment::make_alignment_path lengths: "+this->itos(left_length)+" "+this->itos(right_length),3);

    vector<Matrix_pointer> mvect;
    vector<Matrix_pointer> xvect;
    vector<Matrix_pointer> yvect;

    mvectp = &mvect;
    xvectp = &xvect;
    yvectp = &yvect;

    double small = -HUGE_VAL;

    Matrix_pointer mpm(0.0,0,0,-1);
    mvect.push_back(mpm);

    Matrix_pointer mpg(small,0,0,-1);
    xvect.push_back(mpg);
    yvect.push_back(mpg);

    vector<int> left_child_site_to_path_index;
    vector<int> right_child_site_to_path_index;

    left_child_site_to_path_index_p = &left_child_site_to_path_index;
    right_child_site_to_path_index_p = &right_child_site_to_path_index;

    left_child_site_to_path_index.push_back(0);
    right_child_site_to_path_index.push_back(0);

    vector<int> left_child_site_to_last_path_index;
    vector<int> right_child_site_to_last_path_index;

    left_child_site_to_last_path_index_p = &left_child_site_to_last_path_index;
    right_child_site_to_last_path_index_p = &right_child_site_to_last_path_index;

    left_child_site_to_last_path_index.push_back(0);
    right_child_site_to_last_path_index.push_back(0);


    vector<int> path_to_left_child_site_index;
    vector<int> path_to_right_child_site_index;

    path_to_left_child_site_index_p = &path_to_left_child_site_index;
    path_to_right_child_site_index_p = &path_to_right_child_site_index;

    path_to_left_child_site_index.push_back(0);
    path_to_right_child_site_index.push_back(0);


    int i_ind = 0; int j_ind = 0;

    Site * left_site;
    Site * right_site;

    Matrix_pointer mp;

    int prev_mat = -1;
    int last_m_path_index = -1;

    int prev_i_ind = 0;
    int prev_j_ind = 0;

    bool i_seq_start = true;
    bool j_seq_start = true;


    for(int i=0;i<(int) simple_path->size();i++)
    {

        int j_gap_type = Reference_alignment::normal_gap;
        int i_gap_type = Reference_alignment::normal_gap;

        if( j_ind==0 || j_ind == right_length-1 )
        {
            j_gap_type = Reference_alignment::end_gap;
        }

        if( i_ind==0 || i_ind == left_length-1 )
        {
            i_gap_type = Reference_alignment::end_gap;
        }


        Matrix_pointer mpm;
        Matrix_pointer mpx;
        Matrix_pointer mpy;

        Matrix_pointer *max_x = &mpx;
        Matrix_pointer *max_y = &mpy;
        Matrix_pointer *max_m = &mpm;

        mp = simple_path->at(i);


        if(mp.matrix == Reference_alignment::x_mat)
        {

            i_ind++;

            left_child_site_to_path_index.push_back(i+1);
            left_child_site_to_last_path_index.push_back(i+1);

            left_site = left->get_site_at(i_ind);
            right_site = right->get_site_at(j_ind);

            this->iterate_bwd_edges_for_known_gap(left_site,right_site,&xvect,&yvect,&mvect,max_x,true,j_gap_type,j_seq_start);


            if(max_x->y_ind<0)
                max_x->y_ind = path_to_right_child_site_index_p->at(left_child_site_to_path_index_p->at(max_x->x_ind));

            if(max_x->matrix == Reference_alignment::y_mat)
                max_x->y_ind = path_to_right_child_site_index_p->at(left_child_site_to_last_path_index_p->at(max_x->x_ind));



            i_seq_start = false;

            prev_i_ind = i_ind;

        }

        else if(mp.matrix == Reference_alignment::y_mat)
        {

            j_ind++;

            right_child_site_to_path_index.push_back(i+1);
            right_child_site_to_last_path_index.push_back(i+1);

            left_site = left->get_site_at(i_ind);
            right_site = right->get_site_at(j_ind);

            this->iterate_bwd_edges_for_known_gap(left_site,right_site,&yvect,&xvect,&mvect,max_y,false,i_gap_type,i_seq_start);


            if(max_y->x_ind<0)
                max_y->x_ind = path_to_left_child_site_index_p->at(right_child_site_to_path_index_p->at(max_y->y_ind));

            if(max_y->matrix == Reference_alignment::x_mat)
                max_y->x_ind = path_to_left_child_site_index_p->at(right_child_site_to_last_path_index_p->at(max_y->y_ind));



            j_seq_start = false;
            prev_j_ind = j_ind;
        }

        else if(mp.matrix == Reference_alignment::m_mat)
        {
            i_ind++;
            j_ind++;

            left_child_site_to_path_index.push_back(i+1);
            right_child_site_to_path_index.push_back(i+1);

            left_child_site_to_last_path_index.push_back(i+1);
            right_child_site_to_last_path_index.push_back(i+1);

            left_site = left->get_site_at(i_ind);
            right_site = right->get_site_at(j_ind);

            this->iterate_bwd_edges_for_known_match(left_site,right_site,max_m,last_m_path_index);


            i_seq_start = false;
            j_seq_start = false;

            prev_i_ind = i_ind;
            prev_j_ind = j_ind;

            last_m_path_index = i;

        }

        prev_mat = mp.matrix;

        mvect.push_back(mpm);
        xvect.push_back(mpx);
        yvect.push_back(mpy);

        path_to_left_child_site_index.push_back(i_ind);
        path_to_right_child_site_index.push_back(j_ind);

        left_child_site_to_last_path_index.at(i_ind) = i+1;
        right_child_site_to_last_path_index.at(j_ind) = i+1;

    }


    left_child_site_to_path_index.push_back(simple_path->size());
    right_child_site_to_path_index.push_back(simple_path->size());

    left_child_site_to_last_path_index.push_back(simple_path->size());
    right_child_site_to_last_path_index.push_back(simple_path->size());


    Log_output::write_out("Reference_alignment::make_alignment_path vector filled",3);

    path.clear();

    // Find the incoming edge in the end corner; also, get the full_fwd probability
    //
    left_site = left->get_site_at(left_length-1);
    right_site = right->get_site_at(right_length-1);


    Matrix_pointer max_end;
    this->iterate_bwd_edges_for_vector_end(left_site,right_site,&max_end,mp.matrix);


    if(max_end.score==-HUGE_VAL) {
        stringstream ss;
        ss<<"Reference_alignment: ax_end.score==-HUGE_VAL\n;";
        ss<<"left\n"<<this->left->print_sequence();
        ss<<"right\n"<<this->right->print_sequence();
    }

    Log_output::write_out("Reference_alignment::make_alignment_path end found",3);

    // Backtrack the Viterbi path
    Path_pointer pp(max_end,true);
    this->backtrack_new_vector_path(&path,pp,simple_path);


    Log_output::write_out("Reference_alignment::make_alignment_path path found",3);

    // Now build the sequence forward following the path saved in a vector;
    //
    this->build_ancestral_sequence(ancestral_sequence,&path);

    Log_output::write_out("Reference_alignment::make_alignment_path sequence built",3);

}


/********************************************/

void Reference_alignment::backtrack_new_vector_path(vector<Path_pointer> *path,Path_pointer fp,vector<Matrix_pointer> *simple_path)
{

    vector<Edge> *left_edges = left->get_edges();
    vector<Edge> *right_edges = right->get_edges();

    int path_length = simple_path->size();

    int vit_mat = fp.mp.matrix;
    int x_ind = fp.mp.x_ind;
    int y_ind = fp.mp.y_ind;
    int next_path_index = fp.mp.path_index;

    if(fp.mp.x_edge_ind>=0)
        left_edges->at(fp.mp.x_edge_ind).is_used(true);
    if(fp.mp.y_edge_ind>=0)
        right_edges->at(fp.mp.y_edge_ind).is_used(true);


    // Pre-existing gaps in the end skipped over
    //
    int k = path_length;


    if(vit_mat==Reference_alignment::x_mat)
    {
        y_ind = -1;
    }
    else if(vit_mat==Reference_alignment::y_mat)
    {
        x_ind = -1;
    }

    // Uses temporary stack to decrease time consumption
    stack<Path_pointer> path_stack;

    // Actual alignment path
    //
    while(k>=0)
    {

        if(vit_mat == Reference_alignment::m_mat)
        {
            // Pre-existing gaps in the middle skipped over
            //
            while(next_path_index<k)
            {
                Matrix_pointer *smp = &simple_path->at(k-1);
                Matrix_pointer mp(-1,smp->y_ind,smp->x_ind,smp->matrix);
                Path_pointer pp( mp, false );
                //path->insert(path->begin(),pp);
                path_stack.push(pp);

                k--;
            }
            if(k<1)
                break;

            Matrix_pointer mp(-1,x_ind,y_ind,vit_mat);
            Path_pointer pp( mp, true );
            //path->insert(path->begin(),pp);
            path_stack.push(pp);


            vit_mat = (*mvectp)[k].matrix;
            x_ind = (*mvectp)[k].x_ind;
            y_ind = (*mvectp)[k].y_ind;
            next_path_index = (*mvectp)[k].path_index;

            if((*mvectp)[k].x_edge_ind>=0)
                left_edges->at((*mvectp)[k].x_edge_ind).is_used(true);
            if((*mvectp)[k].y_edge_ind>=0)
                right_edges->at((*mvectp)[k].y_edge_ind).is_used(true);

            if(vit_mat==Reference_alignment::x_mat)
            {
                y_ind = -1;
            }
            else if(vit_mat==Reference_alignment::y_mat)
            {
                x_ind = -1;
            }

            k--;

        }
        else if(vit_mat == Reference_alignment::x_mat)
        {
            // Pre-existing gaps in the middle skipped over
            //
            while(next_path_index<k)
            {
                Matrix_pointer *smp = &simple_path->at(k-1);
                Matrix_pointer mp(-1,smp->y_ind,smp->x_ind,smp->matrix);
                Path_pointer pp( mp, false );
                //path->insert(path->begin(),pp);
                path_stack.push(pp);

                k--;
            }
            if(k<1)
                break;

            Matrix_pointer mp(-1,x_ind,y_ind,vit_mat);
            Path_pointer pp( mp, true );
            //path->insert(path->begin(),pp);
            path_stack.push(pp);


            vit_mat = (*xvectp)[k].matrix;
            x_ind = (*xvectp)[k].x_ind;
            y_ind = (*xvectp)[k].y_ind;
            next_path_index = (*xvectp)[k].path_index;


            if((*xvectp)[k].x_edge_ind>=0)
                left_edges->at((*xvectp)[k].x_edge_ind).is_used(true);

            if(vit_mat==Reference_alignment::x_mat)
            {
                y_ind = -1;
            }
            else if(vit_mat==Reference_alignment::y_mat)
            {
                x_ind = -1;
            }

            k--;


        }
        else if(vit_mat == Reference_alignment::y_mat)
        {
            // Pre-existing gaps in the middle skipped over
            //
            while(next_path_index<k)
            {
                Matrix_pointer *smp = &simple_path->at(k-1);
                Matrix_pointer mp(-1,smp->y_ind,smp->x_ind,smp->matrix);
                Path_pointer pp( mp, false );
                //path->insert(path->begin(),pp);
                path_stack.push(pp);

                k--;
            }

            if(k<1)
                break;

            Matrix_pointer mp(-1,x_ind,y_ind,vit_mat);
            Path_pointer pp( mp, true );
            //path->insert(path->begin(),pp);
            path_stack.push(pp);


            vit_mat = (*yvectp)[k].matrix;
            x_ind = (*yvectp)[k].x_ind;
            y_ind = (*yvectp)[k].y_ind;
            next_path_index = (*yvectp)[k].path_index;


            if((*yvectp)[k].y_edge_ind>=0)
                right_edges->at((*yvectp)[k].y_edge_ind).is_used(true);

            if(vit_mat==Reference_alignment::x_mat)
            {
                y_ind = -1;
            }
            else if(vit_mat==Reference_alignment::y_mat)
            {
                x_ind = -1;
            }

            k--;


        }
        else
        {
            Log_output::write_out("Refrence_alignment::incorrect backward pointer: "+this->itos(vit_mat)+".\n",0);
            exit(1);
        }

        if(k<1)
            break;

    }

    path->reserve(path_stack.size());
    while(!path_stack.empty()){
        path->push_back(path_stack.top());
        path_stack.pop();
    }

}

/********************************************/


void Reference_alignment::iterate_bwd_edges_for_known_match(Site * left_site,Site * right_site,Matrix_pointer *max,int last_m_path_index)
{

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        // match score & gap close scores for this match
        //
        double log_match_score = model->log_score(left_site->character_state,right_site->character_state);
        double m_log_match = 2*model->log_non_gap() + log_match_score;

        double x_log_match = this->get_log_gap_close_penalty(left_edge->get_end_site_index(), true) + model->log_non_gap() + log_match_score;
        double y_log_match = this->get_log_gap_close_penalty(right_edge->get_end_site_index(), false) + model->log_non_gap() + log_match_score;

        if( left_child_site_to_path_index_p->at( left_edge->get_start_site_index() ) >= last_m_path_index )
        {
            this->score_m_match_v(left_edge,right_edge,m_log_match,max);
            this->score_y_match_v(left_edge,right_edge,y_log_match,max);
            this->score_x_match_v(left_edge,right_edge,x_log_match,max);
        }

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            if( left_child_site_to_path_index_p->at( left_edge->get_start_site_index() ) >= last_m_path_index )
            {
                this->score_m_match_v(left_edge,right_edge,m_log_match,max);
                this->score_y_match_v(left_edge,right_edge,y_log_match,max);
                this->score_x_match_v(left_edge,right_edge,x_log_match,max);
            }
        }

        // left site extra edges first
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            if( left_child_site_to_path_index_p->at( left_edge->get_start_site_index() ) >= last_m_path_index )
            {
                this->score_m_match_v(left_edge,right_edge,m_log_match,max);
                this->score_y_match_v(left_edge,right_edge,y_log_match,max);
                this->score_x_match_v(left_edge,right_edge,x_log_match,max);
            }

            while(right_site->has_next_bwd_edge())
            {
                right_edge = right_site->get_next_bwd_edge();

                if( left_child_site_to_path_index_p->at( left_edge->get_start_site_index() ) >= last_m_path_index )
                {
                    this->score_m_match_v(left_edge,right_edge,m_log_match,max);
                    this->score_y_match_v(left_edge,right_edge,y_log_match,max);
                    this->score_x_match_v(left_edge,right_edge,x_log_match,max);
                }
            }
        }

    }
}

/********************************************/


void Reference_alignment::iterate_bwd_edges_for_known_gap(Site * left_site,Site * right_site,vector<Matrix_pointer> *z_slice,vector<Matrix_pointer> *w_slice,
                                                 vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type,bool alignment_end)
{
    if(alignment_end)
    {
        Site *site;
        Edge *edge;

        if(is_x_matrix)
        {
            site = left_site;
        }
        else
        {
            site = right_site;
        }

        edge = site->get_first_bwd_edge();

        this->score_gap_open_v(edge,edge,m_slice,max,is_x_matrix,alignment_end);
        this->score_gap_ext_v(edge,edge,z_slice,max,is_x_matrix,gap_type,alignment_end);

        // then right site extra edges
        //
        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();

            this->score_gap_open_v(edge,edge,m_slice,max,is_x_matrix,alignment_end);
            this->score_gap_ext_v(edge,edge,z_slice,max,is_x_matrix,gap_type,alignment_end);
        }

    }
    else
    {
        if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
        {

            Edge * left_edge = left_site->get_first_bwd_edge();
            Edge * right_edge = right_site->get_first_bwd_edge();


            this->score_gap_double_v(left_edge,right_edge,w_slice,max,is_x_matrix);
            this->score_gap_open_v(left_edge,right_edge,m_slice,max,is_x_matrix,alignment_end);
            this->score_gap_ext_v(left_edge,right_edge,z_slice,max,is_x_matrix,gap_type,alignment_end);


            // then right site extra edges
            //
            while(right_site->has_next_bwd_edge())
            {
                right_edge = right_site->get_next_bwd_edge();
                left_edge = left_site->get_first_bwd_edge();

                this->score_gap_double_v(left_edge,right_edge,w_slice,max,is_x_matrix);
                this->score_gap_open_v(left_edge,right_edge,m_slice,max,is_x_matrix,alignment_end);
                this->score_gap_ext_v(left_edge,right_edge,z_slice,max,is_x_matrix,gap_type,alignment_end);
            }

            // left site extra edges first
            //
            while(left_site->has_next_bwd_edge())
            {
                left_edge = left_site->get_next_bwd_edge();
                right_edge = right_site->get_first_bwd_edge();

                this->score_gap_double_v(left_edge,right_edge,w_slice,max,is_x_matrix);
                this->score_gap_open_v(left_edge,right_edge,m_slice,max,is_x_matrix,alignment_end);
                this->score_gap_ext_v(left_edge,right_edge,z_slice,max,is_x_matrix,gap_type,alignment_end);

                while(right_site->has_next_bwd_edge())
                {
                    right_edge = right_site->get_next_bwd_edge();

                    this->score_gap_double_v(left_edge,right_edge,w_slice,max,is_x_matrix);
                    this->score_gap_open_v(left_edge,right_edge,m_slice,max,is_x_matrix,alignment_end);
                    this->score_gap_ext_v(left_edge,right_edge,z_slice,max,is_x_matrix,gap_type,alignment_end);
                }
            }
        }
    }
}

/********************************************/

void Reference_alignment::iterate_bwd_edges_for_vector_end(Site * left_site,Site * right_site,Matrix_pointer *max,int last_matrix)
{
    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        double best_score = -HUGE_VAL;

        // match score & gap close scores for this match
        //
        double m_log_match = model->log_non_gap();

        this->score_m_match_v(left_edge,right_edge,m_log_match,max);

        best_score = max->score;

        this->score_gap_close_v(left_edge,right_edge,xvectp,max,true);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->y_ind = right->sites_length()-2;
        }

        this->score_gap_close_v(left_edge,right_edge,yvectp,max,false);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->x_ind = left->sites_length()-2;
        }

        // first right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->score_m_match_v(left_edge,right_edge,m_log_match,max);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            this->score_gap_close_v(left_edge,right_edge,xvectp,max,true);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->y_ind = right->sites_length()-2;
            }

            this->score_gap_close_v(left_edge,right_edge,yvectp,max,false);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->x_ind = left->sites_length()-2;
            }
        }

        // left site extra edges then
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->score_m_match_v(left_edge,right_edge,m_log_match,max);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            this->score_gap_close_v(left_edge,right_edge,yvectp,max,false);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->x_ind = left->sites_length()-2;
            }

            this->score_gap_close_v(left_edge,right_edge,xvectp,max,true);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->y_ind = right->sites_length()-2;
            }

            while(right_site->has_next_bwd_edge())
            {
                right_edge = right_site->get_next_bwd_edge();

                this->score_m_match_v(left_edge,right_edge,m_log_match,max);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                }

                this->score_gap_close_v(left_edge,right_edge,xvectp,max,true);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                    max->y_ind = right->sites_length()-2;
                }

                this->score_gap_close_v(left_edge,right_edge,yvectp,max,false);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                    max->x_ind = left->sites_length()-2;
                }
            }
        }

        stringstream ss;
        ss<<"viterbi score: "<<setprecision(8)<<max->score<<" ["<<exp(max->score)<<"] ("<<max->matrix<<")\n";
        Log_output::write_out(ss.str(),3);

    }
}


/********************************************/

void Reference_alignment::score_m_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max, bool allow_any)
{

    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();
    int left_path_index = left_child_site_to_path_index_p->at(left_prev_index);

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();
    int right_path_index = right_child_site_to_path_index_p->at(right_prev_index);

    if( left_path_index == right_path_index )
    {
        double this_score = mvectp->at( left_path_index ).score + m_log_match + left_edge_wght + right_edge_wght;

        if(this->first_is_bigger(this_score,max->score) )
        {
            max->score = this_score;
            max->path_index = left_path_index;
            max->x_ind = left_prev_index;
            max->y_ind = right_prev_index;
            max->x_edge_ind = left_edge->get_index();
            max->y_edge_ind = right_edge->get_index();
            max->matrix = Reference_alignment::m_mat;
        }
    }
}

void Reference_alignment::score_x_match_v(Edge * left_edge,Edge * right_edge,double x_log_match,Matrix_pointer *max)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();
    int left_path_index = left_child_site_to_path_index_p->at(left_prev_index);

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    if( path_to_right_child_site_index_p->at( left_path_index ) == right_prev_index )
    {
        double this_score = xvectp->at(left_path_index).score + x_log_match + left_edge_wght + right_edge_wght;

        if(this->first_is_bigger(this_score,max->score) )
        {
            max->score = this_score;
            max->path_index = left_path_index;
            max->x_ind = left_prev_index;
            max->y_ind = right_prev_index;
            max->x_edge_ind = left_edge->get_index();
            max->y_edge_ind = right_edge->get_index();
            max->matrix = Reference_alignment::x_mat;
        }
    }
}

void Reference_alignment::score_y_match_v(Edge * left_edge,Edge * right_edge,double y_log_match,Matrix_pointer *max)
{

    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();
    int right_path_index = right_child_site_to_path_index_p->at(right_prev_index);

    if( path_to_left_child_site_index_p->at( right_path_index ) == left_prev_index )
    {
        double this_score = yvectp->at(right_path_index).score + y_log_match + left_edge_wght + right_edge_wght;

        if(this->first_is_bigger(this_score,max->score) )
        {
            max->score = this_score;
            max->path_index = right_path_index;
            max->x_ind = left_prev_index;
            max->y_ind = right_prev_index;
            max->x_edge_ind = left_edge->get_index();
            max->y_edge_ind = right_edge->get_index();
            max->matrix = Reference_alignment::y_mat;
        }
    }
}

/************************************************/

void Reference_alignment::score_gap_ext_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix,int gap_type,bool alignment_end)
{

    double edge_wght;
    int path_index;
    Edge *edge;
    int opposite_index;

    if(!alignment_end)
    {
        if(is_x_matrix)
        {
            edge = left_edge;
            path_index = left_child_site_to_path_index_p->at(left_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(left_edge);

            opposite_index = path_to_right_child_site_index_p->at(path_index);
            if(opposite_index != right_edge->get_end_site_index())
            {
                return;
            }
        }
        else
        {
            edge = right_edge;
            path_index = right_child_site_to_path_index_p->at(right_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(right_edge);

            opposite_index = path_to_left_child_site_index_p->at(path_index);
            if(opposite_index != left_edge->get_end_site_index())
            {
                return;
            }
        }
    }
    else
    {
        if(is_x_matrix)
        {
            edge = left_edge;
            path_index = left_child_site_to_path_index_p->at(left_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(left_edge);
        }
        else
        {
            edge = right_edge;
            path_index = right_child_site_to_path_index_p->at(right_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(right_edge);
        }
    }

    double this_score =  (*z_slice)[path_index].score + model->log_gap_ext() + edge_wght;

    // this reduces the terminal and pair-end read gap extension cost.
    //
    if(gap_type != Reference_alignment::normal_gap)
    {
        if(gap_type == Reference_alignment::end_gap)
            this_score =  (*z_slice)[path_index].score + model->log_gap_end_ext() + edge_wght;

        else if(gap_type == Reference_alignment::pair_break_gap)
            this_score =  (*z_slice)[path_index].score + model->log_gap_break_ext() + edge_wght;
    }

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->path_index = path_index;

        if(is_x_matrix)
        {
            max->matrix = Reference_alignment::x_mat;
            max->x_ind = edge->get_start_site_index();
            max->x_edge_ind = edge->get_index();
            if(alignment_end)
            {
                max->y_ind = 0;
            }
        }
        else
        {
            max->matrix = Reference_alignment::y_mat;
            max->y_ind = edge->get_start_site_index();
            max->y_edge_ind = edge->get_index();
            if(alignment_end)
            {
                max->x_ind = 0;
            }
        }
    }
}

void Reference_alignment::score_gap_double_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *w_slice,Matrix_pointer *max,bool is_x_matrix)
{

    double edge_wght;
    int path_index;

    Edge *edge;

    if(is_x_matrix)
    {
        edge = left_edge;
        edge_wght = this->get_log_edge_weight(left_edge);

        path_index = right_child_site_to_path_index_p->at(right_edge->get_end_site_index());
        int opposite_index = path_to_left_child_site_index_p->at(path_index);

        if(opposite_index != left_edge->get_start_site_index())
        {
            return;
        }
    }
    else
    {
        edge = right_edge;
        edge_wght = this->get_log_edge_weight(right_edge);

        path_index = left_child_site_to_path_index_p->at(left_edge->get_end_site_index());
        int opposite_index = path_to_right_child_site_index_p->at(path_index);

        if(opposite_index != right_edge->get_start_site_index())
        {
            return;
        }
    }


    double this_score =  (*w_slice)[path_index].score + model->log_gap_close() + model->log_gap_open() + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->path_index = path_index;

        if(is_x_matrix)
        {
            max->matrix = Reference_alignment::y_mat;
            max->x_ind = left_edge->get_start_site_index();
            max->x_edge_ind = left_edge->get_index();
            max->y_ind = right_edge->get_start_site_index();
            max->y_edge_ind = right_edge->get_index();
        }
        else
        {
            max->matrix = Reference_alignment::x_mat;
            max->x_ind = left_edge->get_start_site_index();
            max->x_edge_ind = left_edge->get_index();
            max->y_ind = right_edge->get_start_site_index();
            max->y_edge_ind = right_edge->get_index();
        }
    }

}

void Reference_alignment::score_gap_open_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix,bool alignment_end)
{

    double edge_wght;
    int path_index;
    Edge *edge;
    int opposite_index;

    if(!alignment_end)
    {
        if(is_x_matrix)
        {
            edge = left_edge;
            path_index = left_child_site_to_path_index_p->at(left_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(left_edge);

            opposite_index = path_to_right_child_site_index_p->at(path_index);
            if(opposite_index != right_edge->get_end_site_index())
            {
                return;
            }
        }
        else
        {
            edge = right_edge;
            path_index = right_child_site_to_path_index_p->at(right_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(right_edge);

            opposite_index = path_to_left_child_site_index_p->at(path_index);
            if(opposite_index != left_edge->get_end_site_index())
            {
                return;
            }
        }

    }
    else
    {
        if(is_x_matrix)
        {
            edge = left_edge;
            path_index = left_child_site_to_path_index_p->at(left_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(left_edge);
        }
        else
        {
            edge = right_edge;
            path_index = right_child_site_to_path_index_p->at(right_edge->get_start_site_index());
            edge_wght = this->get_log_edge_weight(right_edge);
        }
    }

    double this_score =  (*m_slice)[path_index].score + model->log_non_gap()
                + this->get_log_gap_open_penalty(edge->get_start_site_index(),is_x_matrix) + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->matrix = Reference_alignment::m_mat;
        max->path_index = path_index;

        if(is_x_matrix)
        {
            max->x_ind = edge->get_start_site_index();
            max->x_edge_ind = edge->get_index();
            if(alignment_end)
            {
                max->y_ind = 0;
            }
        }
        else
        {
            max->y_ind = edge->get_start_site_index();
            max->y_edge_ind = edge->get_index();
            if(alignment_end)
            {
                max->x_ind = 0;
            }
        }
    }
}

void Reference_alignment::score_gap_close_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix)
{

    double edge_wght;
    int path_index;
    Edge *edge;

    if(is_x_matrix)
    {
        edge = left_edge;
        path_index = left_child_site_to_path_index_p->at(left_edge->get_start_site_index());
        edge_wght = this->get_log_edge_weight(left_edge);

        int opposite_index = path_to_right_child_site_index_p->at(path_index);
        if(opposite_index != right_edge->get_start_site_index())
        {
            return;
        }
    }
    else
    {
        edge = right_edge;
        path_index = right_child_site_to_path_index_p->at(right_edge->get_start_site_index());
        edge_wght = this->get_log_edge_weight(right_edge);

        int opposite_index = path_to_left_child_site_index_p->at(path_index);
        if(opposite_index != left_edge->get_start_site_index())
        {
            return;
        }
    }


    double this_score =  (*z_slice)[path_index].score + this->get_log_gap_close_penalty(edge->get_end_site_index(),is_x_matrix) + edge_wght;


    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->path_index = path_index;

        if(is_x_matrix)
        {
            max->matrix = Reference_alignment::x_mat;
            max->x_ind = edge->get_start_site_index();
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Reference_alignment::y_mat;
            max->y_ind = edge->get_start_site_index();
            max->y_edge_ind = edge->get_index();
        }
    }
}

/********************************************/
