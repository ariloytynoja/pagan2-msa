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

#include "main/basic_alignment.h"
#include "main/sequence.h"
#include "utils/exceptions.h"
#include "main/node.h"
#include <iomanip>
#include <fstream>

using namespace std;
using namespace ppa;

Basic_alignment::Basic_alignment() { }

/********************************************/


void Basic_alignment::build_ancestral_sequence(Sequence *sequence, vector<Path_pointer> *path, bool is_reads_sequence)
{

    // The path is given as input.
    // This will create sites with correct child sites.
    this->create_ancestral_sequence(sequence,path,is_reads_sequence);

    // This will add the edges connecting sites.
    this->create_ancestral_edges(sequence);

    // This will do a check-up and delete edges if needed:
    // " To mimic PRANK+F: one has to scan the sequence again to find boundaries 'Match/Skipped' and 'Skipped/Matched'
    //   and record them on the edges. If they are above a limit, the edges (and all between them) are deleted. "
    this->check_skipped_boundaries(sequence);


    if(Settings_handle::st.noise>=4)
    {
        Log_output::write_out("ANCESTRAL SEQUENCE:\n"+sequence->print_sequence(),4);
        Log_output::write_out(this->print_path(path),6);
    }

    sequence->is_read_sequence(is_reads_sequence);
}

void Basic_alignment::create_ancestral_sequence(Sequence *sequence, vector<Path_pointer> *path, bool is_reads_sequence)
{

    vector<Edge> *edges = sequence->get_edges();

    Site first_site( edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_children(0,0);
    first_site.set_posterior_support(1.0);

    sequence->push_back_site(first_site);

    int l_pos = 1;
    int r_pos = 1;

    bool is_dna = sequence->get_data_type() == Model_factory::dna;

    for(unsigned int i=0;i<path->size();i++)
    {

        Site site( edges );
        site.set_empty_children();
        site.set_posterior_support(path->at(i).mp.full_score);

        // mark the pair-end read1 end so that the edge connecting the pair can be split
        if( pair_end_reads && ( r_pos == y_read1_length || l_pos == x_read1_length ) )
        {
            site.set_site_type(Site::break_start_site);
        }

        if(path->at(i).mp.matrix == Basic_alignment::x_mat)
        {
            int lc = left->get_site_at(l_pos)->get_state();
            site.set_state( lc );

            if(is_reads_sequence && (Settings_handle::st.is("use-consensus") || Settings_handle::st.is("build-contigs")  ))
                this->compute_site_consensus(&site,left,l_pos,right,-1,is_dna);

            if(path->at(i).real_site)
                site.set_path_state( Site::xgapped );
            else
            {
                site.set_path_state( Site::xskipped );
                site.set_branch_count_since_last_used(
                        left->get_site_at(l_pos)->get_branch_count_since_last_used()+1 );
                site.set_branch_distance_since_last_used(
                        left->get_site_at(l_pos)->get_branch_distance_since_last_used()+left_branch_length );
            }
            site.set_children(l_pos,-1);

            l_pos++;
        }
        else if(path->at(i).mp.matrix == Basic_alignment::y_mat)
        {
            int rc = right->get_site_at(r_pos)->get_state();
            site.set_state( rc );

            if(is_reads_sequence && (Settings_handle::st.is("use-consensus") || Settings_handle::st.is("build-contigs") ))
                this->compute_site_consensus(&site,left,-1,right,r_pos, is_dna);

            if(path->at(i).real_site)
                site.set_path_state( Site::ygapped );
            else
            {
                site.set_path_state( Site::yskipped );
                site.set_branch_count_since_last_used(
                        right->get_site_at(r_pos)->get_branch_count_since_last_used()+1 );
                site.set_branch_distance_since_last_used(
                        right->get_site_at(r_pos)->get_branch_distance_since_last_used()+right_branch_length );
            }
            site.set_children(-1,r_pos);

            r_pos++;
        }
        else if(path->at(i).mp.matrix == Basic_alignment::m_mat)
        {
            int lc = left->get_site_at(l_pos)->get_state();
            int rc = right->get_site_at(r_pos)->get_state();
            site.set_state( model->parsimony_state(lc,rc) );

            int sttmp = site.get_state();
            if(is_reads_sequence && (Settings_handle::st.is("use-consensus") || Settings_handle::st.is("build-contigs") ))
                this->compute_site_consensus(&site,left,l_pos,right,r_pos, is_dna);
//            cout<<"b "<<sttmp<<"; a "<<site.get_state()<<"\n";

            site.set_path_state( Site::matched );

            site.set_children(l_pos,r_pos);

            l_pos++; r_pos++;
        }


        sequence->push_back_site(site);

    }

    Site last_site( edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_children(left->sites_length()-1,right->sites_length()-1);
    last_site.set_posterior_support(1.0);
    sequence->push_back_site(last_site);

}

void Basic_alignment::create_ancestral_edges(Sequence *sequence)
{

    vector<Site> *sites = sequence->get_sites();

    vector<int> left_child_index;
    vector<int> right_child_index;

    // First create an index for the child's sites in the parent
    //
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        Site *lsite;
        Site *rsite;

        if(offspring->left_index>=0)
        {
            lsite = left->get_site_at(offspring->left_index);
            left_child_index.push_back(i);
        }
        if(offspring->right_index>=0)
        {
            rsite = right->get_site_at(offspring->right_index);
            right_child_index.push_back(i);
        }
    }

    if(Settings::noise>4)
    {
        stringstream ss;
        ss<<"Child sequence site indeces:"<<endl;
        for(unsigned int i=0;i<left_child_index.size();i++)
            ss<<left_child_index.at(i)<<" ";

        ss<<endl;
        for(unsigned int i=0;i<right_child_index.size();i++)
            ss<<right_child_index.at(i)<<" ";
        ss<<endl;

        Log_output::write_out(ss.str(),5);
    }

    // Then copy the edges of child sequences in their parent.
    // Additionally, create edges for cases where skipped gap is flanked by a new gap.
    //
    Edge_history prev(-1,-1);

    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *psite = &sites->at(i);
        int pstate = psite->get_path_state();


        Site_children *offspring = psite->get_children();

        // left sequence is matched
        if(offspring->left_index>=0)
        {
            Site *tsite = left->get_site_at(offspring->left_index);

            if( tsite->has_bwd_edge() )
            {
                Edge *child = tsite->get_first_bwd_edge();

                this->transfer_child_edge(sequence, child, &left_child_index, left_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(sequence, child, &left_child_index, left_branch_length );
                }
            }

            // these create edges to/from skipped sites flanked by gaps.
            if( (pstate == Site::matched || pstate == Site::ends_site ) && prev.left_skip_site_index >= 0 && edges_for_skipped_flanked_by_gaps)
            {
                // edge from the skipped site to the *next* site
                // as no better info is available, *this* edge is "copied" to one coming to the current site
                Edge query(prev.left_skip_site_index,prev.left_skip_site_index+1);
                int ind = left->get_fwd_edge_index_at_site(prev.left_skip_site_index,&query);

                if(ind>=0)
                {
                    Edge *child = &left->get_edges()->at(ind);
                    Edge edge( left_child_index.at(prev.left_skip_site_index), i );
                    this->transfer_child_edge(sequence, edge, child, left_branch_length );
                }

                prev.left_skip_site_index = -1;
            }
            else if(pstate == Site::xskipped && ( prev.path_state == Site::xgapped || prev.path_state == Site::ygapped ) && edges_for_skipped_flanked_by_gaps)
            {
                // the same here: use this as a template for an extra edge
                Edge query(offspring->left_index-1, offspring->left_index);
                int ind = left->get_bwd_edge_index_at_site(offspring->left_index,&query);

                if(ind>=0)
                {
                    Edge *child = &left->get_edges()->at(ind);
                    Edge edge( prev.match_site_index, i );
                    this->transfer_child_edge(sequence, edge, child, left_branch_length );
                }

            }

            if((pstate == Site::xgapped || pstate == Site::xskipped ) && (prev.path_state == Site::ygapped || prev.path_state == Site::yskipped )) // && edges_for_skipped_flanked_by_gaps)
            {
                Edge edge(i-1,i,1.0);
                sequence->push_back_edge(edge);

                sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( sequence->get_current_edge_index() );
                sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( sequence->get_current_edge_index() );

            }

            if(pstate == Site::xskipped)
                prev.left_skip_site_index = offspring->left_index;
            else
                prev.left_real_site_index = offspring->left_index;

            if(pstate == Site::matched)
                prev.match_site_index = i;
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right->get_site_at(offspring->right_index);

            if( tsite->has_bwd_edge() )
            {
                Edge *child = tsite->get_first_bwd_edge();
                this->transfer_child_edge(sequence, child, &right_child_index, right_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(sequence, child, &right_child_index, right_branch_length );
                }
            }

            if( (pstate == Site::matched || pstate == Site::ends_site ) && prev.right_skip_site_index >= 0 && edges_for_skipped_flanked_by_gaps)
            {
                Edge query(prev.right_skip_site_index,prev.right_skip_site_index+1);
                int ind = right->get_fwd_edge_index_at_site(prev.right_skip_site_index,&query);

                if(ind>=0)
                {
                    Edge *child = &right->get_edges()->at(ind);
                    Edge edge( right_child_index.at(prev.right_skip_site_index), i );
                    this->transfer_child_edge(sequence,edge, child, right_branch_length );
                }

                prev.right_skip_site_index = -1;
            }
            else if(pstate == Site::yskipped && ( prev.path_state == Site::xgapped || prev.path_state == Site::ygapped ) && edges_for_skipped_flanked_by_gaps)
            {
                Edge query(offspring->right_index-1, offspring->right_index);
                int ind = right->get_bwd_edge_index_at_site(offspring->right_index,&query);

                if(ind>=0)
                {
                    Edge *child = &right->get_edges()->at(ind);
                    Edge edge( prev.match_site_index, i );
                    this->transfer_child_edge(sequence, edge, child, right_branch_length );
                }

            }

            if((pstate == Site::ygapped || pstate == Site::yskipped) && (prev.path_state == Site::xgapped || prev.path_state == Site::xskipped)) // && edges_for_skipped_flanked_by_gaps)
            {
                Edge edge(i-1,i,1.0);
                sequence->push_back_edge(edge);

                sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( sequence->get_current_edge_index() );
                sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( sequence->get_current_edge_index() );
            }

            if(pstate == Site::yskipped)
                prev.right_skip_site_index = offspring->right_index;
            else
                prev.right_real_site_index = offspring->right_index;

        }
        prev.path_state = pstate;
    }
}

void Basic_alignment::check_skipped_boundaries(Sequence *sequence)
{

    vector<Site> *sites = sequence->get_sites();

    // First, find 'Match/Skipped' and 'Skipped/Matched' boundaries and update the counts
    //
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if( tsite->has_bwd_edge() )
        {
            Edge *edge = tsite->get_first_bwd_edge();

            while( tsite->has_next_bwd_edge() )
            {
                Edge *another = tsite->get_next_bwd_edge();
                if( another->get_start_site_index() > edge->get_start_site_index() )
                    edge = another;
            }

            Site *psite = &sites->at( edge->get_start_site_index() );

            if( ( psite->get_path_state()==Site::matched || psite->get_path_state()==Site::start_site )
                && ( tsite->get_path_state()==Site::xskipped || tsite->get_path_state()==Site::yskipped ) )
            {
                edge->increase_branch_count_as_skipped_edge();
            }
        }

        if( tsite->has_fwd_edge() )
        {
            Edge *edge = tsite->get_first_fwd_edge();

            while( tsite->has_next_fwd_edge() )
            {
                Edge *another = tsite->get_next_fwd_edge();
                if( another->get_start_site_index() < edge->get_start_site_index() )
                    edge = another;
            }

            Site *nsite = &sites->at( edge->get_end_site_index() );

            if( ( tsite->get_path_state()==Site::xskipped || tsite->get_path_state()==Site::yskipped )
                && ( nsite->get_path_state()==Site::matched || nsite->get_path_state()==Site::ends_site ) )
            {
                edge->increase_branch_count_as_skipped_edge();
            }
        }
    }

    // Then, see if any pair of boundaries (covering a skipped gap) is above the limit. Delete the range.
    //
    bool non_skipped = true;
    int skip_start = -1;
    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);
        int tstate = tsite->get_path_state();

        if( non_skipped && ( tstate == Site::xskipped || tstate == Site::yskipped ) )
        {

            if( tsite->has_bwd_edge() )
            {
                Edge *edge = tsite->get_first_bwd_edge();
                while( tsite->has_next_bwd_edge() )
                {
                    Edge *another = tsite->get_next_bwd_edge();
                    if( another->get_start_site_index() > edge->get_start_site_index() )
                        edge = another;
                }
                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                {
                    skip_start = i;
                }
            }

            non_skipped = false;
        }

        if(!non_skipped && skip_start>=0 && tstate == Site::matched)
        {

            int edge_ind = -1;
            if( tsite->has_bwd_edge() )
            {
                Edge *edge = tsite->get_first_bwd_edge();

                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                    edge_ind = edge->get_index();

                while( tsite->has_next_bwd_edge() )
                {
                    edge = tsite->get_next_bwd_edge();

                    if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                        edge_ind = edge->get_index();

                }
            }

            if(edge_ind>=0)
            {
                Log_output::write_out("Basic_alignemnmt: delete range: "+Log_output::itos(edge_ind)+" "+Log_output::itos(skip_start)+" "+Log_output::itos(i)+"\n",4);
                this->delete_edge_range(sequence,edge_ind,skip_start);
            }

            non_skipped = true;
            skip_start = -1;
        }

        if(tstate == Site::xgapped || tstate == Site::ygapped || tstate == Site::matched)
        {
            non_skipped = true;
            skip_start = -1;
        }
    }
}

void Basic_alignment::delete_edge_range(Sequence *sequence,int edge_ind,int skip_start_site)
{
    vector<Edge> *edges = sequence->get_edges();

    Edge *edge = &edges->at(edge_ind);

    int this_site_index = edge->get_start_site_index();

    // if not the start of the range, delete further
    while(this_site_index >= skip_start_site)
    {
        sequence->get_site_at(this_site_index)->set_site_type(Site::non_real);
        sequence->delete_all_bwd_edges_at_site(this_site_index);
        sequence->delete_all_fwd_edges_at_site(this_site_index);
        --this_site_index;
    }

}

void Basic_alignment::transfer_child_edge(Sequence *sequence,Edge *child, vector<int> *child_index, float branch_length,
                                           bool adjust_posterior_weight, float branch_weight)
{
    float edge_weight = 1.0;
    if(weight_edges)
    {
        float weight1 = sequence->get_site_at( child_index->at( child->get_start_site_index() ) )->get_posterior_support();
        float weight2 = sequence->get_site_at( child_index->at( child->get_end_site_index() ) )->get_posterior_support();
        weight1 = this->get_transformed_edge_weight( weight1 );
        weight2 = this->get_transformed_edge_weight( weight2 );

        edge_weight =  weight1 * weight2;
    }

    Edge edge( child_index->at( child->get_start_site_index() ), child_index->at( child->get_end_site_index() ), edge_weight );

    if( reduced_terminal_gap_penalties )
    {
        if( sequence->get_site_at( edge.get_start_site_index() )->get_site_type() == Site::start_site &&
            edge.get_end_site_index() - edge.get_start_site_index() > 1 )
        {
            if(child->get_end_site_index()-child->get_start_site_index()==1)
                edge.set_start_site_index(edge.get_end_site_index()-1);
        }

        if( sequence->get_site_at( edge.get_end_site_index() )->get_site_type() == Site::stop_site &&
            edge.get_end_site_index() - edge.get_start_site_index() > 1 )
        {
            if(child->get_end_site_index()-child->get_start_site_index()==1)
                edge.set_end_site_index(edge.get_start_site_index()+1);
        }
    }

    if( pair_end_reads )
    {
        if( sequence->get_site_at( edge.get_start_site_index() )->get_site_type() == Site::break_start_site &&
            edge.get_end_site_index() - edge.get_start_site_index() > 1 )
        {
            // remove the mark
            //
            sequence->get_site_at( edge.get_start_site_index() )->set_site_type(Site::real_site);

            // split the edge to two
            //
            edge.set_end_site_index(edge.get_start_site_index()+1);

            this->transfer_child_edge(sequence, edge, child, branch_length, adjust_posterior_weight, branch_weight);

            edge.set_start_site_index( child_index->at( child->get_end_site_index() )-1 );
            edge.set_end_site_index( child_index->at( child->get_end_site_index() ) );

            this->transfer_child_edge(sequence, edge, child, branch_length, adjust_posterior_weight, branch_weight);

            return;
        }
    }


    this->transfer_child_edge(sequence, edge, child, branch_length, adjust_posterior_weight, branch_weight);
}


void Basic_alignment::transfer_child_edge(Sequence *sequence, Edge edge, Edge *child, float branch_length,
                                           bool adjust_posterior_weight, float branch_weight)
{

//    cout<<edge.get_start_site_index()<<" "<<edge.get_end_site_index()<<": used "<<child->is_used()<<", since last used "<<child->get_branch_count_since_last_used()<<endl;

    // No identical copies
    if(sequence->get_site_at( edge.get_end_site_index() )->contains_bwd_edge( &edge ) )
    {
        sequence->get_site_at( edge.get_end_site_index() )->update_bwd_edge_details( &edge );
        return;
    }

    // Limits for copying old edges:
    //  first, number of nodes since last used
    if(!child->is_used() && child->get_branch_count_since_last_used()+1 > max_allowed_skip_branches)
        return;

    //  then, total branch distance since last used
    if(!child->is_used() && child->get_branch_distance_since_last_used()+branch_length > max_allowed_skip_distance)
        return;

    // Comparison of distance and node count since last used to find boundaries of path branches.
    // Only start and end of an alternative path should be penalised; continuation on a path not.
    //
    float dist_start = sequence->get_site_at(edge.get_start_site_index())->get_branch_distance_since_last_used();
    float dist_end   = sequence->get_site_at(edge.get_end_site_index()  )->get_branch_distance_since_last_used();

    int count_start = sequence->get_site_at(edge.get_start_site_index())->get_branch_count_since_last_used();
    int count_end   = sequence->get_site_at(edge.get_end_site_index()  )->get_branch_count_since_last_used();

    // Sites on the two ends of an edge have different history: branch point that should be penalised
    if( dist_start != dist_end || count_start != count_end )
    {
        edge.set_branch_distance_since_last_used( max(dist_start,dist_end) );
        edge.set_branch_count_since_last_used( max(count_start,count_end) );

        if(adjust_posterior_weight)
            if(weighted_branch_skip_penalty)
                edge.multiply_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_weight * ( 1.0 - exp( -1.0 * branch_length ) ) );
            else
                edge.multiply_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_probability );
        else
            edge.multiply_weight(child->get_posterior_weight());

    }
    // Edge is not used: just update the history
    else if(!child->is_used() && count_start == 0 &&  count_end == 0)
    {
        edge.set_branch_distance_since_last_used( child->get_branch_distance_since_last_used()+branch_length );
        edge.set_branch_count_since_last_used( child->get_branch_count_since_last_used()+1 );

        if(adjust_posterior_weight)
            if(weighted_branch_skip_penalty)
                edge.multiply_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_weight * ( 1.0 - exp( -1.0 * branch_length ) ) );
            else
                edge.multiply_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_probability );
        else
            edge.multiply_weight(child->get_posterior_weight());

    }
    else if(!child->is_used())
    {
        edge.set_branch_distance_since_last_used( child->get_branch_distance_since_last_used()+branch_length );
        edge.set_branch_count_since_last_used( child->get_branch_count_since_last_used()+1 );
    }

    if(!sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
    {

// ->
        if(!child->is_used())
            edge.set_branch_count_as_skipped_edge( child->get_branch_count_as_skipped_edge() );
        else
            edge.set_branch_count_as_skipped_edge( 0 );
// <-
        sequence->push_back_edge(edge);

        sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( sequence->get_current_edge_index() );
        sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( sequence->get_current_edge_index() );
    }
}

/********************************************/

