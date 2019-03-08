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

#include "main/viterbi_alignment.h"
#include "main/sequence.h"
#include "utils/exceptions.h"
#include "utils/codon_translation.h"
#include "utils/find_anchors.h"
#include "utils/exonerate_queries.h"
#include "utils/tunnel_matrix.h"
#include "main/node.h"
#include <iomanip>
#include <fstream>
#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string.hpp>

#ifdef NCBI_TOOLKIT
#include "utils/ncbi_blast.h"
#endif

using namespace std;
using namespace ppa;


Viterbi_alignment::Viterbi_alignment() { tunnel_defined = false; reduced_terminal_gap_penalties = false; }

float Viterbi_alignment::define_tunnel(Sequence *left_sequence,Sequence *right_sequence, string left_sequence_id, string right_sequence_id, Evol_model *evol_model,bool compute_coverage)
{

    string s1 = left_sequence->get_sequence_string(false);
    string s2 = right_sequence->get_sequence_string(false);

    Codon_translation ct;

    bool is_dna = evol_model->get_data_type() == Model_factory::dna;

    if( ( evol_model->get_data_type() == Model_factory::dna && Settings_handle::st.is("codons") ) || evol_model->get_data_type() == Model_factory::codon)
    {
        ct.define_translation_tables();
        s1 = ct.gapped_DNA_to_protein(&s1);
        s2 = ct.gapped_DNA_to_protein(&s2);
        is_dna = false;
    }

    vector<Substring_hit> hits;
    int p_len = Settings_handle::st.get("prefix-hit-length").as<int>();

//        struct timespec tcpu_start, tcpu_finish;
//        clock_gettime(CLOCK_MONOTONIC, &tcpu_start);

    Find_anchors fa;
    if(Settings_handle::st.is("use-prefix-anchors"))
    {
        fa.find_long_substrings(&s1,&s2,&hits,p_len);
    }
    else if(Settings_handle::st.is("hmmer-anchors"))
    {
        fa.find_hmmer_anchors(&s1,&s2,&hits);
    }
#ifdef NCBI_TOOLKIT
    else if(Settings_handle::st.is("ncbi"))
    {
        BankBlaster blaster;
        Blast_options opt;
        mol_type mol;

        if(is_dna)
        {
            mol = dna;
        }else{
            mol = aminoa;
        }

        if( ( Settings_handle::st.is("both-strands") ) && is_dna )
        {
            opt.strand_opt = strand_both;
        }

        opt.wordsize = Settings_handle::st.get("blast-wordsize").as<int>();
        if(mol != aminoa){
            opt.match_reward = Settings_handle::st.get("blast-match-reward").as<int>();
            opt.mismatch_penalty = Settings_handle::st.get("blast-mismatch-penalty").as<int>();
        }else{
            opt.word_threshold_score = Settings_handle::st.get("blast-word-threshold").as<int>();

            string scoring_matrix = boost::to_upper_copy(Settings_handle::st.get("blast-scoring-matrix").as<string>());

            if(scoring_matrix == "BLOSUM45"){
                opt.scoring_matrix = blosum45;
            }else if(scoring_matrix == "BLOSUM50"){
                opt.scoring_matrix = blosum50;
            }else if(scoring_matrix == "BLOSUM62"){
                opt.scoring_matrix = blosum62;
            }else if(scoring_matrix == "BLOSUM80"){
                opt.scoring_matrix = blosum80;
            }else if(scoring_matrix == "BLOSUM90"){
                opt.scoring_matrix = blosum90;
            }else{
                Log_output::write_warning("Invalid scoring matrix given. Using default (BLOSUM62) instead.",0);
                opt.scoring_matrix = blosum62;
            }
        }

        blaster.runSingleBlast(s1, s2, left_sequence_id, right_sequence_id, hits, mol, &opt);

    }
#endif
    else
    {
        Exonerate_queries er;
        if(er.test_executable())
            er.local_pairwise_alignment(&s1,&s2,&hits);
    }

//        clock_gettime(CLOCK_MONOTONIC, &tcpu_finish);
//        double elapsed;
//        elapsed = (tcpu_finish.tv_sec - tcpu_start.tv_sec);
//        elapsed += (tcpu_finish.tv_nsec - tcpu_start.tv_nsec) / 1000000000.0;
//        cout<<"\n\ntime "<<elapsed<<endl<<endl;

    s1 = left_sequence->get_sequence_string(true);
    s2 = right_sequence->get_sequence_string(true);

    if( ( evol_model->get_data_type() == Model_factory::dna && Settings_handle::st.is("codons") ) || evol_model->get_data_type() == Model_factory::codon)
    {
        s1 = ct.gapped_DNA_to_protein(&s1);
        s2 = ct.gapped_DNA_to_protein(&s2);
    }

#ifdef NCBI_TOOLKIT
    if(Settings_handle::st.is("ncbi"))
    {

        int overlap_total = Settings_handle::st.get("ncbi-threshold-overlap-total").as<int>();
        int overlap_partly = Settings_handle::st.get("ncbi-threshold-overlap-partly").as<int>();

        fa.eliminate_bad_hits(hits, overlap_total, overlap_partly);
        fa.define_tunnel_with_overlapping_hits(hits, upper_bound, lower_bound, s1, s2, Settings_handle::st.get("anchors-offset").as<int>(), empty_tunnel_blocks);

    }
#endif
    else
    {
        fa.check_hits_order_conflict(&s1,&s2,&hits);

        fa.define_tunnel(&hits,&upper_bound,&lower_bound,&s1,&s2);
    }

    float coverage = 0;

    if(compute_coverage)
    {
        long long int sum = 0;
        for(int i=0;i<(int)s1.length();i++)
            sum += lower_bound.at(i)-upper_bound.at(i);

        coverage = ((float)sum/((float)s1.length()*(float)s2.length()));

        stringstream s;
        s.precision(4);
        s<<"Anchoring: Computing "<<coverage*100.0<<"% of DP matrix.\n";
        Log_output::write_out(s.str(),2);
    }

    tunnel_defined = true;

    return coverage;
}

void Viterbi_alignment::align(Sequence *left_sequence,Sequence *right_sequence,Evol_model *evol_model,
                              float l_branch_length,float r_branch_length,bool is_reads_sequence)
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

    if( is_reads_sequence || Settings_handle::st.is("keep-all-edges") )
        this->set_reads_alignment_settings();

    this->set_additional_settings();

    // mark sites where no gap open/close penalty. Important for paired reads
    if(reduced_terminal_gap_penalties)
    {
        this->mark_no_gap_penalty_sites(left, right);
    }


    // Set the edge weighting scheme, define the dynamic-programming matrices
    //
    log_edge_weight = &ppa::Viterbi_alignment::edge_log_posterior_weight;
    edge_weight = &ppa::Viterbi_alignment::edge_posterior_weight;

    transform_edge_weight = &ppa::Viterbi_alignment::square_root_edge_weight;
    if( Settings_handle::st.is("no-weight-transform") )
        transform_edge_weight = &ppa::Viterbi_alignment::plain_edge_weight;
    if( Settings_handle::st.is("cuberoot-weight-transform") )
        transform_edge_weight = &ppa::Viterbi_alignment::cube_root_edge_weight;

    int left_length = left->sites_length();    int right_length = right->sites_length();

    Log_output::write_out("Viterbi_alignment: lengths: "+this->itos(left_length)+" "+this->itos(right_length),3);

    /*
    align_array M(boost::extents[left_length-1][right_length-1]);
    align_array X(boost::extents[left_length-1][right_length-1]);
    align_array Y(boost::extents[left_length-1][right_length-1]);
    */
    Matrix_pointer empty;
    empty.bwd_score = -HUGE_VAL;
    empty.full_score = -HUGE_VAL;
    empty.fwd_score = -HUGE_VAL;

    align_array M(left_length-1, right_length-1, upper_bound, lower_bound, &empty);
    align_array X(left_length-1, right_length-1, upper_bound, lower_bound, &empty);
    align_array Y(left_length-1, right_length-1, upper_bound, lower_bound, &empty);

    match = &M;    xgap = &X;    ygap = &Y;

    Log_output::write_out("Viterbi_alignment: matrix created",3);

    // Dynamic programming loop
    // Note that fwd full probability is integrated in
    // the main dp loop but bwd comutation is done separately
    //
    this->initialise_array_corner();

    int j_max = match->getSize_y();
    int i_max = match->getSize_x();

    if( tunnel_defined )
    {
        for(int i=0;i<i_max;i++)
        {
            for(int j=0;j<j_max;j++)
            {
                if(j>=upper_bound[i] && j<=lower_bound[i])
                {
                    this->compute_fwd_scores(i,j);
                }
            }
        }
    }
    else
    {
        for(int j=0;j<j_max;j++)
        {
            for(int i=0;i<i_max;i++)
            {
                this->compute_fwd_scores(i,j);
            }
        }
    }
    Log_output::write_out("Viterbi_alignment: matrix filled",3);



    // Find the incoming edge in the end corner; also, get the full_fwd probability
    //
    Site * left_site = left->get_site_at(i_max);
    Site * right_site = right->get_site_at(j_max);


    Matrix_pointer max_end;
    this->iterate_bwd_edges_for_end_corner(left_site,right_site,&max_end);
    max_end.bwd_score = 1.0;
    max_end.full_score = 1.0;

    if(max_end.score==-HUGE_VAL && tunnel_defined)
    {
        Log_output::write_msg("anchored alignment failed: trying again",1);

        for(int j=0;j<j_max;j++)
        {
            for(int i=0;i<i_max;i++)
            {
                this->compute_fwd_scores(i,j);
            }
        }

        left_site = left->get_site_at(i_max);
        right_site = right->get_site_at(j_max);

        this->iterate_bwd_edges_for_end_corner(left_site,right_site,&max_end);
        max_end.bwd_score = 1.0;
        max_end.full_score = 1.0;

    }

    if(max_end.score==-HUGE_VAL)
    {
        Log_output::write_out("\nViterbi_alignment: max_end.score==-HUGE_VAL\nleft\n"+this->left->print_sequence()+
        "\nright\n"+this->right->print_sequence()+"\nmatrix\n"+this->print_matrices(),1);
    }
    Log_output::write_out("Viterbi_alignment: corner found",3);


    // If needed, compute the bwd posterior probabilities
    //
    if(compute_full_score)
    {
        // Backward dynamic-programming loop
        //
        this->initialise_array_corner_bwd();

        for(int j=j_max-1;j>=0;j--)
        {
            for(int i=i_max-1;i>=0;i--)
            {
                this->compute_bwd_full_score(i,j);
            }
        }


        // Check that full probability computation worked
        //
        Matrix_pointer max_start = (*match)[0][0];

        Log_output::write_out(" bwd full probability: "+this->ftos(log(max_start.bwd_score))
                        +" ["+this->ftos(max_start.bwd_score)+"]",3);

        double ratio = max_end.fwd_score/max_start.bwd_score;

        if(ratio<0.99 || ratio>1.01)
            Log_output::write_out("Problem in computation? fwd: "+this->ftos(max_end.fwd_score)
                            +", bwd: "+this->ftos(max_start.bwd_score),1);


        // Compute psoterior probbailities
        //
        for(int j=0;j<j_max;j++)
        {
            for(int i=0;i<i_max;i++)
            {
                this->compute_posterior_score(i,j,max_end.fwd_score);
            }
        }

        if(Settings::noise>5)
            Log_output::write_out(this->print_matrices(),5);

    }



    // It is more convenient to build the sequence forward;
    // thus, backtrace the path into a vector ('path') and use that in the next step
    //
    // Backtrack the Viterbi path
    if( ! Settings_handle::st.is("sample-path") )
    {
        Path_pointer pp(max_end,true);
        this->backtrack_new_path(&path,pp);

        Log_output::write_out("Viterbi_alignment: path found",3);


        // Now build the sequence forward following the path saved in a vector;
        //
        ancestral_sequence = new Sequence(path.size(),model->get_data_type());
        this->build_ancestral_sequence(ancestral_sequence,&path,is_reads_sequence);

        Log_output::write_out("Viterbi_alignment: sequence built",3);
    }

    // Instead of Viterbi path, sample a path from the posterior probabilities
    else
    {
        Matrix_pointer sample_end;
        this->iterate_bwd_edges_for_sampled_end_corner(left_site,right_site,&sample_end);
        sample_end.bwd_score = 1.0;
        sample_end.full_score = 1.0;

        Path_pointer sp(sample_end,true);
        vector<Path_pointer> sample_path;

        this->sample_new_path(&sample_path,sp);

        Log_output::write_out("Viterbi_alignment: path sampled",3);


        // Now build the sequence forward following the path saved in a vector;
        //
        ancestral_sequence = new Sequence(sample_path.size(),model->get_data_type());
        this->build_ancestral_sequence(ancestral_sequence,&sample_path,is_reads_sequence);

        Log_output::write_out("Viterbi_alignment: sequence sampled and built",3);
    }


    // Sample additional paths [not finished!]
    //
    if( Settings_handle::st.get("sample-additional-paths").as<int>() > 0 )
    {
        int iter = Settings_handle::st.get("sample-additional-paths").as<int>();
        for(int i=0;i<iter;i++)
        {
            Matrix_pointer sample_end;
            this->iterate_bwd_edges_for_sampled_end_corner(left_site,right_site,&sample_end);
            sample_end.bwd_score = 1.0;
            sample_end.full_score = 1.0;

            Path_pointer sp(sample_end,true);
            vector<Path_pointer> sample_path;

            this->sample_new_path(&sample_path,sp);

            Log_output::write_out("Viterbi_alignment: additional path sampled",3);

            // Now build the sequence forward following the path saved in a vector;
            //
            Sequence *sampled_sequence = new Sequence(sample_path.size(),model->get_data_type());
            this->build_ancestral_sequence(sampled_sequence,&sample_path,is_reads_sequence);

            Log_output::write_out("Viterbi_alignment: additional sequence sampled and built",3);
            Log_output::write_out("Viterbi_alignment: sampled alignment",3);

            this->merge_sampled_sequence(ancestral_sequence, sampled_sequence);

            delete sampled_sequence;
        }
    }


    // Make the posterior probability plots
    //
    if( Settings_handle::st.is("mpost-posterior-plot-file") )
    {

        if(Settings_handle::st.is("plot-slope-up") )
            this->plot_posterior_probabilities_up();
        else
            this->plot_posterior_probabilities_down();
    }
    
}

bool Viterbi_alignment::replace_largest_tunnel_block_with_gap_tunnel()
{
    int tunnel_end = (int)lower_bound.size() - 1;
    int remove_threshold = Settings_handle::st.get("force-gap-threshold").as<int>();
    int tunnel_width = Settings_handle::st.get("anchors-offset").as<int>();
    bool wide_tunnel = Settings_handle::st.is("force-gap-wide-tunnel");

    /*
    std::cout << "--- EMPTY BLOCKS ---" << endl;
    for(int i = 0; i < (int)empty_tunnel_blocks.size(); i++){
        empty_tunnel_blocks[i].print();
    }*/

    if(empty_tunnel_blocks.size() < 1 || empty_tunnel_blocks.back().size() < remove_threshold){
        stringstream ss;
        ss << "No large enough blocks to replace." << endl;
        Log_output::write_out(ss.str(),2);
        return false;
    }

    Tunnel_block& largest_block = empty_tunnel_blocks.back();

    if(wide_tunnel){
        for(int i = largest_block.start.x; i < largest_block.end.x - tunnel_width; i++){
            lower_bound.at(i) = largest_block.start.y + tunnel_width;
        }
        for(int i = largest_block.start.x-1; i >= 0; i--){
            if(lower_bound.at(i) > lower_bound.at(i+1)){
                lower_bound.at(i) = lower_bound.at(i+1);
            }
            else{
                break;
            }
        }

    }else{
        int a = 0;
        for(int i = largest_block.start.x; i < largest_block.end.x; i++){
            lower_bound.at(i) = largest_block.start.y;
            upper_bound.at(i) = std::min(largest_block.start.y, upper_bound.at(i)+a);
            a++;
        }

        upper_bound.at(largest_block.end.x) = largest_block.start.y;


        for(int i = largest_block.start.x-1; i >= 0; i--){
            if(lower_bound.at(i) > lower_bound.at(i+1)){
                lower_bound.at(i) = lower_bound.at(i+1);
            }
            else{
                break;
            }
        }

        int last_i = std::min(largest_block.end.x+tunnel_width+1, tunnel_end);
        int b = 0;
        for(int i = last_i; i >= largest_block.end.x+1; i--){
            upper_bound.at(i) = std::max(upper_bound.at(last_i)-b, largest_block.start.y);
            b++;
            /*
            if(upper_bound.at(i) < largest_block.end.y){
                upper_bound.at(i) = largest_block.end.y;
            }
            else{
                break;
            }*/
        }
    }

    stringstream ss;
    ss << "Largest block (" << largest_block.end.x - largest_block.start.x << "x" << largest_block.end.y - largest_block.start.y << ") replaced with a gap." << endl;
    Log_output::write_out(ss.str(),2);


    empty_tunnel_blocks.pop_back();


    /*
    vector<pair<int, int> > v;
    Find_anchors fa;
    fa.plotRFile(v, upper_bound, lower_bound, upper_bound.size()-1, lower_bound[lower_bound.size()-1], "remove_block.R");
    */

    return true;

}

long long int Viterbi_alignment::get_predicted_memory_consumption(int sequence_1_length, int sequence_2_length)
{
    if(tunnel_defined == false)
    {
        return (long long int)sequence_1_length * (long long int)sequence_2_length * 64 * 3 + 50000000;
    }

    long long int sum = 0;
    for(int i=0;i<(int)upper_bound.size();i++)
        sum += lower_bound.at(i)-upper_bound.at(i);

    return sum*65*3 + 50000000; //cell size at maxtrix = ~64 bytes (3 matrixes). 50MB reserve for other memory allocation.

}

/********************************************/

void Viterbi_alignment::merge_sampled_sequence(Sequence *ancestral_sequence, Sequence *sampled_sequence)
{
    ancestral_sequence->initialise_unique_index();
    sampled_sequence->initialise_unique_index();

    if(Settings::noise>4)
    {
        stringstream ss;
        ss<<"\nancestral\n"<<ancestral_sequence->print_path();
        ss<<"\nsampled\n"<<sampled_sequence->print_path();
        ss<<"\nancestral\n"<<ancestral_sequence->print_sequence();
        ss<<"\nsampled\n"<<sampled_sequence->print_sequence();
        Log_output::write_out(ss.str(),5);
    }

    vector<int> sample_index_for_added;
    vector<int> sample_index_in_ancestral;

    vector<Unique_index> *index = sampled_sequence->get_unique_index();
    for(int i=0;i<(int) index->size();i++)
    {
        int anc_index = ancestral_sequence->unique_index_of_term( &index->at(i) );

        if( anc_index>=0 )
        {
            sample_index_in_ancestral.push_back(anc_index);
        }
        else
        {
            // Make a copy of the site in sampled sequence that is missing from
            // ancestral sequence and add it there; can't make a carbon copy as
            // indeces (neighbour sites, edges) aren't correct
            //
            Site *sample_site = sampled_sequence->get_site_at( index->at(i).site_index );

            Site site( ancestral_sequence->get_edges() );
            ancestral_sequence->copy_site_details(sample_site, &site);
            ancestral_sequence->push_back_site( site );

            ancestral_sequence->add_term_in_unique_index( index->at(i) );

            // Keep track of sites added so that edges can be added below
            //
            sample_index_in_ancestral.push_back( ancestral_sequence->get_current_site_index() );
            sample_index_for_added.push_back( sample_index_in_ancestral.size()-1 );
        }

    }
if( sample_index_for_added.size()>0 )

    // If sites were added, add also their edges
    //
    for(int i=0;i<(int) sample_index_for_added.size();i++)
    {
        Site *sample_site = sampled_sequence->get_site_at( sample_index_for_added.at( i ) );

        if( sample_site->has_bwd_edge() )
        {
            Edge *sample_edge = sample_site->get_first_bwd_edge();

            int edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
            int edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

            Edge edge( edge_start_site, edge_end_site );

            if(!ancestral_sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
            {

                ancestral_sequence->copy_edge_details( sample_edge, &edge );
                ancestral_sequence->push_back_edge(edge);

                ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
            }

            while( sample_site->has_next_bwd_edge() )
            {
                sample_edge = sample_site->get_next_bwd_edge();

                edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
                edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

                Edge edge( edge_start_site, edge_end_site );

                if(!ancestral_sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
                {

                    ancestral_sequence->copy_edge_details( sample_edge, &edge );
                    ancestral_sequence->push_back_edge(edge);

                    ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                    ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
                }
            }
        }

        // The same for fwd edges
        if( sample_site->has_fwd_edge() )
        {
            Edge *sample_edge = sample_site->get_first_fwd_edge();

            int edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
            int edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

            Edge edge( edge_start_site, edge_end_site );

            if(!ancestral_sequence->contains_this_fwd_edge_at_site(edge.get_start_site_index(),&edge))
            {
                ancestral_sequence->copy_edge_details( sample_edge, &edge );
                ancestral_sequence->push_back_edge(edge);

                ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
            }

            while( sample_site->has_next_fwd_edge() )
            {
                sample_edge = sample_site->get_next_fwd_edge();

                edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
                edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

                Edge edge( edge_start_site, edge_end_site );

                if(!ancestral_sequence->contains_this_fwd_edge_at_site(edge.get_start_site_index(),&edge))
                {
                    ancestral_sequence->copy_edge_details( sample_edge, &edge );
                    ancestral_sequence->push_back_edge(edge);

                    ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                    ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
                }
            }
        }
    }

    if( sample_index_for_added.size()>0 )
    {
        // Now reorder the ancestral sequence

        ancestral_sequence->sort_sites_vector();

        ancestral_sequence->remap_edges_vector();

        if(Settings::noise>4)
        {
            Log_output::write_out("\nmerged ancestral\n"+ancestral_sequence->print_sequence(),5);
        }
    }
}

/********************************************/

void Viterbi_alignment::initialise_array_corner()
{
    double small = -HUGE_VAL;

    ((*match)[0][0]).score = 0.0;
    ((*match)[0][0]).fwd_score = 1.0;
    ((*xgap)[0][0]).score = small;
    ((*ygap)[0][0]).score = small;


    /*MISSING FEATURE: allow extending existing gap (for fragments) */
}

/********************************************/

void Viterbi_alignment::initialise_array_corner_bwd()
{
    int i_max = match->shape()[0];
    int j_max = match->shape()[1];

    (*match)[i_max-1][j_max-1].bwd_score = model->non_gap();

    Site * left_site = left->get_site_at(i_max);
    Site * right_site = right->get_site_at(j_max);

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {
        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        double left_edge_wght = this->get_edge_weight(left_edge);
        int left_index = left_edge->get_start_site_index();

        double right_edge_wght = this->get_edge_weight(right_edge);
        int right_index = right_edge->get_start_site_index();

        (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                right_edge_wght = this->get_edge_weight(right_edge);
                right_index = right_edge->get_start_site_index();

                (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;
            }
        }

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

            while(left_site->has_next_bwd_edge())
            {
                left_edge = left_site->get_next_bwd_edge();

                left_edge_wght = this->get_edge_weight(left_edge);
                left_index = left_edge->get_start_site_index();

                (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;
            }
        }
    }

    if( left_site->has_bwd_edge() )
    {
        Edge * left_edge = left_site->get_first_bwd_edge();

        double left_edge_wght = this->get_edge_weight(left_edge);
        int left_index = left_edge->get_start_site_index();

        (*xgap)[left_index][j_max-1].bwd_score = model->gap_close();// * left_edge_wght;

        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            (*xgap)[left_index][j_max-1].bwd_score = model->gap_close();// * left_edge_wght;
        }
    }

    if(right_site->has_bwd_edge())
    {
        Edge * right_edge = right_site->get_first_bwd_edge();

        double right_edge_wght = this->get_edge_weight(right_edge);
        int right_index = right_edge->get_start_site_index();

        (*ygap)[i_max-1][right_index].bwd_score = model->gap_close();// * right_edge_wght;

        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*ygap)[i_max-1][right_index].bwd_score = model->gap_close();// * right_edge_wght;
        }
    }
}

void Viterbi_alignment::compute_fwd_scores(int i,int j)
{
    if(i==0 && j==0)
        return;

    int j_gap_type = Viterbi_alignment::normal_gap;
    int i_gap_type = Viterbi_alignment::normal_gap;

    if( j==0 || j == (int) match->shape()[1]-1 )
    {
        if(!Settings_handle::st.is("no-terminal-edges"))
            j_gap_type = Viterbi_alignment::end_gap;
    }

    if(pair_end_reads && j==y_read1_length)
    {
        j_gap_type = Viterbi_alignment::pair_break_gap;
    }

    if( i==0 || i == (int) match->shape()[0]-1 )
    {
        if(!Settings_handle::st.is("no-terminal-edges"))
            i_gap_type = Viterbi_alignment::end_gap;
    }

    if(pair_end_reads && i==x_read1_length)
    {
        i_gap_type = Viterbi_alignment::pair_break_gap;
    }


    int left_index = i;
    int right_index = j;

    Site * left_site;
    Site * right_site;

    Matrix_pointer *max_x = &(*xgap)[i][j];
    Matrix_pointer *max_y = &(*ygap)[i][j];
    Matrix_pointer *max_m = &(*match)[i][j];


    if(left_index>0)
    {
        left_site = left->get_site_at(left_index);

        /*
        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][j] ];
        align_slice y_slice = (*ygap)[ indices[ range( 0,ygap->shape()[0] ) ][j] ];
        align_slice m_slice = (*match)[ indices[ range( 0,match->shape()[0] ) ][j] ];
        */

        align_slice& x_slice = (*xgap).slice_y(j);
        align_slice& y_slice = (*ygap).slice_y(j);
        align_slice& m_slice = (*match).slice_y(j);

        this->iterate_bwd_edges_for_gap(left_site,&x_slice,&y_slice,&m_slice,max_x,true,j_gap_type);
        max_x->y_ind = j;

    }
    else
    {
        max_x->x_ind = -1;
        max_x->y_ind = -1;
        max_x->matrix = -1;

        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

    if(right_index>0)
    {
        right_site = right->get_site_at(right_index);

        /*
        align_slice x_slice = (*xgap)[ indices[i][ range( 0,xgap->shape()[1] ) ] ];
        align_slice y_slice = (*ygap)[ indices[i][ range( 0,ygap->shape()[1] ) ] ];
        align_slice m_slice = (*match)[ indices[i][ range( 0,match->shape()[1] ) ] ];
        */

        align_slice& x_slice = (*xgap).slice_x(i);
        align_slice& y_slice = (*ygap).slice_x(i);
        align_slice& m_slice = (*match).slice_x(i);

        this->iterate_bwd_edges_for_gap(right_site,&y_slice,&x_slice,&m_slice,max_y,false,i_gap_type);
        max_y->x_ind = i;

    }
    else
    {
        max_y->x_ind = -1;
        max_y->y_ind = -1;
        max_y->matrix = -1;

        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

    if(left_index>0 && right_index>0)
    {
        left_site = left->get_site_at(left_index);
        right_site = right->get_site_at(right_index);

        this->iterate_bwd_edges_for_match(left_site,right_site,max_m);

    }
    else
    {
        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

}

/********************************************/

void Viterbi_alignment::compute_bwd_full_score(int i,int j)
{
    int i_end = match->shape()[0];
    int j_end = match->shape()[1];

    int left_index = i;
    int right_index = j;


    Site * left_site;
    Site * right_site;

    Matrix_pointer *max_x = &(*xgap)[i][j];
    Matrix_pointer *max_y = &(*ygap)[i][j];
    Matrix_pointer *max_m = &(*match)[i][j];

    if(left_index==i_end && right_index==j_end)
    {
        return;
    }

    if(left_index<i_end)
    {
        left_site = left->get_site_at(left_index);

        //align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][j] ];
        align_slice& x_slice = (*xgap).slice_y(j);

        this->iterate_fwd_edges_for_gap(left_site,&x_slice,max_x,max_y,max_m);

    }

    if(right_index<j_end)
    {
        right_site = right->get_site_at(right_index);

        //align_slice y_slice = (*ygap)[ indices[i][ range( 0,ygap->shape()[1] ) ] ];
        align_slice& y_slice = (*ygap).slice_x(i);

        this->iterate_fwd_edges_for_gap(right_site,&y_slice,max_y,max_x,max_m);

    }

    if(left_index<i_end && right_index<j_end)
    {
        left_site = left->get_site_at(left_index);
        right_site = right->get_site_at(right_index);

        this->iterate_fwd_edges_for_match(left_site,right_site,max_x,max_y,max_m);

    }
}


void Viterbi_alignment::compute_posterior_score(int i,int j,double full_score)
{
    (*match)[i][j].full_score = (*match)[i][j].fwd_score * (*match)[i][j].bwd_score / full_score;
    (*xgap)[i][j].full_score = (*xgap)[i][j].fwd_score * (*xgap)[i][j].bwd_score / full_score;
    (*ygap)[i][j].full_score = (*ygap)[i][j].fwd_score * (*ygap)[i][j].bwd_score / full_score;
}

/********************************************/

void Viterbi_alignment::backtrack_new_path(vector<Path_pointer> *path,Path_pointer fp)
{

    //temporary use of stack instead of vector to reduce time consumption
    stack<Path_pointer> path_stack;

    vector<Edge> *left_edges = left->get_edges();
    vector<Edge> *right_edges = right->get_edges();

    int vit_mat = fp.mp.matrix;
    int x_ind = fp.mp.x_ind;
    int y_ind = fp.mp.y_ind;

    bool first_x_site = true;
    bool first_y_site = true;

    if(fp.mp.x_edge_ind>=0)
        left_edges->at(fp.mp.x_edge_ind).is_used(true);
    if(fp.mp.y_edge_ind>=0)
        right_edges->at(fp.mp.y_edge_ind).is_used(true);

    int j = match->shape()[1]-1;
    int i = match->shape()[0]-1;

    int max_j = j+1;
    int max_i = i+1;

    // Pre-existing gaps in the end skipped over
    //
    this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

    this->insert_new_path_pointer(&path_stack,&i,&j,fp);

    // Actual alignment path
    //
    while(j>=0)
    {
        while(i>=0)
        {
            if(vit_mat == Viterbi_alignment::m_mat)
            {
                if(first_x_site)
                {
                    Edge e(x_ind,max_i);
                    int in = left->get_fwd_edge_index_at_site(x_ind,&e);
                    if(in>=0)
                        left_edges->at(in).is_used(true);
                    first_x_site = false;
                }
                if(first_y_site)
                {
                    Edge e(y_ind,max_j);
                    int in = right->get_fwd_edge_index_at_site(y_ind,&e);
                    if(in>=0)
                        right_edges->at(in).is_used(true);
                    first_y_site = false;
                }

                vit_mat = (*match)[i][j].matrix;
                x_ind = (*match)[i][j].x_ind;
                y_ind = (*match)[i][j].y_ind;

                left_edges->at((*match)[i][j].x_edge_ind).is_used(true);
                right_edges->at((*match)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*match)[i][j], true );

                i--; j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else if(vit_mat == Viterbi_alignment::x_mat)
            {
                if(first_x_site)
                {
                    Edge e(x_ind,max_i);
                    int in = left->get_fwd_edge_index_at_site(x_ind,&e);
                    if(in>=0)
                        left_edges->at(in).is_used(true);
                    first_x_site = false;
                }

                vit_mat = (*xgap)[i][j].matrix;
                x_ind = (*xgap)[i][j].x_ind;
                y_ind = (*xgap)[i][j].y_ind;

                left_edges->at((*xgap)[i][j].x_edge_ind).is_used(true);

                Path_pointer pp( (*xgap)[i][j], true );

                i--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else if(vit_mat == Viterbi_alignment::y_mat)
            {
                if(first_y_site)
                {
                    Edge e(y_ind,max_j);
                    int in = right->get_fwd_edge_index_at_site(y_ind,&e);
                    if(in>=0)
                        right_edges->at(in).is_used(true);
                    first_y_site = false;
                }

                vit_mat = (*ygap)[i][j].matrix;
                x_ind = (*ygap)[i][j].x_ind;
                y_ind = (*ygap)[i][j].y_ind;

                right_edges->at((*ygap)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*ygap)[i][j], true );

                j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else
            {
                Log_output::write_out("Viterbi_alignment: incorrect backward pointer: "+Log_output::itos(vit_mat)+"\n",0);
                exit(1);
            }

            if(i<1 && j<1)
                break;

        }

        if(i<1 && j<1)
            break;

    }

    path->reserve(path_stack.size());
    while(!path_stack.empty()){
        path->push_back(path_stack.top());
        path_stack.pop();
    }

}

/********************************************/

void Viterbi_alignment::sample_new_path(vector<Path_pointer> *path,Path_pointer fp)
{
    //temporary use of stack instead of vector to reduce time consumption
    stack<Path_pointer> path_stack;

    vector<Edge> *left_edges = left->get_edges();
    vector<Edge> *right_edges = right->get_edges();

    int vit_mat = fp.mp.matrix;
    int x_ind = fp.mp.x_ind;
    int y_ind = fp.mp.y_ind;

    left_edges->at(fp.mp.x_edge_ind).is_used(true);
    right_edges->at(fp.mp.y_edge_ind).is_used(true);

    int j = match->shape()[1]-1;
    int i = match->shape()[0]-1;

    // Pre-existing gaps in the end skipped over
    //
    this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

    this->insert_new_path_pointer(&path_stack,&i,&j,fp);

    Matrix_pointer bwd_p;

    // Actual alignment path
    //
    while(j>=0)
    {
        while(i>=0)
        {

            if(vit_mat == Viterbi_alignment::m_mat)
            {
                this->iterate_bwd_edges_for_sampled_match(i,j,&bwd_p);

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

//                bwd_p.fwd_score = (*match)[i][j].fwd_score;
//                bwd_p.bwd_score = (*match)[i][j].bwd_score;
//                bwd_p.full_score = (*match)[i][j].full_score;

                left_edges->at(bwd_p.x_edge_ind).is_used(true);
                right_edges->at(bwd_p.y_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                i--; j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else if(vit_mat == Viterbi_alignment::x_mat)
            {
                this->iterate_bwd_edges_for_sampled_gap(i,j,&bwd_p,true);
                bwd_p.y_ind = j;

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

//                bwd_p.fwd_score = (*xgap)[i][j].fwd_score;
//                bwd_p.bwd_score = (*xgap)[i][j].bwd_score;
//                bwd_p.full_score = (*xgap)[i][j].full_score;

                left_edges->at(bwd_p.x_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                i--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else if(vit_mat == Viterbi_alignment::y_mat)
            {
                this->iterate_bwd_edges_for_sampled_gap(j,i,&bwd_p,false);
                bwd_p.x_ind = i;

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

//                bwd_p.fwd_score = (*ygap)[i][j].fwd_score;
//                bwd_p.bwd_score = (*ygap)[i][j].bwd_score;
//                bwd_p.full_score = (*ygap)[i][j].full_score;

                right_edges->at(bwd_p.y_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(&path_stack, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(&path_stack,&i,&j,pp);
            }
            else
            {
                Log_output::write_out("Viterbi_alignment: incorrect backward pointer: "+Log_output::itos(vit_mat)+"\n",0);
                exit(1);
            }

            if(i<1 && j<1)
                break;

        }

        if(i<1 && j<1)
            break;

    }

    path->reserve(path_stack.size());
    while(!path_stack.empty()){
        path->push_back(path_stack.top());
        path_stack.pop();
    }
}

/********************************************/

/********************************************/

void Viterbi_alignment::iterate_bwd_edges_for_gap(Site * site,align_slice *z_slice,align_slice *w_slice,
                                                 align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type)
{
    if(site->has_bwd_edge()) {

        Edge * edge = site->get_first_bwd_edge();

        this->score_gap_ext(edge,z_slice,max,is_x_matrix,gap_type);
        this->score_gap_double(edge,w_slice,max,is_x_matrix);
        this->score_gap_open(edge,m_slice,max,is_x_matrix);

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();

            this->score_gap_ext(edge,z_slice,max,is_x_matrix,gap_type);
            this->score_gap_double(edge,w_slice,max,is_x_matrix);
            this->score_gap_open(edge,m_slice,max,is_x_matrix);
        }

    }
}

/********************************************/

void Viterbi_alignment::iterate_bwd_edges_for_match(Site * left_site,Site * right_site,Matrix_pointer *max)
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

        // start corner
//        if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//        {
////            m_log_match -= model->log_non_gap();
//        }
        // error in correction for bwd computation so let's forget it now.

        double m_match = 0;
        double x_match = 0;
        double y_match = 0;

        if(compute_full_score)
        {
            double match_score = model->score(left_site->character_state,right_site->character_state);
            m_match = model->non_gap() * model->non_gap() * match_score;

            // start corner
//            if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//            {
////                m_match /= model->non_gap();
//            }

            x_match = model->gap_close() * model->non_gap() * match_score;
            y_match = model->gap_close() * model->non_gap() * match_score;

        }

        this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
        this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
        this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
            this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
            this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

        }

        // left site extra edges first
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
            this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
            this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
                this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
                this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);
            }
        }

    }
}

/********************************************/

void Viterbi_alignment::iterate_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max)
{

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        // match score & gap close scores for this match
        //
        double m_log_match = model->log_non_gap();
        double m_match = model->non_gap();

        this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
        double best_score = max->score;

        //align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
        align_slice& x_slice = (*xgap).slice_y(xgap->getSize_y()-1);
        this->score_gap_close(left_edge,&x_slice,max,true);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->y_ind = xgap->shape()[1]-1;
        }

        //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
        align_slice& y_slice = (*ygap).slice_x(xgap->getSize_x()-1);
        this->score_gap_close(right_edge,&y_slice,max,false);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->x_ind = ygap->shape()[0]-1;
        }

        // first right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
            align_slice& y_slice = (*ygap).slice_x(xgap->getSize_x()-1);
            this->score_gap_close(right_edge,&y_slice,max,false);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->x_ind = ygap->shape()[0]-1;
            }
        }

        // left site extra edges then
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            //align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
            align_slice& x_slice = (*xgap).slice_y(xgap->getSize_y()-1);
            this->score_gap_close(left_edge,&x_slice,max,true);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->y_ind = xgap->shape()[1]-1;
            }


            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                }

                //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
                align_slice& y_slice = (*ygap).slice_x(xgap->getSize_x()-1);
                this->score_gap_close(right_edge,&y_slice,max,false);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                    max->x_ind = ygap->shape()[0]-1;
                }

            }
        }

    }

    if(!compute_full_score && Settings::noise>5)
        Log_output::write_out(this->print_matrices(),5);


    if(Settings::noise>1)
    {
        stringstream ss;
        ss<<"viterbi score: "<<setprecision(8)<<max->score<<" ["<<exp(max->score)<<"] ("<<max->matrix<<")";
        if(compute_full_score)
            ss<<", full probability: "<<log(max->fwd_score)<<" ["<<max->fwd_score<<"]";
        ss<<endl;
        Log_output::write_out(ss.str(),1);
    }
}

/********************************************/

void Viterbi_alignment::iterate_fwd_edges_for_gap(Site * site,align_slice *g_slice,
                                   Matrix_pointer *max_s,Matrix_pointer *max_d,Matrix_pointer *max_m)
{
    if(site->has_fwd_edge()) {

        Edge * edge = site->get_first_fwd_edge();
        //int slice_end = (int)g_slice->shape()[0];
        int slice_end = (int)g_slice->getTotalSize();

        if(edge->get_end_site_index() < slice_end )
        {
            this->score_gap_ext_bwd(edge,g_slice,max_s);
            this->score_gap_double_bwd(edge,g_slice,max_d);
            this->score_gap_open_bwd(edge,g_slice,max_m);
        }
        while(site->has_next_fwd_edge())
        {
            edge = site->get_next_fwd_edge();

            if(edge->get_end_site_index() < slice_end )
            {
                this->score_gap_ext_bwd(edge,g_slice,max_s);
                this->score_gap_double_bwd(edge,g_slice,max_d);
                this->score_gap_open_bwd(edge,g_slice,max_m);
            }
        }
    }
}

/********************************************/

void Viterbi_alignment::iterate_fwd_edges_for_match(Site * left_site,Site * right_site,
                                     Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m)
{
    if(left_site->has_fwd_edge() && right_site->has_fwd_edge())
    {

        Edge * left_edge = left_site->get_first_fwd_edge();
        Edge * right_edge = right_site->get_first_fwd_edge();

        // PROBLEM: match score has to be computed for the cell where coming from, not where now!

        int left_end = (int)match->shape()[0];
        int right_end = (int)match->shape()[1];

        if(left_edge->get_end_site_index() < left_end
           && right_edge->get_end_site_index() < right_end )
        {
            this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
        }

        // then right site extra edges
        //
        while(right_site->has_next_fwd_edge())
        {
            right_edge = right_site->get_next_fwd_edge();
            left_edge = left_site->get_first_fwd_edge();

            if(left_edge->get_end_site_index() < left_end
               && right_edge->get_end_site_index() < right_end )
            {
                this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
            }
        }

        // left site extra edges first
        //
        while(left_site->has_next_fwd_edge())
        {
            left_edge = left_site->get_next_fwd_edge();
            right_edge = right_site->get_first_fwd_edge();

            if(left_edge->get_end_site_index() < left_end
               && right_edge->get_end_site_index() < right_end )
            {
                this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
            }

            while(right_site->has_next_fwd_edge())
            {

                right_edge = right_site->get_next_fwd_edge();

                if(left_edge->get_end_site_index() < left_end
                   && right_edge->get_end_site_index() < right_end )
                {
                    this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
                }
            }
        }
    }
}

/********************************************/

void Viterbi_alignment::iterate_bwd_edges_for_sampled_gap(int site_index1,int site_index2,Matrix_pointer *sample_p,bool is_x_matrix)
{

    align_slice *z_slice;
    align_slice *w_slice;
    align_slice *m_slice;

    Site * site;

    if(is_x_matrix)
    {
        /*
        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][site_index2] ];
        align_slice y_slice = (*ygap)[ indices[ range( 0,ygap->shape()[0] ) ][site_index2] ];
        align_slice M_slice = (*match)[ indices[ range( 0,match->shape()[0] ) ][site_index2] ];
        */

        align_slice& x_slice = (*xgap).slice_y(site_index2);
        align_slice& y_slice = (*ygap).slice_y(site_index2);
        align_slice& M_slice = (*match).slice_y(site_index2);

        z_slice = &x_slice;
        w_slice = &y_slice;
        m_slice = &M_slice;

        site = left->get_site_at(site_index1);
    }
    else
    {
        /*
        align_slice x_slice = (*xgap)[ indices[site_index2][ range( 0,xgap->shape()[1] ) ] ];
        align_slice y_slice = (*ygap)[ indices[site_index2][ range( 0,ygap->shape()[1] ) ] ];
        align_slice M_slice = (*match)[ indices[site_index2][ range( 0,match->shape()[1] ) ] ];
        */

        align_slice& x_slice = (*xgap).slice_x(site_index2);
        align_slice& y_slice = (*ygap).slice_x(site_index2);
        align_slice& M_slice = (*match).slice_x(site_index2);

        z_slice = &y_slice;
        w_slice = &x_slice;
        m_slice = &M_slice;

        site = right->get_site_at(site_index1);
    }

    if(site->has_bwd_edge()) {

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;

        Edge * edge = site->get_first_bwd_edge();        

        Matrix_pointer bwd_p;
        this->add_sample_gap_ext(edge,z_slice,&bwd_p,is_x_matrix);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_gap_double(edge,w_slice,&bwd_p,is_x_matrix);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_gap_open(edge,m_slice,&bwd_p,is_x_matrix);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();

            this->add_sample_gap_ext(edge,z_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_gap_double(edge,w_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_gap_open(edge,m_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*sample_p) = bwd_pointers.at(i);
        if(is_x_matrix)
        {
            sample_p->y_ind = xgap->shape()[1]-1;
            sample_p->fwd_score  = (*xgap)[site_index1][site_index2].fwd_score;
            sample_p->bwd_score  = (*xgap)[site_index1][site_index2].bwd_score;
            sample_p->full_score = (*xgap)[site_index1][site_index2].full_score;
        }
        else
        {
            sample_p->x_ind = ygap->shape()[0]-1;
            sample_p->fwd_score  = (*ygap)[site_index2][site_index1].fwd_score;
            sample_p->bwd_score  = (*ygap)[site_index2][site_index1].bwd_score;
            sample_p->full_score = (*ygap)[site_index2][site_index1].full_score;
        }

    }
}

/********************************************/

void Viterbi_alignment::iterate_bwd_edges_for_sampled_match(int left_index,int right_index,Matrix_pointer *sample_p)
{
    Site * left_site = left->get_site_at(left_index);
    Site * right_site = right->get_site_at(right_index);

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;


        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        double match_score = model->score(left_site->character_state,right_site->character_state);
        double m_match = model->non_gap() * model->non_gap() * match_score;

            // start corner
//            if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//            {
////                m_match /= model->non_gap();
//            }

        double x_match = model->gap_close() * model->non_gap() * match_score;
        double y_match = model->gap_close() * model->non_gap() * match_score;


        Matrix_pointer bwd_p;

        this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

        }

        // left site extra edges first
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

                this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

                this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );
            }
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*sample_p) = bwd_pointers.at(i);

        sample_p->fwd_score = (*match)[left_index][right_index].fwd_score;
        sample_p->bwd_score = (*match)[left_index][right_index].bwd_score;
        sample_p->full_score = (*match)[left_index][right_index].full_score;
    }
}

/********************************************/


void Viterbi_alignment::iterate_bwd_edges_for_sampled_end_corner(Site * left_site,Site * right_site,Matrix_pointer *end_p)
{

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;

        // match score & gap close scores for this match
        //
        double m_match = model->non_gap();

        Matrix_pointer bwd_p;
        this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        //align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
        align_slice& x_slice = (*xgap).slice_y(xgap->getSize_y()-1);
        this->add_sample_gap_close(left_edge,&x_slice,&bwd_p,true);
        bwd_p.y_ind = xgap->shape()[1]-1;

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
        align_slice& y_slice = (*ygap).slice_x(ygap->getSize_x()-1);
        this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
        bwd_p.x_ind = ygap->shape()[0]-1;

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        // first right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
            align_slice& y_slice = (*ygap).slice_x(ygap->getSize_x()-1);
            this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
            bwd_p.x_ind = ygap->shape()[0]-1;

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

        }

        // left site extra edges then
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            //align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
            align_slice& x_slice = (*xgap).slice_y(xgap->getSize_y()-1);
            this->add_sample_gap_close(left_edge,&x_slice,&bwd_p,true);
            bwd_p.y_ind = xgap->shape()[1]-1;

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );


                //align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
                align_slice& y_slice = (*ygap).slice_x(ygap->getSize_x()-1);
                this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
                bwd_p.x_ind = ygap->shape()[0]-1;

                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

            }
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*end_p) = bwd_pointers.at(i);
    }

}

/********************************************/

void Viterbi_alignment::score_m_match(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max,double m_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*match)[left_prev_index][right_prev_index].score + m_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Viterbi_alignment::m_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*match)[left_prev_index][right_prev_index].fwd_score * m_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->fwd_score += this_full_score;
    }

}

void Viterbi_alignment::score_x_match(Edge * left_edge,Edge * right_edge,double x_log_match,Matrix_pointer *max, double x_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*xgap)[left_prev_index][right_prev_index].score + x_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Viterbi_alignment::x_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*xgap)[left_prev_index][right_prev_index].fwd_score * x_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->fwd_score += this_full_score;
    }
}

void Viterbi_alignment::score_y_match(Edge * left_edge,Edge * right_edge,double y_log_match,Matrix_pointer *max, double y_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*ygap)[left_prev_index][right_prev_index].score + y_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Viterbi_alignment::y_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*ygap)[left_prev_index][right_prev_index].fwd_score * y_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->fwd_score += this_full_score;
    }
}

/********************************************/

void Viterbi_alignment::score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix,int gap_type)
{
//    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].score + model->log_gap_ext();// + edge_wght;

    // this reduces the terminal and pair-end read gap extension cost.
    //
    if(gap_type != Viterbi_alignment::normal_gap)
    {
        if(gap_type == Viterbi_alignment::end_gap)
            this_score =  (*z_slice)[prev_index].score + model->log_gap_end_ext();// + edge_wght;

        else if(gap_type == Viterbi_alignment::pair_break_gap)
            this_score =  (*z_slice)[prev_index].score + model->log_gap_break_ext();// + edge_wght;
    }

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Viterbi_alignment::x_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Viterbi_alignment::y_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*z_slice)[prev_index].fwd_score * model->gap_ext();// * this->get_edge_weight(edge);
        max->fwd_score += this_full_score;
    }
}

void Viterbi_alignment::score_gap_double(Edge *edge,align_slice *w_slice,Matrix_pointer *max,bool is_x_matrix)
{
//    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*w_slice)[prev_index].score + model->log_gap_close() + model->log_gap_open();// + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Viterbi_alignment::y_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Viterbi_alignment::x_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*w_slice)[prev_index].fwd_score * model->gap_close() * model->gap_open();// * this->get_edge_weight(edge);
        max->fwd_score += this_full_score;
    }

}

void Viterbi_alignment::score_gap_open(Edge *edge,align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix)
{
//    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*m_slice)[prev_index].score + model->log_non_gap() + this->get_log_gap_open_penalty(prev_index,is_x_matrix);// + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->matrix = Viterbi_alignment::m_mat;
        if(is_x_matrix)
        {
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*m_slice)[prev_index].fwd_score * model->non_gap() * model->gap_open();// * this->get_edge_weight(edge);
        max->fwd_score += this_full_score;
    }

}

void Viterbi_alignment::score_gap_close(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix)
{
//    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();
    int this_index = edge->get_end_site_index();

    double this_score =  (*z_slice)[prev_index].score + this->get_log_gap_close_penalty(this_index,is_x_matrix);// + edge_wght;


    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Viterbi_alignment::x_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
            max->y_edge_ind = -1;
        }
        else
        {
            max->matrix = Viterbi_alignment::y_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
            max->x_edge_ind = -1;
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*z_slice)[prev_index].fwd_score * model->gap_close();// * this->get_edge_weight(edge);
        max->fwd_score += this_full_score;
    }

}

/************************************************/

void Viterbi_alignment::score_match_bwd(Edge * left_edge,Edge * right_edge,
                                       Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m)
{
    int left_prev_index = left_edge->get_end_site_index();
    int right_prev_index = right_edge->get_end_site_index();

    double m_match = model->non_gap() * model->non_gap();// * match_score;
    double x_match = model->gap_close() * model->non_gap();// * match_score;
    double y_match = model->gap_close() * model->non_gap();// * match_score;

    double match_score =
        model->score(left->get_site_at(left_prev_index)->get_state(),right->get_site_at(right_prev_index)->get_state());

    double this_full_score =   (*match)[left_prev_index][right_prev_index].bwd_score * match_score
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    max_x->bwd_score += this_full_score  * x_match;
    max_y->bwd_score += this_full_score  * y_match;
    max_m->bwd_score += this_full_score  * m_match;
}


/********************************************/

void Viterbi_alignment::score_gap_ext_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*z_slice)[prev_index].bwd_score * model->gap_ext();// * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

void Viterbi_alignment::score_gap_double_bwd(Edge *edge,align_slice *w_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*w_slice)[prev_index].bwd_score * model->gap_close() * model->gap_open();// * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

void Viterbi_alignment::score_gap_open_bwd(Edge *edge,align_slice *m_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*m_slice)[prev_index].bwd_score * model->non_gap() * model->gap_open();// * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

/********************************************/

void Viterbi_alignment::add_sample_m_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double m_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*match)[left_prev_index][right_prev_index].fwd_score * m_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Viterbi_alignment::m_mat;
}

void Viterbi_alignment::add_sample_x_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double x_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*xgap)[left_prev_index][right_prev_index].fwd_score * x_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Viterbi_alignment::x_mat;
}

void Viterbi_alignment::add_sample_y_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double y_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*ygap)[left_prev_index][right_prev_index].fwd_score * y_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Viterbi_alignment::y_mat;
}

/************************************************/

void Viterbi_alignment::add_sample_gap_ext(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].fwd_score * model->gap_ext() * this->get_edge_weight(edge);
    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Viterbi_alignment::x_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Viterbi_alignment::y_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }
}

void Viterbi_alignment::add_sample_gap_double(Edge * edge,align_slice *w_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*w_slice)[prev_index].fwd_score * model->gap_close() * model->gap_open() * this->get_edge_weight(edge);

    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Viterbi_alignment::y_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Viterbi_alignment::x_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}

void Viterbi_alignment::add_sample_gap_open(Edge * edge,align_slice *m_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*m_slice)[prev_index].fwd_score * model->non_gap() * model->gap_open() * this->get_edge_weight(edge);

    bwd_p->score = this_score;
    bwd_p->matrix = Viterbi_alignment::m_mat;

    if(is_x_matrix)
    {
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}

void Viterbi_alignment::add_sample_gap_close(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].fwd_score * model->gap_close() * this->get_edge_weight(edge);

    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Viterbi_alignment::x_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Viterbi_alignment::y_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}


/********************************************/




/********************************************/
int Viterbi_alignment::plot_number = 1;

void Viterbi_alignment::plot_posterior_probabilities_down()
{
    string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

    string path = file;
    path.append(".mp");

    string path2 = file;
    path2.append(".tex");

    ofstream output(path.c_str(), (ios::out|ios::app));
    if (! output) { throw IOException ("Viterbi_alignment::plot_posterior_probabilities. Failed to open file"); }

    ofstream output2(path2.c_str(), (ios::out|ios::app));
    if (! output2) { throw IOException ("Viterbi_alignment::plot_posterior_probabilities. Failed to open file"); }

    output<<"beginfig("<<plot_number<<");\n"
            "u := 1mm;\ndefaultscale := 1.5pt/fontsize(defaultfont);\n"
            "path l[]; path r[]; path sqr; sqr := unitsquare scaled u;\n"
            "pickup pencircle scaled 0.1; labeloffset := 1bp; \n";

    for(unsigned int j=1;j<match->shape()[1];j++)
    {
        for(unsigned int i=1;i<match->shape()[0];i++)
        {
            if((*match)[i][j].full_score>0){
                int score = int( abs( log( (*match)[i][j].full_score ) ) );
                float red = 1;

                float green = score*7;
                if(green>255) green = 255;
                green /= 255;

                float blue = (score-39)*7;
                if(blue<0) blue = 0;
                if(blue>255) blue = 255;
                blue /= 255;

                output<<"fill sqr shifted ("<<i<<"*u,-"<<j<<"*u)\n";
                output<<"withcolor ("<<red<<","<<green<<","<<blue<<");\n";
            }
            output<<"draw sqr shifted ("<<i<<"*u,-"<<j<<"*u);\n";

        }
    }

    output<<"label.lft(\"#"<<plot_number<<"#\" infont defaultfont scaled 0.2,(0mm,1.5mm));\n";

    vector<Site> *left_sites = left->get_sites();
    vector<Site> *right_sites = right->get_sites();
    string full_alphabet = left->get_full_alphabet();

    for(unsigned int i=0;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"l"<<i<<" = circle"<<"(("<<0.5+i<<"mm,1.5mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 0;
            string place = "edgetop";
            if(start+1==stop)
                place = "edgebot";
            else if(start+2==stop)
                angle = 45;
            else if(start+3==stop)
                angle = 35;
            else if(start+4<=stop)
                angle = 25;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 0;
                place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 45;
                else if(start+3==stop)
                    angle = 35;
                else if(start+4<=stop)
                    angle = 25;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    for(unsigned int i=0;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"r"<<i<<" = circle"<<"((-0.5mm,"<<0.5-i<<"mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 270;
            string place = "edgelft";
            if(start+1==stop)
                place = "edgergt";
            else if(start+2==stop)
                angle = 225;
            else if(start+3==stop)
                angle = 235;
            else if(start+4<=stop)
                angle = 245;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 270;
                place = "edgelft";
                if(start+1==stop)
                    place = "edgergt";
                else if(start+2==stop)
                    angle = 225;
                else if(start+3==stop)
                    angle = 235;
                else if(start+4<=stop)
                    angle = 245;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    vector<Site> *sites = ancestral_sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if(tsite->get_site_type()==Site::real_site)
        {

            Site_children *offspring = sites->at(i).get_children();
            int lc = offspring->left_index;
            int rc = offspring->right_index;

            if(lc>=0 && rc>=0)
            {
                output<<"draw sqr shifted ("<<lc<<"*u,-"<<rc<<"*u) withcolor (0,0,1) withpen pencircle scaled 0.5;\n";
            }
        }
    }

    output<<"endfig;\n";
    output.close();

    output2<<"\\includegraphics[width=0.9\\columnwidth]{"<<file<<"."<<plot_number<<"}\n\n\\bigskip\n";
    output2<<"~\n\n\\bigskip\n";

    output2.close();

    Viterbi_alignment::plot_number += 1;
}

void Viterbi_alignment::plot_posterior_probabilities_up()
{
    string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

    string path = file;
    path.append(".mp");

    string path2 = file;
    path2.append(".tex");

    ofstream output(path.c_str(), (ios::out|ios::app));
    if (! output) { throw IOException ("Viterbi_alignment::plot_posterior_probabilities. Failed to open file"); }

    ofstream output2(path2.c_str(), (ios::out|ios::app));
    if (! output2) { throw IOException ("Viterbi_alignment::plot_posterior_probabilities. Failed to open file"); }

    output<<"beginfig("<<plot_number<<");\n"
            "u := 1mm;\ndefaultscale := 1.5pt/fontsize(defaultfont);\n"
            "path l[]; path r[]; path sqr; sqr := unitsquare scaled u;\n"
            "pickup pencircle scaled 0.1; labeloffset := 1bp; \n";

    for(unsigned int j=1;j<match->shape()[1];j++)
    {
        for(unsigned int i=1;i<match->shape()[0];i++)
        {
            if((*match)[i][j].full_score>0){
                int score = int( abs( log( (*match)[i][j].full_score ) ) );
                float red = 1;

                float green = score*7;
                if(green>255) green = 255;
                green /= 255;

                float blue = (score-39)*7;
                if(blue<0) blue = 0;
                if(blue>255) blue = 255;
                blue /= 255;

                output<<"fill sqr shifted ("<<i<<"*u,"<<j<<"*u)\n";
                output<<"withcolor ("<<red<<","<<green<<","<<blue<<");\n";
            }
            output<<"draw sqr shifted ("<<i<<"*u,"<<j<<"*u);\n";

        }
    }

    output<<"label.lft(\"#"<<plot_number<<"#\" infont defaultfont scaled 0.2,(0mm,-0.5mm));\n";

    vector<Site> *left_sites = left->get_sites();
    vector<Site> *right_sites = right->get_sites();
    string full_alphabet = left->get_full_alphabet();

    for(unsigned int i=0;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"l"<<i<<" = circle"<<"(("<<0.5+i<<"mm,-0.5mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 0;
            string place = "edgebot";
            if(start+1==stop)
                place = "edgetop";
            else if(start+2==stop)
                angle = 45;
            else if(start+3==stop)
                angle = 35;
            else if(start+4<=stop)
                angle = 25;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 0;
                place = "edgebot";
                if(start+1==stop)
                    place = "edgetop";
                else if(start+2==stop)
                    angle = 45;
                else if(start+3==stop)
                    angle = 35;
                else if(start+4<=stop)
                    angle = 25;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    for(unsigned int i=0;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"r"<<i<<" = circle"<<"((-0.5mm,"<<0.5+i<<"mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 90;
            string place = "edgelft";
            if(start+1==stop)
                place = "edgergt";
            else if(start+2==stop)
                angle = 135;
            else if(start+3==stop)
                angle = 125;
            else if(start+4<=stop)
                angle = 115;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 90;
                place = "edgelft";
                if(start+1==stop)
                    place = "edgergt";
                else if(start+2==stop)
                    angle = 135;
                else if(start+3==stop)
                    angle = 125;
                else if(start+4<=stop)
                    angle = 115;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    vector<Site> *sites = ancestral_sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if(tsite->get_site_type()==Site::real_site)
        {

            Site_children *offspring = sites->at(i).get_children();
            int lc = offspring->left_index;
            int rc = offspring->right_index;

            if(lc>=0 && rc>=0)
            {
                output<<"draw sqr shifted ("<<lc<<"*u,"<<rc<<"*u) withcolor (0,0,1) withpen pencircle scaled 0.5;\n";
            }
        }
    }

    output<<"endfig;\n";
    output.close();

    output2<<"\\includegraphics[width=0.9\\columnwidth]{"<<file<<"."<<plot_number<<"}\n\n\\bigskip\n";
    output2<<"~\n\n\\bigskip\n";

    output2.close();

    Viterbi_alignment::plot_number += 1;
}

/********************************************/

string Viterbi_alignment::print_matrices()
{
    stringstream ss;
    ss << fixed << noshowpos << setprecision (4);

    ss<<"m"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<(*match)[i][j].matrix<<" ";
        }
        ss<<endl;
    }
    ss<<endl;

    ss<<"m"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<setw(8)<<fixed<<(*match)[i][j].score<<" ";
        }
        ss<<endl;
    }
    ss<<endl;

    if(compute_full_score)
    {
        ss<<"m"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*match)[i][j].fwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"m"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*match)[i][j].bwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"m"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*match)[i][j].full_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;
    }

    ss<<"x"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<(*xgap)[i][j].matrix<<" ";
        }
        ss<<endl;
    }
    ss<<endl;


    ss<<"x"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<setw(8)<<fixed<<(*xgap)[i][j].score<<" ";
        }
        ss<<endl;
    }
    ss<<endl;

    if(compute_full_score)
    {
        ss<<"x"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*xgap)[i][j].fwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"x"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*xgap)[i][j].bwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"x"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*xgap)[i][j].full_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;
    }


    ss<<"y"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<(*ygap)[i][j].matrix<<" ";
        }
        ss<<endl;
    }
    ss<<endl;



    ss<<"y"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            ss<<setw(8)<<fixed<<(*ygap)[i][j].score<<" ";
        }
        ss<<endl;
    }
    ss<<endl;

    if(compute_full_score)
    {
        ss<<"y"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*ygap)[i][j].fwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"y"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*ygap)[i][j].bwd_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"y"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
                ss<<setw(8)<<log((*ygap)[i][j].full_score)<<" ";
            }
            ss<<endl;
        }
        ss<<endl;
    }
    return ss.str();
}
