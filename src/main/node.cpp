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

#include "main/node.h"
#include <iostream>

using namespace std;
using namespace ppa;


Node::~Node()
{
    if(this->has_left_child())
        delete left_child;

    if(this->has_right_child())
        delete right_child;

    if(this->node_has_sequence_object)
        delete sequence;

}

int Node::number_of_nodes = 0;
int Node::alignment_number = 0;

boost::mutex Node::list_mutex;
boost::mutex Node::model_mutex;
boost::mutex Node::log_mutex;
//boost::mutex Node::anchor_mutex;


/*******************************************************************************/

void Node::align_sequences_this_node(Model_factory *mf, bool is_reads_sequence)
{

    if(!Settings_handle::st.is("silent"))
    {
        if(!is_reads_sequence)
        {
            stringstream ss;
            ss<<" aligning node "<<this->get_name()<<" ("<<alignment_number<<"/"<<number_of_nodes<<"): "<<left_child->get_name()<<" - "<<right_child->get_name()<<".";
            Log_output::write_msg(ss.str(),0);
        }
//        else
//            Log_output::append_msg(" to node '"+left_child->get_name()+"'.",0);
        alignment_number++;
    }

    clock_t t_start=clock();

    double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();
    Evol_model model = mf->alignment_model(dist);

    stringstream ss;
    ss << "Time node::model: "<<double(clock()-t_start)/CLOCKS_PER_SEC<<"\n";
    Log_output::write_out(ss.str(),"time");

    Viterbi_alignment va;
    float tunnel_coverage = 0;

    ///---- memory usage checks ->
    long long int left_length = left_child->get_sequence()->sites_length();
    long long int right_length = right_child->get_sequence()->sites_length();

    long long int allowed_mem = 1024 * 1024 * (long long int)Settings_handle::st.get("memory-for-single-alignment").as<int>();
    long long int mem;

    if( Settings_handle::st.is("no-anchors") )
    {
        mem = va.get_predicted_memory_consumption(left_length, right_length);

        stringstream sm;
        sm << "Predicted memory usage: ~" << mem << " bytes";
        sm << " (allowed: " << allowed_mem << " bytes)" << endl;

        Log_output::write_out(sm.str(),2);

        if(mem > allowed_mem)
        {
            stringstream se;
            se << "ERROR: Memory usage over limits, allow more RAM or allow anchoring to continue." << endl;
            se << "Aligning terminated." << endl;
            Log_output::write_out(se.str(),2);
            exit(1);
        }

    }
    else
    {
        tunnel_coverage = va.define_tunnel(left_child->get_sequence(),right_child->get_sequence(), left_child->name, right_child->name, &model,true);
        mem = va.get_predicted_memory_consumption(left_length, right_length);

        stringstream sm;
        sm << "Predicted memory usage: ~" << mem << " bytes";
        sm << " (allowed: " << allowed_mem << " bytes)" << endl;
        Log_output::write_out(sm.str(),2);

        if(Settings_handle::st.is("force-gap"))
        {
            while(mem > allowed_mem)
            {
                stringstream se;
                se << "Memory usage over limits, replacing largest poorly aligned block with a gap." << endl;
                Log_output::write_out(se.str(),2);
                if(!va.replace_largest_tunnel_block_with_gap_tunnel()){
                    stringstream se;
                    se << "ERROR: Memory usage over limits (no more gaps to force), allow more RAM or decreace force-gap-threshold to continue." << endl;
                    se << "Aligning terminated." << endl;
                    Log_output::write_warning(se.str(), 0);
                    exit(1);
                }

                mem = va.get_predicted_memory_consumption(left_length, right_length);

                stringstream sm;
                sm << "Predicted memory usage: ~" << mem << " bytes";
                sm << " (allowed: " << allowed_mem << " bytes)" << endl;
                Log_output::write_out(sm.str(),2);
            }
        }
        else
        {
            if(mem > allowed_mem)
            {
                stringstream se;
                se << "ERROR: Memory usage over limits, allow more RAM or allow force-gap to continue." << endl;
                se << "Aligning terminated." << endl;
                Log_output::write_warning(se.str(), 0);
                exit(1);
            }
        }

    }
    ///---- memory usage checks end! ----

    float threshold = Settings::tunneling_coverage;

    if(tunnel_coverage <= threshold)
    {
        va.align(left_child->get_sequence(),right_child->get_sequence(),&model,
                 left_child->get_distance_to_parent(),right_child->get_distance_to_parent(), is_reads_sequence);

        ss.str(string());
        ss << "Time node::viterbi: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
        Log_output::write_out(ss.str(),"time");

        this->add_ancestral_sequence( va.get_simple_sequence() );

        if(is_reads_sequence)
        {
            this->get_sequence()->is_read_descendants(true);
            right_child->set_use_for_exonerate(true);
        }

        if(Settings::noise>2)
            this->print_alignment();

        if( Settings_handle::st.is("check-valid-graphs") )
            this->check_valid_graph();
    }
    else
    {
        ss.str(string());
        ss<<" anchoring coverage "<<tunnel_coverage<<" is above the threshold. Skipping the full alignment.";
        Log_output::write_msg(ss.str(),1);
    }

    ss.str(string());
    ss << "Time node::exit: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
    Log_output::write_out(ss.str(),"time");

}

/*******************************************************************************/

void Node::start_threaded_alignment(Model_factory *mf, int n_threads)
{

    this->number_of_nodes = this->get_number_of_leaves()-1;
    this->alignment_number = 1;

    vector<Node*> wait_nodes;
    vector<Node*> run_nodes;

    build_queues(wait_nodes, run_nodes);

    boost::thread_group threads;

    for (int i = 0; i < n_threads; ++i) {
        threads.create_thread(boost::bind(&ppa::Node::threaded_function, \
                                          this, \
                                          boost::ref(mf), \
                                          boost::ref(wait_nodes), \
                                          boost::ref(run_nodes)));
    }
        
    threads.join_all();


    if(Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors"))
        this->reconstruct_parsimony_ancestor(mf);

}

/*******************************************************************************/

void Node::start_openmp_alignment(Model_factory *mf, int n_threads)
{

    this->number_of_nodes = this->get_number_of_leaves()-1;
    this->alignment_number = 1;

    vector<Node*> wait_nodes;
    vector<Node*> run_nodes;

    build_queues(wait_nodes, run_nodes);

    omp_set_num_threads(n_threads);

    while(true) {
        #pragma omp parallel for
        for(int i = 0; i < int(run_nodes.size()); ++i) {
            run_nodes[i]->align_sequences_this_node_openmp(mf);
        }

        if (wait_nodes.empty()) {
            break;
        }

        run_nodes.clear();

        vector<Node*>::iterator it = wait_nodes.begin();
        while(it != wait_nodes.end()) {
            if((*it)->left_child->node_has_sequence_object and \
               (*it)->right_child->node_has_sequence_object) {

                run_nodes.push_back(*it);
                wait_nodes.erase(it);
                continue;
            }
            ++it;
        }

    }

    if(Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors"))
        this->reconstruct_parsimony_ancestor(mf);

}

/*******************************************************************************/

void Node::build_queues(vector<Node*>& wait_nodes, vector<Node*>& run_nodes) {
    if(not this->is_leaf()) {
        if(left_child->node_has_sequence_object and right_child->node_has_sequence_object) {
            run_nodes.push_back(this);
        }
        else {
            left_child->build_queues(wait_nodes, run_nodes);
            right_child->build_queues(wait_nodes, run_nodes);

            wait_nodes.push_back(this);
        }
    }
}

/*******************************************************************************/

void Node::threaded_function(Model_factory *mf, vector<Node*>& w_nodes, vector<Node*>& r_nodes)
{
    Node* n = NULL;
    
    try
    {
        while (true)
        {
            list_mutex.lock();
            
            if(r_nodes.empty() and w_nodes.empty()) {
                list_mutex.unlock();
                break;
            }
            
            if (not r_nodes.empty())
            {
                n = r_nodes.back();
                r_nodes.pop_back();
            }
            else {
                n = NULL;
            }
            
            list_mutex.unlock();
            
            
            if(n != NULL) {
                n->align_sequences_this_node_threaded(mf);
            }
            
            
            list_mutex.lock();
            
            if (not w_nodes.empty())
            {
                vector<Node*>::iterator li = w_nodes.begin();
                while(li != w_nodes.end())
                {
                    if((*li)->left_child->node_has_sequence_object && (*li)->right_child->node_has_sequence_object)
                    {
                        r_nodes.push_back((*li));
                        w_nodes.erase(li);
                        continue;
                    }
                    li++;
                }
            }
            
            list_mutex.unlock();
        }
    }
    catch (boost::lock_error& le)
    {
      cout << le.what() << endl;
    }
}

/*******************************************************************************/

void Node::align_sequences_this_node_threaded(Model_factory *mf)
{

    if(!Settings_handle::st.is("silent"))
    {
        stringstream ss;
        log_mutex.lock();

        ss<<" aligning node "<<this->get_name()<<" ("<<alignment_number<<"/"<<number_of_nodes<<"): "<<left_child->get_name()<<" - "<<right_child->get_name()<<".";

        Log_output::write_msg(ss.str(),0);
        alignment_number++;

        log_mutex.unlock();

    }

    double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();

    model_mutex.lock();

    Evol_model model = mf->alignment_model(dist);

    model_mutex.unlock();

    Viterbi_alignment va;
    float tunnel_coverage = 0;

    if( ! Settings_handle::st.is("no-anchors") )
    {
//        anchor_mutex.lock();

        tunnel_coverage = va.define_tunnel(left_child->get_sequence(),right_child->get_sequence(), left_child->name, right_child->name,&model,false);

//        anchor_mutex.unlock();
    }

    va.align(left_child->get_sequence(),right_child->get_sequence(),&model,
             left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

    this->add_ancestral_sequence( va.get_simple_sequence() );

}

/*******************************************************************************/

void Node::align_sequences_this_node_openmp(Model_factory *mf)
{
    if(!Settings_handle::st.is("silent"))
    {
            stringstream ss;
            #pragma omp critical(other)
            ss<<" aligning node "<<this->get_name()<<" ("<<alignment_number<<"/"<<number_of_nodes<<"): "<<left_child->get_name()<<" - "<<right_child->get_name()<<".";

            #pragma omp critical(other)
            Log_output::write_msg(ss.str(),0);

            #pragma omp critical(other)
            alignment_number++;
    }

    double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();

    Evol_model model(mf->get_sequence_data_type(), dist);

    #pragma omp critical(other)
    model = mf->alignment_model(dist);

    Viterbi_alignment va;
    float tunnel_coverage = 0;

    /*
    if( ! Settings_handle::st.is("no-anchors") )
        tunnel_coverage = va.define_tunnel(left_child->get_sequence(),right_child->get_sequence(), left_child->name, right_child->name,&model,true);*/

    ///---- memory usage checks ->
    long long int left_length = left_child->get_sequence()->sites_length();
    long long int right_length = right_child->get_sequence()->sites_length();

    long long int allowed_mem = 1024 * 1024 * (long long int)Settings_handle::st.get("memory-for-single-alignment").as<int>();
    long long int mem;

    if( Settings_handle::st.is("no-anchors") )
    {
        mem = va.get_predicted_memory_consumption(left_length, right_length);

        stringstream sm;
        sm << "Predicted memory usage: ~" << mem << " bytes";
        sm << " (allowed: " << allowed_mem << " bytes)" << endl;

        Log_output::write_out(sm.str(),2);

        if(mem > allowed_mem)
        {
            stringstream se;
            se << "ERROR: Memory usage over limits, allow more RAM or allow anchoring to continue." << endl;
            se << "Aligning terminated." << endl;
            Log_output::write_out(se.str(),2);
            exit(1);
        }

    }
    else
    {
        tunnel_coverage = va.define_tunnel(left_child->get_sequence(),right_child->get_sequence(), left_child->name, right_child->name, &model,true);
        mem = va.get_predicted_memory_consumption(left_length, right_length);

        stringstream sm;
        sm << "Predicted memory usage: ~" << mem << " bytes";
        sm << " (allowed: " << allowed_mem << " bytes)" << endl;
        Log_output::write_out(sm.str(),2);

        if(Settings_handle::st.is("force-gap"))
        {
            while(mem > allowed_mem)
            {
                stringstream se;
                se << "Memory usage over limits, replacing largest empty block in tunnel with gap." << endl;
                Log_output::write_out(se.str(),2);
                if(!va.replace_largest_tunnel_block_with_gap_tunnel()){
                    stringstream se;
                    se << "ERROR: Memory usage over limits (no more gaps to force), allow more RAM or decreace force-gap-threshold to continue." << endl;
                    se << "Aligning terminated." << endl;
                    Log_output::write_out(se.str(),2);
                    exit(1);
                }

                mem = va.get_predicted_memory_consumption(left_length, right_length);

                stringstream sm;
                sm << "Predicted memory usage: ~" << mem << " bytes";
                sm << " (allowed: " << allowed_mem << " bytes)" << endl;
                Log_output::write_out(sm.str(),2);
            }
        }
        else
        {
            if(mem > allowed_mem)
            {
                stringstream se;
                se << "ERROR: Memory usage over limits, allow more RAM or allow force-gap to continue." << endl;
                se << "Aligning terminated." << endl;
                Log_output::write_out(se.str(),2);
                exit(1);
            }
        }

    }
    ///---- memory usage checks end! ----

    va.align(left_child->get_sequence(),right_child->get_sequence(),&model,
             left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

    this->add_ancestral_sequence( va.get_simple_sequence() );
}

/*******************************************************************************/

void Node::add_sequence( Fasta_entry seq_entry, int data_type, bool gapped, bool no_trimming, bool turn_revcomp)
{
    sequence = new Sequence(seq_entry, data_type, gapped, no_trimming, turn_revcomp);
    this->node_has_sequence_object= true;
}

void Node::get_node_sequence(Fasta_entry *seq)
{
    Sequence *root = this->get_sequence();
    int root_length = root->sites_length();

    seq->name = this->get_name();

    for(int j=1;j<root_length-1;j++)
    {

        string c = Model_factory::get_ancestral_character_alphabet_at( sequence->get_site_at(j)->get_state() );

        int pstate = sequence->get_site_at(j)->get_path_state();
        int ptype  = sequence->get_site_at(j)->get_site_type();

        if( pstate == Site::xskipped || pstate == Site::yskipped || ptype == Site::non_real)
            c = "";//sequence->get_gap_symbol();

        seq->sequence.append(c);
    }
}

void Node::get_alignment(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes)
{
    vector<Node*> nodes;
    if(include_internal_nodes)
        this->get_all_nodes(&nodes);
    else
        this->get_leaf_nodes(&nodes);

    for(unsigned int i=0;i<nodes.size();i++)
    {
        Fasta_entry entry;
        entry.name = nodes.at(i)->get_name();
        entry.comment = nodes.at(i)->get_name_comment();

        aligned_sequences->push_back(entry);
    }

    get_alignment_for_nodes(aligned_sequences,include_internal_nodes);
}

void Node::get_alignment_for_nodes(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes)
{
    if(!this->sequence_site_index_needs_correcting())
    {
        Sequence *root = this->get_sequence();
        int root_length = root->sites_length();

        for(int j=1;j<root_length-1;j++)
        {
            vector<string> column;
            this->get_alignment_column_at(j,&column,include_internal_nodes);

            for(unsigned int i=0;i<aligned_sequences->size();i++)
            {
                aligned_sequences->at(i).sequence.append(column.at(i));
            }
        }

        if(Settings_handle::st.is("pileup-alignment") && Settings_handle::st.is("use-consensus"))
            this->add_root_consensus(aligned_sequences);
    }
    else
    {
        Sequence *root = this->get_sequence();
        int root_length = root->sites_length();

        for(int j=1;j<root_length;j++)
        {

            vector<Insertion_at_node> addition;
            this->additional_sites_before_alignment_column(j,&addition);

            if((int)addition.size()>0)
            {
//                cout<<"addition "<<addition.size()<<endl;
//                for(vector<Insertion_at_node>::iterator it=addition.begin();it!=addition.end();it++)
//                    cout<<" ins at "<<j<<": "<<it->node_name_wanted<<" "<<it->start_site<<" "<<it->length<<" "<<it->left_child_wanted<<endl;

                for(int l=0;l<(int)addition.size();l++)
                {
                    vector<string> column;

//                    this->get_multiple_alignment_columns_before(j,&columns,addition.at(l).node_name_wanted,
//                                                                addition.at(l).left_child_wanted,include_internal_nodes);

                    this->get_multiple_alignment_columns_before(addition.at(l),&column,include_internal_nodes);

                    for(unsigned int k=0;k<aligned_sequences->size();k++)
                    {
                        aligned_sequences->at(k).sequence.append(column.at(k));
                    }
                }
            }

            if(j<root_length-1)
            {
                vector<string> column;

                this->get_alignment_column_at(j,&column,include_internal_nodes);

                for(unsigned int i=0;i<aligned_sequences->size();i++)
                {
                    aligned_sequences->at(i).sequence.append(column.at(i));
                }
            }
        }
    }
}

void Node::get_alignment_for_reads(vector<Fasta_entry> *aligned_sequences, bool show_ref_insertions)
{
    vector<Node*> nodes;
    this->get_read_nodes_below(&nodes);
    for(unsigned int i=0;i<nodes.size();i++)
    {
        Fasta_entry entry;
        entry.name = nodes.at(i)->get_name();
        entry.comment = nodes.at(i)->get_name_comment();
        aligned_sequences->push_back(entry);
    }

    get_alignment_for_read_nodes(aligned_sequences, show_ref_insertions);
}

void Node::get_alignment_for_read_nodes(vector<Fasta_entry> *aligned_sequences, bool show_ref_insertions)
{
    Sequence *seq = this->get_sequence();
    int seq_length = seq->sites_length();

    for(int j=1;j<seq_length-1;j++)
    {
        int path_state = seq->get_site_at(j)->get_path_state();
        bool included_in_reference = this->site_in_reference(j);

        vector<string> column;
        bool has_characters = false;
        this->get_alignment_column_for_reads_at(j,&column,&has_characters);

        if(has_characters || (included_in_reference && path_state != Site::xskipped && path_state != Site::yskipped))
        {
            for(unsigned int i=0;i<aligned_sequences->size();i++)
            {
                aligned_sequences->at(i).sequence.append(column.at(i));
            }
        }
        else if(show_ref_insertions)
        {
            for(unsigned int i=0;i<aligned_sequences->size();i++)
            {
                aligned_sequences->at(i).sequence.append("-");
            }
        }

    }
}

void Node::get_alignment_column_for_reads_at(int j,vector<string> *column, bool *has_characters)
{

    if(!this->get_sequence()->is_read_sequence())
        return;

    if(this->is_leaf())
    {
        column->push_back(this->get_sequence()->get_site_at(j)->get_symbol());
        *has_characters = true;
    }
    else
    {
        Site_children *offspring = this->get_sequence()->get_site_at(j)->get_children();
        int lj = offspring->left_index;
        if(lj>=0)
        {
            left_child->get_alignment_column_for_reads_at(lj,column,has_characters);
        }
        else
        {
            int nl = left_child->get_number_of_read_leaves();

            for(int i=0;i<nl;i++)
                column->push_back(this->get_sequence()->get_gap_symbol());

        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            right_child->get_alignment_column_for_reads_at(rj,column,has_characters);
        }
        else
        {
            int nl = right_child->get_number_of_read_leaves();
            for(int i=0;i<nl;i++)
                column->push_back(this->get_sequence()->get_gap_symbol());
        }
    }
}


void Node::add_root_consensus(vector<Fasta_entry> *aligned_sequences)
{
    Sequence *root = this->get_sequence();
    int root_length = root->sites_length();
    Fasta_entry entry;
    entry.name = "consensus";
    entry.comment = "";

    int min_num_seqs = int( this->get_weighted_number_of_leaves() * (float)Settings_handle::st.get("consensus-minimum-proportion").as<float>() );

    if(min_num_seqs < Settings_handle::st.get("consensus-minimum").as<int>())
        min_num_seqs = Settings_handle::st.get("consensus-minimum").as<int>();

    for(int j=1;j<root_length-1;j++)
    {
        Site *site = root->get_site_at(j);
        int sA = site->get_sumA();
        int sC = site->get_sumC();
        int sG = site->get_sumG();
        int sT = site->get_sumT();

        if(sA+sC+sG+sT < min_num_seqs)
        {
            entry.sequence.append("-");
        }
        else{
            if(sA>sC && sA>sG && sA>sT)
                entry.sequence.append("A");
            else if(sC>sA && sC>sG && sC>sT)
                entry.sequence.append("C");
            else if(sG>sA && sG>sC && sG>sT)
                entry.sequence.append("G");
            else if(sT>sA && sT>sC && sT>sG)
                entry.sequence.append("T");
            else if(sA>sC && sA==sG && sA>sT)
                entry.sequence.append("R");
            else if(sC>sA && sC>sG && sC==sT)
                entry.sequence.append("Y");
            else if(sA==sC && sA>sG && sA>sT)
                entry.sequence.append("M");
            else if(sG>sA && sG>sC && sG==sT)
                entry.sequence.append("K");
            else if(sA>sC && sA>sG && sA==sT)
                entry.sequence.append("W");
            else if(sC>sA && sC==sG && sC>sT)
                entry.sequence.append("S");
            else if(sC>sA && sC==sG && sC==sT)
                entry.sequence.append("B");
            else if(sA>sC && sA==sG && sA==sT)
                entry.sequence.append("D");
            else if(sA==sC && sA>sG && sA==sT)
                entry.sequence.append("H");
            else if(sA==sC && sA==sG && sA>sT)
                entry.sequence.append("V");
            else if(sA==sC && sA==sG && sA==sT)
                entry.sequence.append("N");
        }
    }
    aligned_sequences->push_back(entry);
}



void Node::get_alignment_column_at(int j,vector<string> *column, bool include_internal_nodes)
{

    if(leaf)
    {
        column->push_back(sequence->get_site_at(j)->get_symbol());
    }
    else
    {
        Site_children *offspring = sequence->get_site_at(j)->get_children();
        int lj = offspring->left_index;

        if(lj>=0)
        {
            left_child->get_alignment_column_at(lj,column,include_internal_nodes);
        }
        else
        {
            int nl = left_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = left_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back(sequence->get_gap_symbol());
        }

        if(include_internal_nodes)
        {
            string c = Model_factory::get_ancestral_character_alphabet_at( sequence->get_site_at(j)->get_state() );

            int pstate = sequence->get_site_at(j)->get_path_state();
            int ptype  = sequence->get_site_at(j)->get_site_type();

            if( pstate == Site::xskipped || pstate == Site::yskipped || ptype == Site::non_real)
                c = sequence->get_gap_symbol();
            column->push_back(c);
        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            right_child->get_alignment_column_at(rj,column,include_internal_nodes);
        }
        else
        {
            int nl = right_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = right_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back(sequence->get_gap_symbol());
        }
    }
}


void Node::get_multiple_alignment_columns_before(Insertion_at_node ins,vector<string> *column,bool include_internal_nodes)
{
    if(this->is_leaf())
    {
        column->push_back(sequence->get_gap_symbol());
        return;
    }


    if(this->get_name() == ins.node_name_wanted)
    {

        if(ins.left_child_wanted)
        {

            int k = 0;
            for(int i = ins.start_site; i <= ins.end_site; i++,k++)
            {
                this->get_left_child()->get_alignment_column_at(i,column,include_internal_nodes);
            }

            if(include_internal_nodes)
            {
                column->push_back(sequence->get_gap_symbol());
            }

            this->get_right_child()->get_multiple_alignment_columns_before(ins,column,include_internal_nodes);
        }
        else
        {
            this->get_left_child()->get_multiple_alignment_columns_before(ins,column,include_internal_nodes);

            if(include_internal_nodes)
            {
                column->push_back(sequence->get_gap_symbol());
            }

            int k = 0;
            for(int i = ins.start_site; i <= ins.end_site; i++,k++)
            {
                this->get_right_child()->get_alignment_column_at(i,column,include_internal_nodes);
            }
        }
    }
    else
    {
        this->get_left_child()->get_multiple_alignment_columns_before(ins,column,include_internal_nodes);

        if(include_internal_nodes)
        {
            column->push_back(sequence->get_gap_symbol());
        }

        this->get_right_child()->get_multiple_alignment_columns_before(ins,column,include_internal_nodes);
    }
}

void Node::get_multiple_alignment_columns_before(int j,vector< vector<string> > *columns, string node_name_wanted, bool left_child_wanted,bool include_internal_nodes)
{

    if(this->is_leaf())
    {
        for(int i=0;i<(int)columns->size();i++)
        {
            columns->at(i).push_back(sequence->get_gap_symbol());
        }
        return;
    }

    int lj = -1; int rj = -1;
    if(j>=0)
    {
        Site_children *offspring = sequence->get_site_at(j)->get_children();
        lj = offspring->left_index;
        rj = offspring->right_index;
    }

    if(this->get_name() == node_name_wanted)
    {
        if(j<0)
            Log_output::write_out("Node: error: wanted node "+node_name_wanted+" but index is "+Log_output::itos(j)+"\n",1);

        if(left_child_wanted)
        {
            if(lj<0)
                Log_output::write_out("Node: error: wanted node "+node_name_wanted+" but index is "+Log_output::itos(lj)+"\n",1);

            int k = 0;
            for(int i = lj-columns->size(); i < lj; i++,k++)
            {
                this->get_left_child()->get_alignment_column_at(i,&columns->at(k),include_internal_nodes);
            }

            if(include_internal_nodes)
            {
                for(int i=0;i<(int)columns->size();i++)
                {
                    columns->at(i).push_back(sequence->get_gap_symbol());
                }
            }

            this->get_right_child()->get_multiple_alignment_columns_before(rj,columns,node_name_wanted,left_child_wanted,include_internal_nodes);
        }
        else
        {
            if(rj<0)
                Log_output::write_out("Node: error: wanted node "+node_name_wanted+" but index is "+Log_output::itos(rj)+"\n",1);

            this->get_left_child()->get_multiple_alignment_columns_before(lj,columns,node_name_wanted,left_child_wanted,include_internal_nodes);

            if(include_internal_nodes)
            {
                for(int i=0;i<(int)columns->size();i++)
                {
                    columns->at(i).push_back(sequence->get_gap_symbol());
                }
            }

            int k = 0;
            for(int i = rj-columns->size(); i < rj; i++,k++)
            {
                this->get_right_child()->get_alignment_column_at(i,&columns->at(k),include_internal_nodes);
            }
        }
    }
    else
    {
        this->get_left_child()->get_multiple_alignment_columns_before(lj,columns,node_name_wanted,left_child_wanted,include_internal_nodes);

        if(include_internal_nodes)
        {
            for(int i=0;i<(int)columns->size();i++)
            {
                columns->at(i).push_back(sequence->get_gap_symbol());
            }
        }

        this->get_right_child()->get_multiple_alignment_columns_before(rj,columns,node_name_wanted,left_child_wanted,include_internal_nodes);
    }
}

void Node::additional_sites_before_alignment_column(int j,vector<Insertion_at_node> *addition)
{
    if(this->is_leaf())
        return;


//    cout<<"additional_sites_before_alignment_column " <<this->get_name()<<" "<<j<<endl;
    Site_children *offspring = sequence->get_site_at(j)->get_children();

    int lj = offspring->left_index;
    int rj = offspring->right_index;

    if(j>0)
    {
        int prev_lj=-1;
        int prev_rj=-1;
        int sj = j;

        while(j>0)
        {
            prev_lj = sequence->get_site_at(j-1)->get_children()->left_index;

            if(prev_lj>=0)
                break;

            j--;
        }

        j = sj;
        while(j>0)
        {
            prev_rj = sequence->get_site_at(j-1)->get_children()->right_index;

            if(prev_rj>=0)
                break;

            j--;
        }

        if(lj>0 && prev_lj>=0 && lj-prev_lj != 1)
        {
            for(int k=prev_lj+1;k<lj;k++)
            {
                left_child->additional_sites_before_alignment_column(k,addition);

                Insertion_at_node ins;
                ins.node_name_wanted = this->get_name();
                ins.length = 1;
                ins.left_child_wanted = true;
                ins.start_site = k;
                ins.end_site = k;

//                cout<<" add L ("<<this->get_name()<<") "<<ins.node_name_wanted<<" "<<ins.start_site<<endl;
                addition->push_back(ins);
            }
        }

        if(rj>0 && prev_rj>=0 && rj-prev_rj != 1)
        {
            for(int k=prev_rj+1;k<rj;k++)
            {
                right_child->additional_sites_before_alignment_column(k,addition);

                Insertion_at_node ins;
                ins.node_name_wanted = this->get_name();
                ins.length = 1;
                ins.left_child_wanted = false;
                ins.start_site = k;
                ins.end_site = k;

//                cout<<" add R ("<<this->get_name()<<") "<<ins.node_name_wanted<<" "<<ins.start_site<<endl;
                addition->push_back(ins);
            }
        }
    }

    if(lj>=0)
        left_child->additional_sites_before_alignment_column(lj,addition);

    if(rj>=0)
        right_child->additional_sites_before_alignment_column(rj,addition);


}

/***************************************************************************/

void Node::write_metapost_sequence_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath c[];\ndefaultscale := 0.5;\n";
    vector<Site> *sites = this->sequence->get_sites();
    string full_alphabet = this->sequence->get_full_alphabet();

    stringstream all_chars;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = this->get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        *output<<"c"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<c<<"\","<<color<<");\n";
    }

    if(leaf)
        *output<<"label.top(btex $"<<this->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = this->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

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
                angle = 40;
            else if(start+3==stop)
                angle = 30;
            else if(start+4<=stop)
                angle = 20;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

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
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    *output<<"endfig;\n";


    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";

    ++*count;

}

void Node::write_metapost_alignment_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    vector<Site> *sites = this->sequence->get_sites();

    vector<int> left_child_index;
    vector<int> right_child_index;

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            left_child_index.push_back(i);
        }
        if(offspring->right_index>=0)
        {
            right_child_index.push_back(i);
        }
    }


    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath l[]; path r[];\ndefaultscale := 0.5;\n";
    string full_alphabet = this->sequence->get_full_alphabet();

    *output<<"l0 = circle((0cm,1.5cm),\"s\",white);\n";
    *output<<"r0 = circle((0cm,0cm),\"s\",white);\n";


    stringstream all_chars;
    for(unsigned int i=1;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *lsite = this->left_child->get_sequence()->get_site_at(offspring->left_index);

            char lc = 's';
            if(lsite->get_site_type()==Site::real_site || lsite->get_site_type()==Site::non_real)
                lc = full_alphabet.at(lsite->get_state());
            else if(lsite->get_site_type()==Site::stop_site)
                lc = 'e';

            string color = this->get_node_fill_color(lc);
            if(lsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"l"<<i<<" = circle"<<"(("<<0.5*i<<"cm,1.5cm),\""<<lc<<"\","<<color<<");\n";
        }
        if(offspring->right_index>=0)
        {
            Site *rsite = this->right_child->get_sequence()->get_site_at(offspring->right_index);

            char rc = 's';
            if(rsite->get_site_type()==Site::real_site || rsite->get_site_type()==Site::non_real)
                rc = full_alphabet.at(rsite->get_state());
            else if(rsite->get_site_type()==Site::stop_site)
                rc = 'e';

            string color = this->get_node_fill_color(rc);
            if(rsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"r"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<rc<<"\","<<color<<");\n";
        }
    }


    if(left_child->is_leaf())
        *output<<"label.top(btex $"<<left_child->get_name()<<"$ etex,(0.125cm,1.75cm));\n";
    else
    {
        string n = left_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,1.75cm));\n";
    }

    if(right_child->is_leaf())
        *output<<"label.top(btex $"<<right_child->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = right_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {

        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *tsite = left_child->get_sequence()->get_site_at(offspring->left_index);

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
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = left_child_index.at( start );
                stop  = left_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

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
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = left_child_index.at( start );
                    stop  = left_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right_child->get_sequence()->get_site_at(offspring->right_index);

            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                int start = tedge->get_start_site_index();
                int stop  = tedge->get_end_site_index();

//                int angle = 0;
//                string place = "edgebot";
//                if(start+1==stop)
//                    place = "edgetop";
//                else if(start+2==stop)
//                    angle = 320;
//                else if(start+3==stop)
//                    angle = 330;
//                else if(start+4<=stop)
//                    angle = 340;
                int angle = 0;
                string place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = right_child_index.at( start );
                stop  = right_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    start = tedge->get_start_site_index();
                    stop  = tedge->get_end_site_index();

//                    angle = 0;
//                    place = "edgebot";
//                    if(start+1==stop)
//                        place = "edgetop";
//                    else if(start+2==stop)
//                        angle = 320;
//                    else if(start+3==stop)
//                        angle = 330;
//                    else if(start+4<=stop)
//                        angle = 340;
                    angle = 0;
                    place = "edgetop";
                    if(start+1==stop)
                        place = "edgebot";
                    else if(start+2==stop)
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = right_child_index.at( start );
                    stop  = right_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }
    }

    *output<<"endfig;\n";

    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";
    *output2<<"~\n\n\\bigskip\n";

    ++*count;
}

void Node::check_valid_graph() const
{
    vector<Site> *sites = sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *ssite = &sites->at(i);
        if( ssite->has_fwd_edge() )
        {
            Edge *edge = ssite->get_first_fwd_edge();
            Site *esite = &sites->at(edge->get_end_site_index());

            if(!esite->contains_bwd_edge(edge,true))
            {
                Log_output::write_out("Node: site "+Log_output::itos(i)+" has fwd edge from "+Log_output::itos(edge->get_start_site_index())+
                                      " to "+Log_output::itos(edge->get_end_site_index())+" but no return\n",1);
            }

            while( ssite->has_next_fwd_edge() )
            {
                edge = ssite->get_next_fwd_edge();
                esite = &sites->at(edge->get_end_site_index());

                if(!esite->contains_bwd_edge(edge,true))
                {
                    Log_output::write_out("Node: site "+Log_output::itos(i)+" has fwd edge from "+Log_output::itos(edge->get_start_site_index())+
                                          " to "+Log_output::itos(edge->get_end_site_index())+" but no return\n",1);
                }
            }
        }

        if( ssite->has_bwd_edge() )
        {
            Edge *edge = ssite->get_first_bwd_edge();
            Site *esite = &sites->at(edge->get_start_site_index());

            if(!esite->contains_fwd_edge(edge,true))
            {
                Log_output::write_out("Node: site "+Log_output::itos(i)+" has bwd edge from "+Log_output::itos(edge->get_start_site_index())+
                                      " to "+Log_output::itos(edge->get_end_site_index())+" but no return\n",1);
            }

            while( ssite->has_next_bwd_edge() )
            {
                edge = ssite->get_next_bwd_edge();
                esite = &sites->at(edge->get_start_site_index());

                if(!esite->contains_fwd_edge(edge,true))
                {
                    Log_output::write_out("Node: site "+Log_output::itos(i)+" has bwd edge from "+Log_output::itos(edge->get_start_site_index())+
                                          " to "+Log_output::itos(edge->get_end_site_index())+" but no return\n",1);
                }
            }
        }

    }
}

string Node::get_sequence_string(bool with_gaps)
{

    Sequence *seq = this->get_sequence();

    int seq_length = seq->sites_length();
    string seq_str;

    bool codonstr = ( ( seq->get_data_type()==Model_factory::dna && Settings_handle::st.is("codons") ) || seq->get_data_type() == Model_factory::codon);

    if(seq->is_terminal_sequence())
    {
        for(int j=1;j<seq_length-1;j++)
            seq_str += seq->get_site_at(j)->get_symbol();
    }
    else
    {
        for(int j=1;j<seq_length-1;j++)
        {
            string c = Model_factory::get_ancestral_character_alphabet_at( seq->get_site_at(j)->get_state() );

            int pstate = seq->get_site_at(j)->get_path_state();
            int ptype  = seq->get_site_at(j)->get_site_type();

            if( pstate != Site::xskipped && pstate != Site::yskipped && ptype != Site::non_real)
                seq_str += c;
            else
                if(with_gaps)
                {
                    if(codonstr)
                        seq_str +=  "---";
                    else
                        seq_str +=  '-';
                }
        }
    }

    return seq_str;

}

void Node::prune_down()
{
    if(this->is_leaf())
        return;

    this->has_left_child(true);
    this->has_right_child(true);

    left_child->prune_down();
    right_child->prune_down();

    if(!left_child->has_sequence()){
        delete this->left_child;
        this->has_left_child(false);
    }

    if(!right_child->has_sequence()){
        delete this->right_child;
        this->has_right_child(false);
    }

    if(this->has_left_child() && !left_child->is_leaf())
    {
        if(!left_child->has_left_child() && left_child->has_right_child())
        {
            Node *new_child = left_child->right_child;
            new_child->set_distance_to_parent (left_child->get_distance_to_parent()+
                                        left_child->right_child->get_distance_to_parent());

            left_child->has_right_child(false);
            delete left_child;
            this->add_left_child(new_child);
        }
        else if(left_child->has_left_child() && !left_child->has_right_child())
        {
            Node *new_child = left_child->left_child;
            new_child->set_distance_to_parent (left_child->get_distance_to_parent()+
                                left_child->left_child->get_distance_to_parent());

            left_child->has_left_child(false);
            delete left_child;
            this->add_left_child(new_child);
        }
    }

    if(this->has_right_child() && !right_child->is_leaf())
    {
        if(!right_child->has_left_child() && right_child->has_right_child())
        {
            Node *new_child = right_child->right_child;
            new_child->set_distance_to_parent (right_child->get_distance_to_parent()+
                                        right_child->right_child->get_distance_to_parent() );

            right_child->has_right_child(false);
            delete right_child;
            this->add_right_child(new_child);
        }
        else if(right_child->has_left_child() && !right_child->has_right_child())
        {
            Node *new_child = right_child->left_child;
            new_child->set_distance_to_parent (right_child->get_distance_to_parent()+
                                            right_child->left_child->get_distance_to_parent());
            right_child->has_left_child(false);
            delete right_child;
            this->add_right_child(new_child);
        }
    }

    if(this->has_left_child() && left_child->has_sequence())
        this->has_sequence(true);

    if(this->has_right_child() && right_child->has_sequence())
        this->has_sequence(true);
}

void Node::prune_up()
{
    if(!this->is_leaf() && !this->has_left_child() && this->has_right_child())
    {
        Node* tmp_child = right_child;
        left_child = tmp_child->left_child;
        right_child = tmp_child->right_child;
        tmp_child->has_left_child(false);
        tmp_child->has_right_child(false);
        delete tmp_child;
    }

    if(!this->is_leaf() && this->has_left_child() && !this->has_right_child())
    {
        Node* tmp_child = left_child;
        left_child = tmp_child->left_child;
        right_child = tmp_child->right_child;
        tmp_child->has_left_child(false);
        tmp_child->has_right_child(false);
        delete tmp_child;
    }
}

