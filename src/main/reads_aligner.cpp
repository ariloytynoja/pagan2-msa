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

#include "main/reads_aligner.h"
#include "utils/exonerate_queries.h"
#include "utils/text_utils.h"
#include "utils/log_output.h"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace ppa;

Reads_aligner::Reads_aligner(){}

void Reads_aligner::align(Node *root, Model_factory *mf, int count)
{

    // Handle also the first one correctly
    //
    if(!Settings_handle::st.is("ref-seqfile"))
    {
        root->get_sequence()->is_read_sequence(true);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    string file = Settings_handle::st.get("queryfile").as<string>();

    Log_output::write_header("Aligning reads ",0);

    Fasta_reader fr;
    vector<Fasta_entry> reads;
    Log_output::write_out("Reads data file: "+file+"\n",1);

    try
    {
        fr.read(file, reads, true);
        fr.remove_gaps(&reads);
    }
    catch (ppa::IOException& e) {
        Log_output::write_out("Error reading the reads file '"+file+"'.\nExiting.\n\n",0);
        exit(1);
    }

    int data_type = fr.check_sequence_data_type(&reads);
    bool is_dna = data_type == Model_factory::dna;

    if(!fr.check_alphabet(&reads,data_type))
        Log_output::write_out(" Warning: Illegal characters in input reads sequences removed!\n",2);




    //////////////////////////////////////////////////////////////////
    //                                                              //
    //     Different read placement options  (until now)            //
    //                                                              //
    //////////////////////////////////////////////////////////////////


    bool has_dna_seqs = true;
    for(vector<Fasta_entry>::iterator it = reads.begin();it!=reads.end();it++)
    {
        if(it->dna_sequence.length()==0)
        {
            has_dna_seqs = false;
            break;
        }
    }


    if( Settings_handle::st.is("pileup-alignment") || Settings_handle::st.is("align-reads-at-root") )
        {
            // translated DNA with ORF search
            if( ( Settings_handle::st.is("find-orfs") ) && has_dna_seqs)
            {
                Log_output::write_header("Aligning reads: pileup with ORF search",0);
                this->translated_pileup_alignment(root,&reads,mf,count);
            }
            // default; can compare reverse strand
            else
            {
                Log_output::write_header("Aligning reads: simple pileup",0);
                this->pileup_alignment(root,&reads,mf,count);
            }
    }

#ifdef NCBI_TOOLKIT

    else if(Settings_handle::st.is("ncbi"))
    {
        mol_type mol;
        if(is_dna){
            mol = dna;
        }else{
            mol = aminoa;
        }
        this->query_placement_one_ncbi(root,&reads,mf,count,mol);
    }

#endif

    // Query placement: search for optimal node or TID tags in the tree
    else
    {
        // translated DNA with ORF search
        if( ( Settings_handle::st.is("find-orfs") ) && has_dna_seqs)
        {
            Log_output::write_header("Aligning reads: placement with ORF search",0);
            if(Settings_handle::st.is("fragments"))
                this->translated_query_placement_all(root,&reads,mf,count);
            else
                this->translated_query_placement_one(root,&reads,mf,count);
        }
        // default; can compare reverse strand
        else
        {
            Log_output::write_header("Aligning reads: simple placement",0);
            if(Settings_handle::st.is("fragments"))
                this->query_placement_all(root,&reads,mf,count,is_dna);
            else
                this->query_placement_one(root,&reads,mf,count,is_dna);
        }
    }

}

/**********************************************************************/


void Reads_aligner::pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    string ref_root_name = root->get_name();
    global_root = root;

    int start_i = 0;
    if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
        start_i = 1;

    int max_attempts = Settings_handle::st.get("query-cluster-attempts").as<int>();
    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;

    for(int j=0; j < max_attempts; j++)
    {

        for(int i=start_i;i<(int)reads->size();i++)
        {

            if(reads->at(i).cluster_attempts >= max_attempts)
                continue;

            stringstream ss;
            ss<<"#"<<count<<"#";

            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            this->create_temp_node(node,ss.str(), global_root, &reads->at(i),false);

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: "+reads->at(i).name+" "+reads->at(i).comment+".",0);

            node->align_sequences_this_node(mf,true);
            this->compute_read_overlap(node,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap,&read_identity);


            Node * node_rc = 0;
            float read_overlap_rc = -1;
            float read_identity_rc = -1;

            if(compare_reverse)
            {
                node_rc = new Node();
                this->create_temp_node(node_rc,ss.str(), global_root, &reads->at(i),true);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read (rc): "+reads->at(i).name+" "+reads->at(i).comment+".",0);

                node_rc->align_sequences_this_node(mf,true);
                this->compute_read_overlap(node_rc,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap_rc,&read_identity_rc);

                Log_output::write_out("forward overlap: "+Log_output::ftos(read_overlap)+"; backward overlap: "+Log_output::ftos(read_overlap_rc)+"\n",1);
            }
            else
            {
                Log_output::write_out("forward overlap: "+Log_output::ftos(read_overlap)+"\n",1);
            }

            reads->at(i).cluster_attempts++;


            if(read_overlap > read_overlap_rc && read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                global_root = node;

                if(compare_reverse)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }

                reads->at(i).cluster_attempts = max_attempts;

                this->fix_branch_lengths(root,node);
            }

            else if( read_overlap_rc > min_overlap && read_identity_rc > min_identity )
            {
                count++;
                global_root = node_rc;

                node->has_left_child(false);
                delete node;

                reads->at(i).cluster_attempts = max_attempts;

                this->fix_branch_lengths(root,node_rc);
            }

            else
            {
                reads->at(i).cluster_attempts++;

                node->has_left_child(false);
                delete node;

                if(compare_reverse)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }
            }
        }
    }
}

void Reads_aligner::translated_pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    this->define_translation_tables();

    string ref_root_name = root->get_name();
    global_root = root;

    int start_i = 0;
    if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
        start_i = 1;

    int max_attempts = Settings_handle::st.get("query-cluster-attempts").as<int>();
    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    for(int j=0; j < max_attempts; j++)
    {

        for(int i=start_i;i<(int)reads->size();i++)
        {

            if(reads->at(i).cluster_attempts >= max_attempts)
                continue;

            Node * best_node = 0;
            Orf *best_orf = 0;

            vector<Orf> open_frames;
            this->find_orfs(&reads->at(i),&open_frames);

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: '"+reads->at(i).name+"'",0);

            if(open_frames.size() > 0)
            {
                float best_overlap = -1;
                float best_identity = -1;

                for(int h=0;h<(int)open_frames.size();h++)
                {
                    Node * node = new Node();
                    float read_overlap = -1;
                    float read_identity = -1;

                    this->create_temp_orf_node(node,global_root, &reads->at(i), &open_frames.at(h));

                    node->align_sequences_this_node(mf,true);
                    this->compute_read_overlap(node,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap,&read_identity);

                    if(read_overlap > best_overlap || (read_overlap == best_overlap && read_overlap > read_identity) )
                    {
                        best_node = node;
                        best_orf = &open_frames.at(h);
                        best_overlap = read_overlap;
                        best_identity = read_identity;
                    }
                    else
                    {
                        node->has_left_child(false);
                        delete node;
                    }
                }

                if(best_overlap > min_overlap && best_identity > min_identity)
                {
                    reads->at(i).cluster_attempts = max_attempts;
                    stringstream cs;
                    cs<<reads->at(i).name<<"_orf1";
                    best_node->get_right_child()->set_name(cs.str());
                    cs.str("");
                    cs<<" ["<<best_orf->frame<<"."<<best_orf->start+1<<"."<<best_orf->end+1<<"]";
                    best_node->get_right_child()->add_name_comment(best_node->get_right_child()->get_name_comment()+cs.str());

                    best_node->set_nhx_tid(best_node->get_left_child()->get_nhx_tid());
                    best_node->get_right_child()->set_nhx_tid(best_node->get_left_child()->get_nhx_tid());

                    best_node->get_right_child()->set_Orf(best_orf);

                    stringstream ss;
                    ss<<"#"<<count<<"#";
                    best_node->set_name(ss.str());

                    count++;
                    global_root = best_node;

                    ss.str(string());
                    ss.precision(2);
                    ss<<"orf "<<best_orf->frame<<" ("<<best_orf->start<<"-"<<best_orf->end<<"): overlap: "<<best_overlap<<", identity: "<<best_identity<<endl;
                    Log_output::write_out(ss.str(),2);
                }
                else
                {
                    reads->at(i).cluster_attempts++;

                    best_node->has_left_child(false);
                    delete best_node;
                }
            }
        }
    }
}

void Reads_aligner::query_placement_all(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count, bool is_dna)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;

    map<string,string> target_sequences;
    this->preselect_target_sequences(root,reads,&target_sequences, is_dna);

    if((int)target_sequences.size()>1)
         this->find_targets_for_queries(root, reads, mf, &target_sequences, is_dna);
    else
        this->find_nodes_for_queries(root, reads, mf, is_dna);


    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    Log_output::write_header("Aligning query sequences",0);


    set<string> unique_nodeset;
    int reads_to_place = 0;
    int reads_discarded = 0;

    for(int i=0;i<(int)reads->size();i++)
    {
        stringstream nodestream;
        nodestream << reads->at(i).node_to_align;
        string val;
        while(nodestream >> val)
        {
            unique_nodeset.insert(val);

            if(val != "discarded_read")
                reads_to_place++;
            else
                reads_discarded++;
        }
    }

    if(unique_nodeset.find("discarded_read") != unique_nodeset.end())
        unique_nodeset.erase(unique_nodeset.find("discarded_read"));

    if(reads_to_place>0)
    {
        stringstream msg;
        msg<<"Aligning reads: "<<(int)reads->size()-reads_discarded<<"/"<<(int)reads->size()<<" with "<<reads_to_place<<" placements";
        Log_output::write_header(msg.str(),0);
    }
    else
    {
        stringstream msg;
        msg<<"None of the queries matched and nothing to align\n";
        Log_output::write_warning(msg.str(),0);
    }

    vector<string> unique_nodes;
    for(set<string>::iterator sit = unique_nodeset.begin(); sit != unique_nodeset.end(); sit++)
    {
        unique_nodes.push_back(*sit);
    }
    sort(unique_nodes.begin(),unique_nodes.end(),Reads_aligner::node_is_smaller);

    map<string,Node*> nodes_map;
    root->get_all_nodes(&nodes_map);

    map<string,int> nodes_number;

    // do one tagged node at time
    //
    for(vector<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
    {
        vector<Fasta_entry> reads_for_this;

        for(int i=0;i<(int)reads->size();i++)
        {
            stringstream nodestream;
            nodestream << reads->at(i).node_to_align;
            string val;
            while(nodestream >> val)
            {
                if(val == *sit)
                    reads_for_this.push_back(reads->at(i));
            }
        }

        this->sort_reads_vector(&reads_for_this);


        string ref_node_name = *sit;

        Node *current_root = nodes_map.find(ref_node_name)->second;
        double orig_dist = current_root->get_distance_to_parent();

        bool alignment_done = false;
        int alignments_done = 0;

        // align the remaining reads to this node
        //
        for(int i=0;i<(int)reads_for_this.size();i++)
        {
            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            if(reads_for_this.at(i).query_strand != Fasta_entry::reverse_strand)
            {
                this->create_temp_node(node,ss.str(), current_root, &reads_for_this.at(i),false);
                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read: '"+reads_for_this.at(i).name+"'",0);

                node->align_sequences_this_node(mf,true);
                this->compute_read_overlap(node,reads_for_this.at(i).name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);
            }

            Node * node_rc = 0;
            float read_overlap_rc = -1;
            float read_identity_rc = -1;

            bool reverse_computed = false;

            if(compare_reverse && reads_for_this.at(i).query_strand != Fasta_entry::forward_strand)
            {
                node_rc = new Node();
                this->create_temp_node(node_rc,ss.str(), current_root, &reads_for_this.at(i),true);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read (rc): "+reads_for_this.at(i).name+".",0);

                node_rc->align_sequences_this_node(mf,true);
                this->compute_read_overlap(node_rc,reads_for_this.at(i).name,ref_node_name,current_root->get_name(),&read_overlap_rc,&read_identity_rc);

                reverse_computed = true;
            }

            Log_output::write_out("forward overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"; backward: "+Log_output::ftos(read_overlap_rc)+"/"+Log_output::ftos(read_identity_rc)+"\n",1);


            if(read_overlap > read_overlap_rc && read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                current_root = node;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }

                // Add suffix to get unique names
                map<string,int>::iterator it = nodes_number.find(reads_for_this.at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                }
                else
                {
                    nodes_number.insert(pair<string,int>(reads_for_this.at(i).name,1));
                }
            }

            else if( read_overlap_rc > min_overlap && read_identity_rc > min_identity )
            {
                count++;
                current_root = node_rc;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                node->has_left_child(false);
                delete node;

                // Add suffix to get unique names
                map<string,int>::iterator it = nodes_number.find(reads_for_this.at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                }
                else
                {
                    nodes_number.insert(pair<string,int>(reads_for_this.at(i).name,1));
                }
            }

            else
            {
                node->has_left_child(false);
                delete node;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }
            }

            current_root->set_distance_to_parent(orig_dist);

            if(alignment_done)
            {
                if(single_ref_sequence)
                {
                    global_root = current_root;
                }
                else
                {
                    bool parent_found = this->correct_sites_index(current_root, ref_node_name, alignments_done, &nodes_map);

                    if(!parent_found)
                    {
                        global_root = current_root;
                    }
                }

                this->fix_branch_lengths(root,current_root);

                ref_node_name = current_root->get_name();
            }

        }

    }
}

void Reads_aligner::query_placement_one(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count, bool is_dna)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;


    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    string discarded_filename = "outfile";
    if(Settings_handle::st.is("outfile"))
        discarded_filename = Settings_handle::st.get("outfile").as<string>();

    discarded_filename.append(".discarded");
    fstream discarded_fstream;

    map<string,string> target_sequences;
    this->preselect_target_sequences(root,reads,&target_sequences, is_dna);

    multimap<string,string> added_sequences;

    Log_output::write_header("Aligning query sequences",0);


    for(int i=0;i<(int)reads->size();i++)
    {

        string org_nodes_to_align = reads->at(i).node_to_align;

        if((int)target_sequences.size()>0 && Settings::placement_preselection)
            this->find_targets_for_query(root, &reads->at(i), mf, &target_sequences, &added_sequences, is_dna, i==0);
        else
            this->find_nodes_for_query(root, &reads->at(i), mf, is_dna, i==0);

        vector<string> unique_nodes;

        stringstream nodestream;
        nodestream << reads->at(i).node_to_align;
        string val;

        while(nodestream >> val)
        {

            if(val != "discarded_read")
                unique_nodes.push_back(val);
        }


        if((int)unique_nodes.size()==0)
        {

            if(Settings_handle::st.is("output-discarded-queries"))
            {
                if( ! discarded_fstream.is_open() )
                    discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                discarded_fstream << ">" << reads->at(i).name << endl << reads->at(i).sequence << endl;
            }
            else
            {
                stringstream msg;
                msg<<"Query  "<< reads->at(i).name<<" has no match";
                Log_output::write_warning(msg.str(),0);
            }
            continue;
        }

        sort(unique_nodes.begin(),unique_nodes.end(),Reads_aligner::node_is_smaller);

        map<string,Node*> nodes_map;
        root->get_all_nodes(&nodes_map);

        map<string,int> nodes_number;

        // do one tagged node at time
        //
        for(vector<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
        {

            string ref_node_name = *sit;

            Node *current_root = nodes_map.find(ref_node_name)->second;
            double orig_dist = current_root->get_distance_to_parent();

            bool alignment_done = false;


            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            if(reads->at(i).query_strand != Fasta_entry::reverse_strand)
            {
                this->create_temp_node(node,ss.str(), current_root, &reads->at(i),false);
                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: '"+reads->at(i).name+"'",0);

                node->align_sequences_this_node(mf,true);
                this->compute_read_overlap(node,reads->at(i).name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);
            }

            Node * node_rc = 0;
            float read_overlap_rc = -1;
            float read_identity_rc = -1;

            bool reverse_computed = false;

            if(compare_reverse && reads->at(i).query_strand != Fasta_entry::forward_strand)
            {
                node_rc = new Node();
                this->create_temp_node(node_rc,ss.str(), current_root, &reads->at(i),true);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read (rc): "+reads->at(i).name+".",0);

                node_rc->align_sequences_this_node(mf,true);
                this->compute_read_overlap(node_rc,reads->at(i).name,ref_node_name,current_root->get_name(),&read_overlap_rc,&read_identity_rc);

                reverse_computed = true;
            }

            Log_output::write_out("forward overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"; backward: "+Log_output::ftos(read_overlap_rc)+"/"+Log_output::ftos(read_identity_rc)+"\n",1);

            string unique_read_name = reads->at(i).name;

            if(read_overlap > read_overlap_rc && read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                current_root = node;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }

                // Add suffix to get unique names
                map<string,int>::iterator it = nodes_number.find(reads->at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                    unique_read_name = n.str();
                }
                else
                {
                    nodes_number.insert(pair<string,int>(reads->at(i).name,1));
                }
            }

            else if( read_overlap_rc > min_overlap && read_identity_rc > min_identity )
            {
                count++;
                current_root = node_rc;
//                cout<<"\n\nroot "<<current_root->get_name()<<endl;
                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;

                node->has_left_child(false);
                delete node;

                // Add suffix to get unique names
                map<string,int>::iterator it = nodes_number.find(reads->at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node_rc->right_child->get_name()<<"."<<it->second;
                    node_rc->right_child->set_name(n.str());
                    it->second = it->second+1;
                }
                else
                {
                    nodes_number.insert(pair<string,int>(reads->at(i).name,1));
                }
            }

            else
            {
                node->has_left_child(false);
                delete node;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }

                if(Settings_handle::st.is("output-discarded-queries"))
                {
                    if( ! discarded_fstream.is_open() )
                        discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                    discarded_fstream << ">" << reads->at(i).name << endl << reads->at(i).sequence << endl;
                }
            }

            current_root->set_distance_to_parent(orig_dist);

            if(alignment_done)
            {
                if(single_ref_sequence)
                {
                    root = current_root;
                    single_ref_sequence = false;
                }
                else
                {
                    bool parent_found = this->correct_sites_index(current_root, ref_node_name, 1, &nodes_map);

                    if(!parent_found)
                    {
                        root = current_root;
                    }
                }


                if(Settings_handle::st.is("tid-for-subroot"))
                {
                    current_root->set_nhx_tid(current_root->get_left_child()->get_nhx_tid());
                    current_root->get_left_child()->set_nhx_tid("");
                    current_root->get_right_child()->set_nhx_tid("");
                }

                this->fix_branch_lengths(root,current_root);

                if(root->get_parent_node(current_root->get_name()) != 0)
                {
                    Node *subroot = root->get_parent_node(current_root->get_name());
                    if(subroot->get_left_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,true);
                    else if(subroot->get_right_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,false);
                }

                if(Settings::placement_preselection)
                {
                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::terminal_nodes)
                    {
                        target_sequences.insert(make_pair(unique_read_name,reads->at(i).sequence));

                        stringstream str(org_nodes_to_align);
                        string nname;

                        while(str >> nname)
                        {
                            added_sequences.insert(make_pair(nname,unique_read_name));
                        }
                    }

                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::internal_nodes)
                    {
                        target_sequences.insert(make_pair(current_root->name,current_root->get_sequence()->get_sequence_string(false)));

                        stringstream str(org_nodes_to_align);
                        string nname;

                        while(str >> nname)
                        {
                            added_sequences.insert(make_pair(nname,current_root->name));
                        }
                    }
                }
            }
            global_root = root;

        }
    }
}

#ifdef NCBI_TOOLKIT

void Reads_aligner::query_placement_one_ncbi(Node *root, vector<Fasta_entry> *iqueries, Model_factory *mf, int count, mol_type molecyle)
{
    BankBlaster blaster;

    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    float min_overlap = max( Settings_handle::st.get("min-query-overlap").as<float>(), float(0) );
    float min_identity = max( Settings_handle::st.get("min-query-identity").as<float>(), float(0) );


    string discarded_filename = "outfile";
    if(Settings_handle::st.is("outfile"))
        discarded_filename = Settings_handle::st.get("outfile").as<string>();

    discarded_filename.append(".discarded");
    fstream discarded_fstream;


    // Preselect targets: keep 'num_ncbi_targets' best hits for each query
    //
    int num_ncbi_targets = 10;
    //this->preselect_target_sequences_ncbi(root,queries,num_ncbi_targets,molecyle, mf);


    Log_output::write_header("Aligning query sequences",0);


    //

    vector<Fasta_entry> *queries = iqueries;
    vector<Fasta_entry> potential_orfs;

    if(Settings_handle::st.is("find-orfs"))
    {
        this->define_translation_tables();

        for(int i=0;i<(int)iqueries->size();i++)
        {

            vector<Orf> open_frames;
            this->find_orfs(&iqueries->at(i),&open_frames);

            for(int j=0;j<(int)open_frames.size();j++)
            {
                Fasta_entry fe;
                stringstream ss;
                ss<<iqueries->at(i).name<<"_orf"<<j+1;
                fe.name = ss.str();
                ss.str("");
                ss<<" ["<<open_frames.at(j).frame<<"."<<open_frames.at(j).start+1<<"."<<open_frames.at(j).end+1<<"]";
                fe.comment = iqueries->at(i).comment+ss.str();
                fe.sequence = open_frames.at(j).translation;
                fe.dna_sequence = open_frames.at(j).dna_sequence;
                fe.data_type = Model_factory::protein;
                fe.tid = iqueries->at(i).tid;

                potential_orfs.push_back(fe);
            }
        }

        queries = &potential_orfs;
    }
    //

    for(int i=0;i<(int)queries->size();i++)
    {

        //vector<Fasta_entry> queries_one;        //just because the
        //queries_one.push_back(queries->at(i));
        if(Settings::placement_target_nodes == Settings::tid_nodes){
            this->preselect_target_sequences_ncbi(root,&(queries->at(i)),num_ncbi_targets,molecyle, mf, &(queries->at(i).tid));
        }else{
            this->preselect_target_sequences_ncbi(root,&(queries->at(i)),num_ncbi_targets,molecyle, mf);
        }


        // Now get 'num_pagan_targets' best hits for this query
        //
        int num_pagan_targets = 3;
        vector<Alignment_set> targets;

        blaster.copyXBestAlignsFromBank(targets, queries->at(i).name, num_pagan_targets);

/*
        // place holder ->
        set<string> leafs;
        root->get_all_terminal_node_names(&leafs);
        for(set<string>::iterator it=leafs.begin();it!=leafs.end();it++)
        {
            Alignment_set hit;
            hit.seq2_id = *it;
            hit.seq1_id = queries->at(i).name;
            hit.total_score = 1;
            hit.score_strand_plus = true;

            targets.push_back(hit);

            // This is a bit problematic: all target sequences are in the same orientation;
            // if the query hits to reverse strand, it should be in that orientation for all targets.
            // We now simply pick the orientation of the last hits and use that.
            //
            queries->at(i).query_strand = Fasta_entry::forward_strand;
            if( not hit.score_strand_plus )
                queries->at(i).query_strand = Fasta_entry::reverse_strand;

            //for testing
            if(queries->at(i).name =="reverse")
                queries->at(i).query_strand = Fasta_entry::reverse_strand;
        }
        // <- place holder
*/






        if((int)targets.size()==0)
        {
            if(Settings_handle::st.is("output-discarded-queries"))
            {
                if( ! discarded_fstream.is_open() )
                    discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                discarded_fstream << ">" << queries->at(i).name << endl << queries->at(i).sequence << endl;
            }
            else
            {
                stringstream msg;
                msg<<"Query  "<< queries->at(i).name<<" has no match";
                Log_output::write_warning(msg.str(),0);
            }
            continue;
        }
        else if((int)targets.size()>1)
        {
            stringstream msg;
            msg << "("<<i+1<<"/"<<queries->size()<<") mapping query: '"<<queries->at(i).name<<" "<<queries->at(i).comment<<"'";
            Log_output::write_msg(msg.str(),0);

            //Determines orientation based on best hit's (alignment_set) orientation
            queries->at(i).query_strand = Fasta_entry::forward_strand;
            if( !targets.at(0).score_strand_plus )
            {
                queries->at(i).query_strand = Fasta_entry::reverse_strand;
                stringstream msg;
                msg << "Query '" << queries->at(i).name << "' has higher reverse strand scores! Sequence orientation will be changed.";
                Log_output::write_warning(msg.str(),2);
                queries->at(i).name += "(reverse)";
            }

            select_node_for_query(root,&targets,&queries->at(i),mf,molecyle);
        }


        map<string,Node*> nodes_map;
        root->get_all_nodes(&nodes_map);

        map<string,int> nodes_number;

        // do one tagged node at time
        //
        for(vector<Alignment_set>::iterator sit = targets.begin(); sit != targets.end(); sit++)
        {

            string ref_node_name = (*sit).seq2_id;
            string unique_query_name = queries->at(i).name;

            Node *current_root = nodes_map.find(ref_node_name)->second;
            double orig_dist = current_root->get_distance_to_parent();

            bool alignment_done = false;


            Node * node = new Node();
            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());


            bool is_reverse = false;
            if( ( Settings_handle::st.is("both-strands") )
                                && mf->get_sequence_data_type()==Model_factory::dna )
            {
                is_reverse = queries->at(i).query_strand==Fasta_entry::reverse_strand;
            }


            this->create_temp_node(node,ss.str(), current_root, &queries->at(i),is_reverse);
            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(queries->size())+") aligning read: '"+queries->at(i).name+"'",0);
            node->align_sequences_this_node(mf,true);

            float read_overlap = -1;
            float read_identity = -1;
            this->compute_read_overlap(node,queries->at(i).name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);

            Log_output::write_out("overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"\n",1);


            if(read_overlap > min_overlap && read_identity > min_identity)
            {
                alignment_done = true;

                count++;
                current_root = node;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();


                // Add suffix to get unique names
                map<string,int>::iterator it = nodes_number.find(queries->at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                    unique_query_name = n.str();
                }
                else
                {
                    nodes_number.insert(pair<string,int>(queries->at(i).name,1));
                }
            }

            else
            {
                node->has_left_child(false);
                delete node;

                if(Settings_handle::st.is("output-discarded-queries"))
                {
                    if( ! discarded_fstream.is_open() )
                        discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                    discarded_fstream << ">" << queries->at(i).name << endl << queries->at(i).sequence << endl;
                }
            }



            current_root->set_distance_to_parent(orig_dist);

            if(alignment_done)
            {
                if(single_ref_sequence)
                {
                    root = current_root;
                    single_ref_sequence = false;
                }
                else
                {
                    bool parent_found = this->correct_sites_index(current_root, ref_node_name, 1, &nodes_map);

                    if(!parent_found)
                    {
                        root = current_root;
                    }
                }

                this->fix_branch_lengths(root,current_root);

                if(root->get_parent_node(current_root->get_name()) != 0)
                {
                    Node *subroot = root->get_parent_node(current_root->get_name());
                    if(subroot->get_left_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,true);
                    else if(subroot->get_right_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,false);
                }

                // Add newly aligned sequence and/or its parent to the targets
                // Function for this complete missing

                if(Settings_handle::st.is("tid-for-subroot"))
                {
                    current_root->set_nhx_tid(current_root->get_left_child()->get_nhx_tid());
                    current_root->get_left_child()->set_nhx_tid("");
                    current_root->get_right_child()->set_nhx_tid("");

                    blaster.insertOneTarget(current_root->get_name(), current_root->get_sequence()->get_sequence_string(false));
                }
                else
                {
                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::terminal_nodes)
                    {
                        blaster.insertOneTarget(unique_query_name, queries->at(i).sequence);
                    }

                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::internal_nodes)
                    {
                        blaster.insertOneTarget(current_root->name,current_root->get_sequence()->get_sequence_string(false));

                    }
                }

            }
            global_root = root;

        }
    }
}

#endif

void Reads_aligner::fix_branch_lengths(Node *root,Node *current_root)
{
    if(root->get_parent_node(current_root->get_name()) != 0)
    {
        Node *subroot = root->get_parent_node(current_root->get_name());

        vector<Fasta_entry> subalignment;
        subroot->get_alignment(&subalignment,true);
        vector<Fasta_entry>::iterator it = subalignment.begin();

        Fasta_entry lnode,rnode,pnode,tnode;

        for(;it!=subalignment.end();it++)
        {
            if(it->name == current_root->left_child->get_name())
                lnode = *it;
            if(it->name == current_root->right_child->get_name())
                rnode = *it;
            if(it->name == current_root->get_name())
                tnode = *it;
            if(it->name == subroot->get_name())
                pnode = *it;
        }

        int share12=0; int share13=0; int share23=0;
        int ident12=0; int ident13=0; int ident23=0;

        for(int i=0;i<(int)pnode.sequence.length();i++)
        {            
            if(pnode.sequence.at(i)!='-' && pnode.sequence.at(i)!='.')
            {
                if(lnode.sequence.at(i)!='-' && lnode.sequence.at(i)!='.')
                {
                    share12++;
                    if(pnode.sequence.at(i)==lnode.sequence.at(i))
                        ident12++;
                }

                if(rnode.sequence.at(i)!='-' && rnode.sequence.at(i)!='.')
                {
                    share13++;
                    if(pnode.sequence.at(i)==rnode.sequence.at(i))
                        ident13++;
                }
            }
            if(lnode.sequence.at(i)!='-' && lnode.sequence.at(i)!='.' &&
                    rnode.sequence.at(i)!='-' && rnode.sequence.at(i)!='.')
            {
                share23++;
                if(lnode.sequence.at(i)==rnode.sequence.at(i))
                    ident23++;
            }
        }

        float d12 = 1-(float)ident12/(float)share12;
        float d13 = 1-(float)ident13/(float)share13;
        float d23 = 1-(float)ident23/(float)share23;

        float l2 = 0.5*d23 + 0.5*(d12-d13);
        float l3 = 0.5*d23 + 0.5*(d13-d12);
        float l1 = 0.5*(d12+d13-d23);

        float mult = 1;
//        if(Settings_handle::st.is("ref-guidetree"))
        if((l1+l2)>0)
            mult = (current_root->dist_to_parent+current_root->left_child->dist_to_parent)/(l1+l2);

        l1*=mult;
        l2*=mult;
        l3*=mult;

        current_root->set_distance_to_parent(l1);
        current_root->left_child->set_distance_to_parent(l2);
        current_root->right_child->set_distance_to_parent(l3);
    }
    else
    {
        vector<Fasta_entry> subalignment;
        current_root->get_alignment(&subalignment,true);
        vector<Fasta_entry>::iterator it = subalignment.begin();

        Fasta_entry lnode,rnode;

        for(;it!=subalignment.end();it++)
        {
            if(it->name == current_root->left_child->get_name())
                lnode = *it;
            if(it->name == current_root->right_child->get_name())
                rnode = *it;
        }

        int share=0;
        int ident=0;

        for(int i=0;i<(int)lnode.sequence.length();i++)
        {
            if(lnode.sequence.at(i)!='-' && lnode.sequence.at(i)!='.' &&
                    rnode.sequence.at(i)!='-' && rnode.sequence.at(i)!='.')
            {
                share++;
                if(lnode.sequence.at(i)==rnode.sequence.at(i))
                    ident++;
            }
        }

        float d = (1-(float)ident/(float)share)/2;

        current_root->left_child->set_distance_to_parent(d);
        current_root->right_child->set_distance_to_parent(d);
    }
}

void Reads_aligner::translated_query_placement_all(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    this->define_translation_tables();

    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    vector<Fasta_entry> orfs;
    vector<Fasta_entry> potential_orfs;

    for(int i=0;i<(int)reads->size();i++)
    {
        vector<Orf> open_frames;
        this->find_orfs(&reads->at(i),&open_frames);

        for(int j=0;j<(int)open_frames.size();j++)
        {
            Fasta_entry fe;
            stringstream ss;
            ss<<reads->at(i).name<<"_orf"<<j+1;
            fe.name = ss.str();
            ss.str("");
            ss<<" ["<<open_frames.at(j).frame<<"."<<open_frames.at(j).start+1<<"."<<open_frames.at(j).end+1<<"]";
            fe.comment = reads->at(i).comment+ss.str();
            fe.sequence = open_frames.at(j).translation;
            fe.dna_sequence = open_frames.at(j).dna_sequence;
            fe.data_type = Model_factory::protein;
            fe.tid = reads->at(i).tid;

            potential_orfs.push_back(fe);
        }
    }

    map<string,string> target_sequences;
    this->preselect_target_sequences(root,&potential_orfs,&target_sequences, false);

    if((int)target_sequences.size()>1)
         this->find_targets_for_queries(root, &potential_orfs, mf, &target_sequences, false);
    else
        this->find_nodes_for_queries(root, &potential_orfs, mf, false);

    for(int j=0;j<(int)potential_orfs.size();j++)
    {
        if(potential_orfs.at(j).node_to_align != "discarded_read")
            orfs.push_back(potential_orfs.at(j));
    }


    Log_output::write_header("Aligning query sequences",0);

    set<string> unique_nodeset;
    int orfs_to_place = 0;

    for(int i=0;i<(int)orfs.size();i++)
    {
        stringstream nodestream;
        nodestream << orfs.at(i).node_to_align;
        string val;
        while(nodestream >> val)
        {
            unique_nodeset.insert(val);

            orfs_to_place++;
        }
    }

    if(unique_nodeset.find("discarded_read") != unique_nodeset.end())
        unique_nodeset.erase(unique_nodeset.find("discarded_read"));

    if(orfs_to_place>0)
    {
        stringstream msg;
        msg<<"\nAligning ORFs: "<<(int)orfs.size()<<"/"<<(int)reads->size()<<" with "<<orfs_to_place<<" placements";
//        Log_output::clean_output();
        Log_output::write_new_header(msg.str(),0);
    }
    else
    {
        stringstream msg;
        if(potential_orfs.size()==0)
            msg<<"\nNo long enough ORFs and nothing to align\n";
        else
            msg<<"\nNone of the ORFs matched and nothing to align\n";

        Log_output::write_warning(msg.str(),0);
    }


    vector<string> unique_nodes;
    for(set<string>::iterator sit = unique_nodeset.begin(); sit != unique_nodeset.end(); sit++)
    {
        unique_nodes.push_back(*sit);
    }
    sort(unique_nodes.begin(),unique_nodes.end(),Reads_aligner::node_is_smaller);

    map<string,Node*> nodes_map;
    root->get_all_nodes(&nodes_map);

    map<string,int> nodes_number;

    // do one tagged node at time
    //
    for(vector<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
    {
        vector<Fasta_entry> reads_for_this;

        for(int i=0;i<(int)orfs.size();i++)
        {
            stringstream nodestream;
            nodestream << orfs.at(i).node_to_align;
            string val;
            while(nodestream >> val)
            {
                if(val == *sit)
                    reads_for_this.push_back(orfs.at(i));
            }
        }

        this->sort_reads_vector(&reads_for_this);


        string ref_node_name = *sit;
//        cout<<"REF "<<ref_node_name<<endl;

        Node *current_root = nodes_map.find(ref_node_name)->second;
        double orig_dist = current_root->get_distance_to_parent();

        bool alignment_done = false;
        int alignments_done = 0;

        // align the remaining reads to this node
        //
        for(int i=0;i<(int)reads_for_this.size();i++)
        {
            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            this->create_temp_node(node,ss.str(), current_root, &reads_for_this.at(i),false);

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read: '"+reads_for_this.at(i).name+"'",0);

            node->align_sequences_this_node(mf,true);
            this->compute_read_overlap(node,reads_for_this.at(i).name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);

            Log_output::write_out("forward overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"\n",1);

            if(read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                current_root = node;

//                cout<<"\nis read "<<current_root->get_right_child()->get_sequence()->is_read_sequence()<<"\n\n";
                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                map<string,int>::iterator it = nodes_number.find(reads_for_this.at(i).name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                }
                else
                {
                    nodes_number.insert(pair<string,int>(reads_for_this.at(i).name,1));
                }

            }

            else
            {
                node->has_left_child(false);
                delete node;
            }


            current_root->set_distance_to_parent(orig_dist);

            if(alignment_done)
            {
                if(single_ref_sequence)
                {
                    global_root = current_root;
                }
                else
                {
                    bool parent_found = this->correct_sites_index(current_root, ref_node_name, alignments_done, &nodes_map);

                    if(!parent_found)
                    {
                        global_root = current_root;
                    }
                }

                this->fix_branch_lengths(root,current_root);

                ref_node_name = current_root->get_name();
            }

            global_root = root;

        }
    }
}


void Reads_aligner::translated_query_placement_one(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    this->define_translation_tables();

    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    string discarded_filename = "outfile";
    if(Settings_handle::st.is("outfile"))
        discarded_filename = Settings_handle::st.get("outfile").as<string>();

    discarded_filename.append(".discarded");
    fstream discarded_fstream;

    vector<Fasta_entry> potential_orfs;

    for(int i=0;i<(int)reads->size();i++)
    {

        vector<Orf> open_frames;
        this->find_orfs(&reads->at(i),&open_frames);

        for(int j=0;j<(int)open_frames.size();j++)
        {
            Fasta_entry fe;
            stringstream ss;
            ss<<reads->at(i).name<<"_orf"<<j+1;
            fe.name = ss.str();
            ss.str("");
            ss<<" ["<<open_frames.at(j).frame<<"."<<open_frames.at(j).start+1<<"."<<open_frames.at(j).end+1<<"]";
            fe.comment = reads->at(i).comment+ss.str();
            fe.sequence = open_frames.at(j).translation;
            fe.dna_sequence = open_frames.at(j).dna_sequence;
            fe.data_type = Model_factory::protein;
            fe.tid = reads->at(i).tid;

            potential_orfs.push_back(fe);
        }
    }

    map<string,string> target_sequences;
    this->preselect_target_sequences(root,&potential_orfs,&target_sequences, false);

    multimap<string,string> added_sequences;

    Log_output::write_header("Aligning query sequences",0);

    map<string,int> nodes_number;

    for(int i=0;i<(int)potential_orfs.size();i++)
    {
        string org_nodes_to_align = potential_orfs.at(i).node_to_align;

        Fasta_entry potential_orf = potential_orfs.at(i);

        if((int)target_sequences.size()>0)
            this->find_targets_for_query(root, &potential_orf, mf, &target_sequences, &added_sequences, false, i==0);
        else
            this->find_nodes_for_query(root, &potential_orf, mf, false);

        if(potential_orf.node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("output-discarded-queries"))
            {
                if( ! discarded_fstream.is_open() )
                    discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                discarded_fstream << ">" << potential_orf.name << endl << potential_orf.sequence << endl;
            }

            continue;
        }


        map<string,Node*> nodes_map;
        root->get_all_nodes(&nodes_map);


        stringstream nodestream;
        nodestream << potential_orf.node_to_align;

        string ref_node_name;

        while(nodestream >> ref_node_name)
        {
//            cout<<"REF "<<ref_node_name<<endl;

            Node *current_root = nodes_map.find(ref_node_name)->second;
            double orig_dist = current_root->get_distance_to_parent();

            bool alignment_done = false;
            int alignments_done = 0;


            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            this->create_temp_node(node,ss.str(), current_root, &potential_orf,false);

            Log_output::write_msg(" aligning read: '"+potential_orf.name+"'",0);

            node->align_sequences_this_node(mf,true);
            this->compute_read_overlap(node,potential_orf.name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);

            Log_output::write_out("forward overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"\n",1);

            string unique_potential_orf_name = potential_orf.name;

            if(read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                current_root = node;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                map<string,int>::iterator it = nodes_number.find(potential_orf.name);
                if(it!=nodes_number.end())
                {
                    stringstream n;
                    n<<node->right_child->get_name()<<"."<<it->second;
                    node->right_child->set_name(n.str());
                    it->second = it->second+1;
                    unique_potential_orf_name = n.str();
                }
                else
                {
                    nodes_number.insert(pair<string,int>(potential_orf.name,1));
                }

            }

            else
            {
                node->has_left_child(false);
                delete node;

                if(Settings_handle::st.is("output-discarded-queries"))
                {
                    if( ! discarded_fstream.is_open() )
                        discarded_fstream.open(discarded_filename.c_str(), fstream::out);
                    discarded_fstream << ">" << potential_orf.name << endl << potential_orf.sequence << endl;
                }
            }



            current_root->set_distance_to_parent(orig_dist);

            if(alignment_done)
            {
                if(single_ref_sequence)
                {
                    root = current_root;
                    single_ref_sequence = false;
                }
                else
                {
                    bool parent_found = this->correct_sites_index(current_root, ref_node_name, alignments_done, &nodes_map);

                    if(!parent_found)
                    {
                        root = current_root;
                    }
                }

                /* untested here */
                if(Settings_handle::st.is("tid-for-subroot"))
                {
                    current_root->set_nhx_tid(current_root->get_left_child()->get_nhx_tid());
                    current_root->get_left_child()->set_nhx_tid("");
                    current_root->get_right_child()->set_nhx_tid("");
                }
                /* untested here */

                this->fix_branch_lengths(root,current_root);

                /* untested here */
                if(root->get_parent_node(current_root->get_name()) != 0)
                {
                    Node *subroot = root->get_parent_node(current_root->get_name());
                    if(subroot->get_left_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,true);
                    else if(subroot->get_right_child()->get_name() == current_root->get_name())
                        subroot->reconstruct_one_parsimony_ancestor(mf,false);
                }
                /* untested here */

                if(Settings::placement_preselection)
                {
                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::terminal_nodes)
                    {
                        target_sequences.insert(make_pair(unique_potential_orf_name,potential_orf.sequence));

                        stringstream str(org_nodes_to_align);
                        string nname;

                        while(str >> nname)
                        {
                            added_sequences.insert(make_pair(nname,unique_potential_orf_name));
                        }
                    }

                    if(Settings::placement_target_nodes == Settings::all_nodes ||
                           Settings::placement_target_nodes == Settings::internal_nodes)
                    {
                        target_sequences.insert(make_pair(current_root->name,current_root->get_sequence()->get_sequence_string(false)));

                        stringstream str(org_nodes_to_align);
                        string nname;

                        while(str >> nname)
                        {
                            added_sequences.insert(make_pair(nname,current_root->name));
                        }
                    }
                }
            }
            global_root = root;
        }
    }
}


/**********************************************************************/

void Reads_aligner::find_targets_for_query(Node *root, Fasta_entry *read, Model_factory *mf,map<string,string> *target_sequences, multimap<string,string> *added_sequences, bool is_dna, bool warnings)
{

    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        if(warnings)
            Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }


    // Get the original target nodes
    //
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;

    // Handle one read at time
    //
    read->node_score = -1.0;
    read->query_strand = Fasta_entry::unknown_strand;

    string tid = read->tid;
    if( tid=="" )
        tid = "<empty>";

    // Targets for this includes the original hit plus all new sequences that are aligned to those hits
    //
    map<string,string> targets_for_this;

//    map<string,string>::iterator it = target_sequences->begin();
//    for(;it!=target_sequences->end();it++)
//        cout<<"target: "<<it->first<<endl;

    stringstream str(read->node_to_align);
    string nodename;
    while(str >> nodename)
    {
        map<string,string>::iterator it = target_sequences->find(nodename);
        if(it!=target_sequences->end())
        {
            targets_for_this.insert(make_pair(it->first,it->second));
            pair <multimap<string,string>::iterator, multimap<string,string>::iterator> ret = added_sequences->equal_range(nodename);
            for (multimap<string,string>::iterator mit=ret.first; mit!=ret.second; mit++)
            {
//                cout<<"find "<<mit->second<<endl;
                it = target_sequences->find(mit->second);
//                cout<<"crash "<<it->first<<" "<<it->second<<endl;
                targets_for_this.insert(make_pair(it->first,it->second));
            }
        }
    }


    Log_output::write_msg("mapping query: '"+read->name+" "+read->comment+"'",0);

    // Call Exonerate to reduce the search space
    //
    map<string,hit> exonerate_hits;

    if(has_exonerate && Settings::exonerate_gapped_keep_best > 0 )
    {
        er.local_alignment(&targets_for_this,read,&exonerate_hits,false,is_dna);
    }


    // Handle (or reject) reads discarded by Exonerate
    //
    if(read->node_to_align == "discarded_read" || exonerate_hits.size()==0)
    {
        stringstream ss;
        ss<<"Read "<<read->name<<" with the tid "<<tid<<" was discarded by Exonerate.\n";
        Log_output::write_out(ss.str(),2);

        return;
    }


    // Has TID or exhaustive search: now find the one with the best match
    //
    int matches = exonerate_hits.size();


    // Has only one matching node
    //
    if(matches == 1)
    {
        string target_node = exonerate_hits.begin()->first;

        stringstream ss;
        ss<<"Read "<<read->name<<" with the tid "<<tid<<" only matches the node "<<target_node<<"."<<endl;
        Log_output::write_out(ss.str(),2);

        read->node_to_align = target_node;
    }


    // Has several matching nodes
    //
    else
    {
        double best_score = -HUGE_VAL;
        string best_node = root->get_name();
        int query_strand = Fasta_entry::forward_strand;

        map<string,Node*> nodes;
        root->get_all_nodes(&nodes);

        int matching_nodes = exonerate_hits.size();

        stringstream ss;
        ss<<"Read "<<read->name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
        Log_output::write_out(ss.str(),2);

        map<string,hit>::iterator it = exonerate_hits.begin();

        for(;it != exonerate_hits.end(); it++)
        {
            string target_node = it->first;

            map<string,Node*>::iterator nit = nodes.find(target_node);
            double score = this->read_match_score( nit->second, read, mf);

            stringstream ss;
            ss<<target_node<<" with score "<<score<<"\n";
            Log_output::write_out(ss.str(),2);

            if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
            {
                best_score = score;
                best_node.append(" "+target_node);
            }
            else if(score>=best_score)
            {
                best_score = score;
                best_node = target_node;
            }

            if(compare_reverse)
            {
                Fasta_entry rev_seq = *read;
                rev_seq.sequence = this->reverse_complement(rev_seq.sequence);

                double score = this->read_match_score( nit->second, &rev_seq, mf);

                stringstream ss;
                ss<<target_node<<"(rc) with score "<<score<<"\n";
                Log_output::write_out(ss.str(),2);

                if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                {
                    best_score = score;
                    best_node.append(" "+target_node);
                    query_strand = Fasta_entry::reverse_strand;
                }
                else if(score>=best_score)
                {
                    best_score = score;
                    best_node = target_node;
                    query_strand = Fasta_entry::reverse_strand;
                }
            }

            matching_nodes--;
        }

        if(best_score<0.05)
        {
            Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);
            read->node_to_align = "discarded_read";
        }
        else
        {
            stringstream ss;
            ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
            Log_output::write_out(ss.str(),2);

            read->node_score = best_score;
            read->node_to_align = best_node;
            read->query_strand = query_strand;
        }
    }
    // All done; continue
}



void Reads_aligner::find_nodes_for_query(Node *root, Fasta_entry *read, Model_factory *mf,bool is_dna,bool warnings)
{
    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }


    // Get the original target nodes
    //
    multimap<string,string> tid_nodes;
    bool ignore_tid_tags = true;
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;

    bool preselected_nodes = false;
    if(read->node_to_align != "")
    {
        stringstream str(read->node_to_align);
        string nodename;
        while(str >> nodename)
        {
            tid_nodes.insert(make_pair(nodename,nodename));
        }

        preselected_nodes = true;
    }
    else
    {
        this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
    }

    if( !Settings_handle::st.is("all-nodes") &&
        !Settings_handle::st.is("internal-nodes") &&
        !Settings_handle::st.is("terminal-nodes") &&
        !Settings_handle::st.is("test-every-internal-node") &&
        !Settings_handle::st.is("test-every-terminal-node") &&
        !Settings_handle::st.is("test-every-node") && ignore_tid_tags)
    {
        if(warnings)
            Log_output::write_warning("No tagged nodes found. Considering all nodes!",0);
    }

    // Handle one read at time
    //
    read->node_score = -1.0;
    read->query_strand = Fasta_entry::unknown_strand;

    string tid = read->tid;
    if( ignore_tid_tags )
        tid = "<empty>";

    Log_output::write_msg("mapping query: '"+read->name+" "+read->comment+"'",0);

    // Call Exonerate to reduce the search space
    //
    map<string,hit> exonerate_hits;

    if(has_exonerate && Settings::exonerate_local_keep_best > 0 && !preselected_nodes)
    {
        tid_nodes.clear();
        this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

        if(tid_nodes.size()>1)
            er.local_alignment(root,read,&tid_nodes,&exonerate_hits, true,is_dna,ignore_tid_tags);

        if(tid_nodes.size()>1 && Settings::exonerate_gapped_keep_best > 0 )
            er.local_alignment(root,read,&tid_nodes,&exonerate_hits,false,is_dna,ignore_tid_tags);
    }
    else if(has_exonerate && Settings::exonerate_gapped_keep_best > 0 )
    {
        if(!preselected_nodes)
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
        }
        if(tid_nodes.size()>1)
            er.local_alignment(root,read,&tid_nodes,&exonerate_hits,false,is_dna,ignore_tid_tags);
    }


    // Handle (or reject) reads discarded by Exonerate
    //
    if(read->node_to_align == "discarded_read")
    {
        if(Settings_handle::st.is("exhaustive-placement"))
        {
            read->node_to_align = "";
            tid_nodes.clear();

            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
        }
        else
        {
            stringstream ss;
            ss<<"Read "<<read->name<<" with the tid "<<tid<<" was discarded by Exonerate.\n";
            Log_output::write_out(ss.str(),2);

            return;
        }
    }


    // Has TID or exhaustive search: now find the one with the best match
    //
    if(tid != "")
    {
        int matches = tid_nodes.count(tid);

        if( ignore_tid_tags )
            matches = tid_nodes.size();


        // Has TID but no matching node
        //
        if(matches == 0)
        {
            stringstream ss;
            ss<<"Read "<<read->name<<" with the tid "<<tid<<" has no matching node. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            read->node_to_align = root->get_name();
        }
        // All done; continue


        // Has only one matching node and no need for ranking
        //
        else if(matches == 1)
        {
            multimap<string,string>::iterator tit = tid_nodes.find(tid);

            if( ignore_tid_tags )
                tit = tid_nodes.begin();

            if(tit != tid_nodes.end())
            {
                stringstream ss;
                ss<<"Read "<<read->name<<" with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                Log_output::write_out(ss.str(),2);

                read->node_to_align = tit->second;
            }
        }
        // All done; continue


        // Has TID and matching nodes, or exhaustive search
        //
        else
        {
            double best_score = -HUGE_VAL;
            string best_node = root->get_name();
            int query_strand = Fasta_entry::forward_strand;

            map<string,Node*> nodes;
            root->get_all_nodes(&nodes);

            multimap<string,string>::iterator tit;
            int matching_nodes = 0;

            if( ignore_tid_tags )
            {
                tit = tid_nodes.begin();
                matching_nodes = tid_nodes.size();
            }
            else
            {
                tit = tid_nodes.find(tid);
                matching_nodes = tid_nodes.count(tid);
            }

            if(tit != tid_nodes.end())
            {

                stringstream ss;
                ss<<"Read "<<read->name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                Log_output::write_out(ss.str(),2);

                while(tit != tid_nodes.end() && matching_nodes>0)
                {
                    map<string,Node*>::iterator nit = nodes.find(tit->second);
                    double score = this->read_match_score( nit->second, read, mf);

                    stringstream ss;
                    ss<<tit->second<<" with score "<<score<<"\n";
                    Log_output::write_out(ss.str(),2);

                    if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                    {
                        best_score = score;
                        best_node.append(" "+tit->second);
                        query_strand = Fasta_entry::forward_strand;
                    }
                    else if(score>=best_score)
                    {
                        best_score = score;
                        best_node = tit->second;
                        query_strand = Fasta_entry::forward_strand;
                    }

                    if(compare_reverse)
                    {
                        Fasta_entry rev_seq = *read;
                        rev_seq.sequence = this->reverse_complement(rev_seq.sequence);

                        double score = this->read_match_score( nit->second, &rev_seq, mf);

                        stringstream ss;
                        ss<<tit->second<<"(rc) with score "<<score<<"\n";
                        Log_output::write_out(ss.str(),2);

                        if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                        {
                            best_score = score;
                            best_node.append(" "+tit->second);
                            query_strand = Fasta_entry::reverse_strand;
                        }
                        else if(score>=best_score)
                        {
                            best_score = score;
                            best_node = tit->second;
                            query_strand = Fasta_entry::reverse_strand;
                        }
                    }

                    tit++;
                    matching_nodes--;
                }
            }

            if(best_score<0.05)
            {
                if(Settings_handle::st.is("align-bad-reads-at-root"))
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                    read->node_to_align = root->get_name();

                }
                else
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                    read->node_to_align = "discarded_read";
                }
            }
            else
            {
                stringstream ss;
                ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                Log_output::write_out(ss.str(),2);

                read->node_score = best_score;
                read->node_to_align = best_node;
                read->query_strand = query_strand;
//                cout<<"Read "<<read->name<<" strand "<<read->query_strand<<"\n";
            }
        }
        // All done; continue

    }


    // no TID, aligning at root

    else
    {
        stringstream ss;
        ss<<"Read "<<read->name<<" has no tid. Aligned to root.\n";
        Log_output::write_out(ss.str(),2);

        read->node_to_align = root->get_name();
    }
}

/**********************************************************************/

#ifdef NCBI_TOOLKIT

void Reads_aligner::select_node_for_query(Node *root, vector<Alignment_set> *targets, Fasta_entry *query, Model_factory *mf, mol_type molecyle)
{

    double best_score = -HUGE_VAL;
    string best_node = root->get_name();
    vector<Alignment_set> best_targets;

    map<string,Node*> nodes;
    root->get_all_nodes(&nodes);

    stringstream ss;
    ss<<"Read "<<query->name<<" has "<<targets->size()<<" target nodes.\n";
    Log_output::write_out(ss.str(),1); //2);

    for(int i=0;i<(int)targets->size();i++)
    {
        map<string,Node*>::iterator nit = nodes.find(targets->at(i).seq2_id);
        double score = this->query_match_score(nit->second, query, mf);

        stringstream ss;
        ss<<"matches "<<nit->first<<" with score "<<score<<"\n";
        Log_output::write_out(ss.str(),1); //2);

        if(score==best_score && !Settings_handle::st.is("one-placement-only"))
        {
            best_score = score;
            best_node.append(" "+nit->first);
            best_targets.push_back(targets->at(i));
        }
        else if(score>=best_score)
        {
            best_score = score;
            best_node = nit->first;

            best_targets.clear();
            best_targets.push_back(targets->at(i));
        }
    }

    if(best_score<0.05)
    {
        Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",1); //2);
        query->node_score = -1.0;
        query->node_to_align = "discarded_read";
    }
    else
    {
        stringstream ss;
        ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
        Log_output::write_out(ss.str(),1); //2);

        query->node_score = best_score;
        query->node_to_align = best_node;

        targets->clear();
        targets->insert(targets->begin(),best_targets.begin(),best_targets.end());
    }
}

#endif

/**********************************************************************/

double Reads_aligner::query_match_score(Node *node, Fasta_entry *query, Model_factory *mf)
{


    double org_dist = node->get_distance_to_parent();
    node->set_distance_to_parent(0.001);

    Node * query_node = new Node();
    double r_dist = Settings_handle::st.get("query-distance").as<float>();
    query_node->set_distance_to_parent(r_dist);

    bool revcomp = query->query_strand==Fasta_entry::reverse_strand;
    this->copy_node_details(query_node,query,revcomp);

    Node * tmpnode = new Node();
    tmpnode->set_name("(tmp)");

    tmpnode->add_left_child(query_node);
    tmpnode->add_right_child(node);

    tmpnode->align_sequences_this_node(mf,true);

    node->set_distance_to_parent(org_dist);

    double score = 0;

    if(tmpnode->node_has_sequence_object)
    {

        // For scoring (below)
        Evol_model model = mf->alignment_model(r_dist+0.001);

        int matching = 0;
        int aligned = 0;

        float subst_score = 0;
        float max_subst_score_l = 0;
        float max_subst_score_r = 0;

        for( int k=1; k < tmpnode->get_sequence()->sites_length()-1; k++ )
        {
            Site *site = tmpnode->get_sequence()->get_site_at(k);

            if(site->get_children()->right_index>=0 && site->get_children()->left_index>=0)
            {

                Site *site1 = tmpnode->get_left_child()->get_sequence()->get_site_at(site->get_children()->left_index);
                Site *site2 = tmpnode->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);

                if(site1->get_state() == site2->get_state())
                    matching++;

                subst_score += model.score(site1->get_state(),site2->get_state());
                max_subst_score_l += model.score(site2->get_state(),site2->get_state());

                aligned++;
            }

            if(site->get_children()->left_index>=0)
            {
                Site *site1 = tmpnode->get_left_child()->get_sequence()->get_site_at(site->get_children()->left_index);
                max_subst_score_r += model.score(site1->get_state(),site1->get_state());
            }

        }

        double score_s = (double) matching/ (double) query_node->get_sequence()->sites_length();
        double score_l = (double) subst_score/ (double) max_subst_score_l;
        double score_r = (double) subst_score/ (double) max_subst_score_r;

        score = score_r;

        if(Settings_handle::st.is("use-identity-score"))
            score = score_s;
        else if(Settings_handle::st.is("use-target-normalised-score"))
            score = score_l;
    }

    tmpnode->has_right_child(false);
    delete tmpnode;

    return score;
}

/**********************************************************************/


void Reads_aligner::find_targets_for_queries(Node *root, vector<Fasta_entry> *reads, Model_factory *mf,map<string,string> *target_sequences, bool is_dna)
{
    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }


    // Get the original target nodes
    //
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;



    // Handle one read at time
    //
    for(int i=0;i<(int)reads->size();i++)
    {

        if(reads->at(i).node_to_align == "" || reads->at(i).node_to_align == "discarded_read")
        {
            Log_output::write_warning("Query "+reads->at(i).name+" has no match. Skipping.",1);
            continue;
        }

        reads->at(i).node_score = -1.0;
        reads->at(i).query_strand = Fasta_entry::unknown_strand;

        string tid = reads->at(i).tid;
        if( tid=="" )
            tid = "<empty>";

        map<string,string> targets_for_this;

        stringstream str(reads->at(i).node_to_align);
        string nodename;
        while(str >> nodename)
        {
            map<string,string>::iterator it = target_sequences->find(nodename);
            if(it!=target_sequences->end())
            {
                targets_for_this.insert(make_pair(it->first,it->second));
            }
        }


        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping query: '"+reads->at(i).name+" "+reads->at(i).comment+"'",0);

        // Call Exonerate to reduce the search space
        //
        map<string,hit> exonerate_hits;

        if(has_exonerate && Settings::exonerate_gapped_keep_best > 0 )
        {
            er.local_alignment(&targets_for_this,&reads->at(i),&exonerate_hits,false,is_dna);
        }


        // Handle (or reject) reads discarded by Exonerate
        //
        if(reads->at(i).node_to_align == "discarded_read" || exonerate_hits.size()==0)
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
            Log_output::write_out(ss.str(),2);

            continue;
        }

        // Has TID or exhaustive search: now find the one with the best match
        //
        int matches = exonerate_hits.size();


        // Has only one matching node
        //
        if(matches == 1 )
        {
            string target_node = exonerate_hits.begin()->first;

            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<target_node<<"."<<endl;
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = target_node;
        }


        // Has several matching nodes
        //
        else
        {
            double best_score = -HUGE_VAL;
            string best_node = root->get_name();
            int query_strand = Fasta_entry::forward_strand;

            map<string,Node*> nodes;
            root->get_all_nodes(&nodes);

            int matching_nodes = exonerate_hits.size();

            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
            Log_output::write_out(ss.str(),2);

            map<string,hit>::iterator it = exonerate_hits.begin();

            for(;it != exonerate_hits.end(); it++)
            {
                string target_node = it->first;

                map<string,Node*>::iterator nit = nodes.find(target_node);
                double score = this->read_match_score( nit->second, &reads->at(i), mf);

                stringstream ss;
                ss<<target_node<<" with score "<<score<<"\n";
                Log_output::write_out(ss.str(),2);

                if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                {
                    best_score = score;
                    best_node.append(" "+target_node);
                }
                else if(score>=best_score)
                {
                    best_score = score;
                    best_node = target_node;
                }

                if(compare_reverse)
                {
                    Fasta_entry rev_seq = reads->at(i);
                    rev_seq.sequence = this->reverse_complement(rev_seq.sequence);

                    double score = this->read_match_score( nit->second, &rev_seq, mf);

                    stringstream ss;
                    ss<<target_node<<"(rc) with score "<<score<<"\n";
                    Log_output::write_out(ss.str(),2);

                    if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                    {
                        best_score = score;
                        best_node.append(" "+target_node);
                        query_strand = Fasta_entry::reverse_strand;
                    }
                    else if(score>=best_score)
                    {
                        best_score = score;
                        best_node = target_node;
                        query_strand = Fasta_entry::reverse_strand;
                    }
                }

                matching_nodes--;
            }

            if(best_score<0.05)
            {
                Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);
                reads->at(i).node_to_align = "discarded_read";
            }
            else
            {
                stringstream ss;
                ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_score = best_score;
                reads->at(i).node_to_align = best_node;
                reads->at(i).query_strand = query_strand;
            }
        }
        // All done; continue

    }
}

/**********************************************************************/

void Reads_aligner::find_nodes_for_queries(Node *root, vector<Fasta_entry> *reads, Model_factory *mf,bool is_dna)
{
    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }


    // Get the original target nodes
    //
    multimap<string,string> tid_nodes;
    bool ignore_tid_tags = true;
    bool compare_reverse = ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna;

    bool preselected_nodes = false;
    for(int i=0;i<(int)reads->size();i++)
    {
        if(reads->at(i).node_to_align != "")
        {
            stringstream str(reads->at(i).node_to_align);
            string nodename;
            while(str >> nodename)
            {
                tid_nodes.insert(make_pair(nodename,nodename));
            }
            preselected_nodes = true;
        }
    }
    if(!preselected_nodes)
    {
        this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
    }

    if( !Settings_handle::st.is("all-nodes") &&
        !Settings_handle::st.is("internal-nodes") &&
        !Settings_handle::st.is("terminal-nodes") &&
        !Settings_handle::st.is("test-every-internal-node") &&
        !Settings_handle::st.is("test-every-terminal-node") &&
        !Settings_handle::st.is("test-every-node") && ignore_tid_tags)
    {
        Log_output::write_warning("No tagged nodes found. Considering all nodes!",0);
    }

    // Handle one read at time
    //
    for(int i=0;i<(int)reads->size();i++)
    {
        reads->at(i).node_score = -1.0;
        reads->at(i).query_strand = Fasta_entry::unknown_strand;

        string tid = reads->at(i).tid;
        if( ignore_tid_tags )
            tid = "<empty>";

        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping query: '"+reads->at(i).name+" "+reads->at(i).comment+"'",0);

        // Call Exonerate to reduce the search space
        //
        map<string,hit> exonerate_hits;

        if(has_exonerate && Settings::exonerate_local_keep_best > 0  && !preselected_nodes)
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits, true,is_dna,ignore_tid_tags);

            if(tid_nodes.size()>1 && Settings::exonerate_gapped_keep_best > 0 )
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,is_dna,ignore_tid_tags);
        }
        else if(has_exonerate && Settings::exonerate_gapped_keep_best > 0 )
        {
            if(!preselected_nodes)
            {
                tid_nodes.clear();
                this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
            }
            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,is_dna,ignore_tid_tags);
        }

//        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping query: '"+reads->at(i).name+" "+reads->at(i).comment+"'",0);

        // Handle (or reject) reads discarded by Exonerate
        //
        if(reads->at(i).node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("exhaustive-placement"))
            {
                reads->at(i).node_to_align = "";
                tid_nodes.clear();

                this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
            }
            else
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
                Log_output::write_out(ss.str(),2);

                continue;
            }
        }


        // Has TID or exhaustive search: now find the one with the best match
        //
        if(tid != "")
        {
            int matches = tid_nodes.count(tid);

            if( ignore_tid_tags )
                matches = tid_nodes.size();


            // Has TID but no matching node
            //
            if(matches == 0)
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_to_align = root->get_name();
            }
            // All done; continue


            // Has only one matching node and no need for ranking
            //
            else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
            {
                multimap<string,string>::iterator tit = tid_nodes.find(tid);

                if( ignore_tid_tags )
                    tit = tid_nodes.begin();

                if(tit != tid_nodes.end())
                {
                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_to_align = tit->second;
                }
            }
            // All done; continue


            // Has TID and matching nodes, or exhaustive search
            //
            else
            {
                double best_score = -HUGE_VAL;
                string best_node = root->get_name();
                int query_strand = Fasta_entry::forward_strand;

                map<string,Node*> nodes;
                root->get_all_nodes(&nodes);

                multimap<string,string>::iterator tit;
                int matching_nodes = 0;

                if( ignore_tid_tags )
                {
                    tit = tid_nodes.begin();
                    matching_nodes = tid_nodes.size();
                }
                else
                {
                    tit = tid_nodes.find(tid);
                    matching_nodes = tid_nodes.count(tid);
                }

                if(tit != tid_nodes.end())
                {

                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                    Log_output::write_out(ss.str(),2);

                    while(tit != tid_nodes.end() && matching_nodes>0)
                    {
                        map<string,Node*>::iterator nit = nodes.find(tit->second);
                        double score = this->read_match_score( nit->second, &reads->at(i), mf);

                        stringstream ss;
                        ss<<tit->second<<" with score "<<score<<"\n";
                        Log_output::write_out(ss.str(),2);

                        if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                        {
                            best_score = score;
                            best_node.append(" "+tit->second);
                        }
                        else if(score>=best_score)
                        {
                            best_score = score;
                            best_node = tit->second;
                        }

                        if(compare_reverse)
                        {
                            Fasta_entry rev_seq = reads->at(i);
                            rev_seq.sequence = this->reverse_complement(rev_seq.sequence);

                            double score = this->read_match_score( nit->second, &rev_seq, mf);

                            stringstream ss;
                            ss<<tit->second<<"(rc) with score "<<score<<"\n";
                            Log_output::write_out(ss.str(),2);

                            if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                            {
                                best_score = score;
                                best_node.append(" "+tit->second);
                                query_strand = Fasta_entry::reverse_strand;
                            }
                            else if(score>=best_score)
                            {
                                best_score = score;
                                best_node = tit->second;
                                query_strand = Fasta_entry::reverse_strand;
                            }
                        }

                        tit++;
                        matching_nodes--;
                    }
                }

                if(best_score<0.05)
                {
                    if(Settings_handle::st.is("align-bad-reads-at-root"))
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                        reads->at(i).node_to_align = root->get_name();

                    }
                    else
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                        reads->at(i).node_to_align = "discarded_read";
                    }
                }
                else
                {
                    stringstream ss;
                    ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_score = best_score;
                    reads->at(i).node_to_align = best_node;
                    reads->at(i).query_strand = query_strand;
                }
            }
            // All done; continue

        }


        // no TID, aligning at root

        else
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") has no tid. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = root->get_name();

        }
    }
}

void Reads_aligner::preselect_target_sequences(Node *root, vector<Fasta_entry> *reads, map<string,string> *target_sequences, bool is_dna)
{
    if(!Settings::placement_preselection)
        return;


    Log_output::write_header("Preselecting target sequences",0);

    // Get target node sequences for all nodes; done only once
    //
    set<string> node_names;

    if(Settings::placement_target_nodes == Settings::tid_nodes)
    {
        root->get_node_names_with_tid_tag(&node_names);
        if((int)node_names.size()==0)
        {
            Settings::placement_target_nodes = Settings::all_nodes;
            Log_output::write_warning("No TID tags found. Considering all nodes for placement.",0);
        }
    }

    if(Settings::placement_target_nodes == Settings::internal_nodes)
        root->get_internal_node_names(&node_names);

    else if(Settings::placement_target_nodes == Settings::all_nodes)
        root->get_node_names(&node_names);


    bool is_now_dna = is_dna;

    map<string,string> unaligned_sequences;

    if(Settings::placement_target_nodes == Settings::terminal_nodes)
    {
        if(Settings_handle::st.is("translate") && Settings_handle::st.is("score-as-dna"))
        {
            root->get_dna_sequences(&unaligned_sequences);
            is_now_dna = true;
        }
        else
            root->get_unaligned_sequences(&unaligned_sequences);
    }
    else
    {
        if(Settings_handle::st.is("translate") && Settings_handle::st.is("score-as-dna"))
        {
            Log_output::write_warning("combination '--translate' and '--score-as-dna'' can only be used with option '--terminal-nodes'!",0);
            root->get_dna_sequences(&unaligned_sequences);
            is_now_dna = true;

            Settings::placement_target_nodes = Settings::terminal_nodes;
        }
        else
        {
            vector<Fasta_entry> aligned_sequences;
            root->get_alignment(&aligned_sequences,true);

            vector<Fasta_entry>::iterator it = aligned_sequences.begin();
            for(;it!=aligned_sequences.end();it++)
            {
                if(node_names.find(it->name) != node_names.end())
                {
                    string seq = it->sequence;
                    for (string::iterator si = seq.begin();si != seq.end();)
                        if(*si == '-')
                            seq.erase(si);
                        else
                            si++;

                    unaligned_sequences.insert(make_pair(it->name,seq));
                }
            }
        }
    }

    //    multimap<string,string> all_tid_nodes;

//TÄSTÄ ASTI
    // Call Exonerate to reduce the number of target nodes
    //
    map<string, multimap<string,hit> > exonerate_hits;

    Exonerate_queries er;
    if(!er.test_executable())
        Log_output::write_out("The executable for Exonerate not found! Preselection not done!",0);

    else if(unaligned_sequences.size()>0)
        er.preselect_targets(&unaligned_sequences,reads,target_sequences,&exonerate_hits,is_now_dna);

    set<string> keep_nodes;

    multimap<string, string> query_target;

    for(map<string, multimap<string,hit> >::iterator it = exonerate_hits.begin();it != exonerate_hits.end();it++)
    {
        Log_output::write_out("Preselect_targets: "+it->first+" has "+Log_output::itos(it->second.size())+" hits\n",2);

        multimap<string,hit>::iterator it2 = it->second.begin();

        for( ;it2 != it->second.end(); it2++ )
        {
            Log_output::write_out("  "+it->first+" matches "+it2->first+" with score "+Log_output::itos(it2->second.score)+"\n",2);
            keep_nodes.insert(it2->first);
            query_target.insert(make_pair(it2->second.query,it2->second.node));
            ///cout << "query target++: query: " << it2->second.query << ", node: " << it2->second.node << endl;
        }

    }
    Log_output::write_out("Preselect_targets: keeping "+Text_utils::to_string(keep_nodes.size())+" targets\n",1);
//TÄHÄN ASTI KAI

    for(set<string>::iterator it = keep_nodes.begin();it != keep_nodes.end();it++)
    {
        Log_output::write_out(" keep "+*it+"\n",2);
    }

    root->set_node_names_for_exonerate(&keep_nodes);


    for(vector<Fasta_entry>::iterator it = reads->begin();it!=reads->end();it++)
    {
        pair <multimap<string,string>::iterator, multimap<string,string>::iterator> it1 = query_target.equal_range(it->name);

        stringstream node_to_align;
        for (multimap<string,string>::iterator it2=it1.first; it2!=it1.second; it2++)
        {
            node_to_align <<it2->second<<" ";
        }
        it->node_to_align = node_to_align.str();
    }
}

#ifdef NCBI_TOOLKIT

void Reads_aligner::preselect_target_sequences_ncbi(Node *root, Fasta_entry *query, int num_ncbi_targets, mol_type molecyle, Model_factory *mf, string *tid)
{

    Log_output::write_header("Preselecting target sequences",0);

    map<string,string> target_sequences;

    if(Settings::placement_target_nodes == Settings::terminal_nodes)
    {
        root->get_unaligned_sequences(&target_sequences);
    }
    else
    {
        set<string> node_names;

        if(Settings::placement_target_nodes == Settings::tid_nodes)
        {
            root->get_node_names_with_tid_exact(&node_names, tid);
            if((int)node_names.size()==0)
            {
                Settings::placement_target_nodes = Settings::all_nodes;
                Log_output::write_warning("No TID tags found. Considering all nodes for placement.",0);
            }
        }

        if(Settings::placement_target_nodes == Settings::internal_nodes)
            root->get_internal_node_names(&node_names);

        else if(Settings::placement_target_nodes == Settings::all_nodes)
            root->get_node_names(&node_names);

        vector<Fasta_entry> aligned_sequences;
        root->get_alignment(&aligned_sequences,true);

        vector<Fasta_entry>::iterator it = aligned_sequences.begin();
        for(;it!=aligned_sequences.end();it++)
        {
            if(node_names.find(it->name) != node_names.end())
            {
                string seq = it->sequence;
                for (string::iterator si = seq.begin();si != seq.end();)
                    if(*si == '-')
                        seq.erase(si);
                    else
                        si++;

                target_sequences.insert(make_pair(it->name,seq));
            }
        }
    }

    map<string,string> query_sequences;

    query_sequences.insert(make_pair(query->name,query->sequence));


    /// Run Blast 'query_sequences' vs. 'target_sequences', record best hits
    BankBlaster blaster;
    Blast_options opt;

    if( ( Settings_handle::st.is("both-strands") )
                        && mf->get_sequence_data_type()==Model_factory::dna )
    {
        opt.strand_opt = strand_both;
    }

    opt.wordsize = Settings_handle::st.get("blast-wordsize").as<int>();
    if(molecyle != aminoa){
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


    blaster.initBankAlignments(query_sequences, target_sequences, num_ncbi_targets, molecyle, &opt);


}

#endif

void Reads_aligner::find_orfs(Fasta_entry *read,vector<Orf> *open_frames)
{

    int min_orf_length = Settings_handle::st.get("min-orf-length").as<int>();
    if(Settings_handle::st.is("min-orf-coverage"))
        min_orf_length = int (Settings_handle::st.get("min-orf-coverage").as<float>() * read->dna_sequence.length() / 3);


    string dna = read->dna_sequence;
    int length = dna.length()-1;

    if(length/3<min_orf_length)
    {
        Log_output::write_warning("Warning: query sequence '"+read->name+"' is shorter than the minimum ORF length.",0);
        return;
    }

    for(int i=0;i<3;i++)
    {
        string prot;
        string sequence = dna.substr(i);
        int start_site = i;
        int end_site = i+2;

        for (unsigned int j=0; j<sequence.length(); j+=3)
        {
            string codon = sequence.substr(j,3);
            if (codon_to_aa.find(codon) == codon_to_aa.end()
                    || codon_to_aa.find(codon)->second == "X" || codon_to_aa.find(codon)->second == "-" )
            {
                if((int)prot.length() >= min_orf_length)
                {
                    Orf o;
                    o.translation = prot;
                    o.frame = i+1;
                    o.start = start_site;
                    o.end = end_site;
                    o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

                    open_frames->push_back(o);
                }

                prot = "";
                start_site = j+i+3;
            }
            else
            {
                prot += codon_to_aa.find(codon)->second;
            }
            end_site = j+i+2;
        }


        if((int)prot.length() >= min_orf_length)
        {
            Orf o;
            o.translation = prot;
            o.frame = i+1;
            o.start = start_site;
            o.end = end_site;
            o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

            open_frames->push_back(o);
        }
    }

    dna = this->reverse_complement(read->dna_sequence);
    for(int i=0;i<3;i++)
    {
        string prot;
        string sequence = dna.substr(i);
        int start_site = i;
        int end_site = i+2;

        for (unsigned int j=0; j<sequence.length(); j+=3)
        {
            string codon = sequence.substr(j,3);
            if (codon_to_aa.find(codon) == codon_to_aa.end()
                    || codon_to_aa.find(codon)->second == "X" || codon_to_aa.find(codon)->second == "-" )
            {
                if((int)prot.length() >= min_orf_length)
                {
                    Orf o;
                    o.translation = prot;
                    o.frame = -1*(i+1);
                    o.start = length - end_site;
                    o.end = length - start_site;
                    o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

                    open_frames->push_back(o);
                }

                prot = "";
                start_site = j+i+3;
            }
            else
            {
                prot += codon_to_aa.find(codon)->second;
            }
            end_site = j+i+2;
        }

        if((int)prot.length() >= min_orf_length)
        {
            Orf o;
            o.translation = prot;
            o.frame = -1*(i+1);
            o.start = length - end_site;
            o.end = length - start_site;
            o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

            open_frames->push_back(o);
        }
    }

    Log_output::write_msg("ORF search: '"+read->name+"' has "+Text_utils::to_string(open_frames->size())+" ORF(s)",1);
}

string Reads_aligner::reverse_complement(string dna)
{
    string rev = dna;
    reverse(rev.begin(), rev.end());

    Text_utils tu;

    tu.replace_all(rev,'A','Z');
    tu.replace_all(rev,'T','A');
    tu.replace_all(rev,'Z','T');
    tu.replace_all(rev,'C','Z');
    tu.replace_all(rev,'G','C');
    tu.replace_all(rev,'Z','G');
    tu.replace_all(rev,'R','Z');
    tu.replace_all(rev,'Y','R');
    tu.replace_all(rev,'Z','Y');
    tu.replace_all(rev,'K','Z');
    tu.replace_all(rev,'M','K');
    tu.replace_all(rev,'Z','M');
    tu.replace_all(rev,'B','Z');
    tu.replace_all(rev,'V','B');
    tu.replace_all(rev,'Z','V');
    tu.replace_all(rev,'D','Z');
    tu.replace_all(rev,'H','D');
    tu.replace_all(rev,'Z','H');

    return rev;
}

void Reads_aligner::define_translation_tables()
{
    string codon[66] = {"TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                        "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                        "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                        "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                        "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
                        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
                        "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
                        "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
                        "NNN", "---"
                       };
    string unaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "X", "W", "R", "R", "R", "R", "S", "S", "R", "R", "G", "G", "G", "G",
                       "X", "-"
                      };
    string mtaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "M", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "W", "W", "R", "R", "R", "R", "S", "S", "X", "X", "G", "G", "G", "G",
                       "X", "-"
                      };

    if(Settings_handle::st.is("mt-translate"))
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],mtaa[i]));
            aa_to_codon.insert(make_pair(mtaa[i],codon[i]));
        }
    }
    else
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],unaa[i]));
            aa_to_codon.insert(make_pair(unaa[i],codon[i]));
        }
    }
}



/**********************************************************************/

void Reads_aligner::read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity)
{
    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;
    int matched = 0;

    string ref_dna_string = "";
    string read_dna_string = "";
    bool as_dna = false;
    int ref_pos = 0;
    int read_pos = 0;
    int step = 1;

    if(Settings_handle::st.is("score-as-dna"))
    {
        ref_dna_string = *(node->get_dna_sequence_for_node(ref_node_name));
        read_dna_string = *(node->get_dna_sequence_for_node(read_name));

        Text_utils tu;
        tu.replace_all(ref_dna_string,"-","");
//       cout<<endl<<ref_node_name<<" "<<ref_dna_string<<endl<<read_name<<" "<<read_dna_string<<endl;

        if((int)ref_dna_string.length()>0 && (int)read_dna_string.length()>0)
        {
            as_dna = true;
            step = 3;
        }
    }

    if(Settings_handle::st.is("overlap-with-any"))
    {
        for( int j=1; j < node_sequence->sites_length(); j++ )
        {
            bool read_has_site = node->has_site_at_alignment_column(j,read_name);
            bool any_other_has_site = node->any_other_has_site_at_alignment_column(j,read_name);
            bool ref_root_has_site = node->has_site_at_alignment_column(j,ref_node_name);

            if(read_has_site && any_other_has_site)
            {

                int state_read = node->get_state_at_alignment_column(j,read_name);
                int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);

                if(state_read>=0 && state_read == state_ref)
                {
                    if(as_dna)
                    {
                        if(ref_pos+3 <= (int)ref_dna_string.length() && read_pos+3 <= (int)read_dna_string.length())
                        {
                            if(ref_dna_string.at(ref_pos) == read_dna_string.at(read_pos))
                                matched++;
                            if(ref_dna_string.at(ref_pos+1) == read_dna_string.at(read_pos+1))
                                matched++;
                            if(ref_dna_string.at(ref_pos+2) == read_dna_string.at(read_pos+2))
                                matched++;
                        }
                    }
                    else
                    {
                        matched++;
                    }
                }
                if(as_dna)
                    aligned += step;
                else
                    aligned++;
            }

            if(read_has_site)
            {
                read_length += step;
                if(as_dna)
                    read_pos += step;
            }

            if(ref_root_has_site && as_dna)
                    ref_pos += step;

        }
    }

    else
    {
        for( int j=1; j < node_sequence->sites_length(); j++ )
        {
            bool read_has_site = node->has_site_at_alignment_column(j,read_name);
            bool ref_root_has_site = node->has_site_at_alignment_column(j,ref_node_name);


            if(read_has_site && ref_root_has_site)
            {
                int state_read = node->get_state_at_alignment_column(j,read_name);
                int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);

                if(state_read>=0 && state_read == state_ref)
                {
                    if(as_dna)
                    {
                        if(ref_pos+3 <= (int)ref_dna_string.length() && read_pos+3 <= (int)read_dna_string.length())
                        {
                            if(ref_dna_string.at(ref_pos) == read_dna_string.at(read_pos))
                                matched++;
                            if(ref_dna_string.at(ref_pos+1) == read_dna_string.at(read_pos+1))
                                matched++;
                            if(ref_dna_string.at(ref_pos+2) == read_dna_string.at(read_pos+2))
                                matched++;
                        }
                    }
                    else
                    {
                        matched++;
                    }
                }
                if(as_dna)
                    aligned += step;
                else
                    aligned++;
            }

            if(read_has_site)
            {
                read_length += step;
                if(as_dna)
                    read_pos += step;
            }

            if(ref_root_has_site && as_dna)
                    ref_pos += step;
        }
    }

    stringstream ss;
    ss<<"aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"]"<<endl;
    Log_output::write_out(ss.str(),2);
    ss.str(string());
    ss<<"matched positions "<<(float)matched/(float)aligned<<" ["<<matched<<"/"<<aligned<<"]"<<endl;
    Log_output::write_out(ss.str(),2);

    *overlap  = (float)aligned/(float)read_length;
    *identity = (float)matched/(float)aligned;
}

double Reads_aligner::read_match_score(Node *node, Fasta_entry *read, Model_factory *mf)
{

    double r_dist = Settings_handle::st.get("query-distance").as<float>();

    double org_dist = node->get_distance_to_parent();
    node->set_distance_to_parent(0.001);

    Node * reads_node1 = new Node();
    reads_node1->set_distance_to_parent(r_dist);
    reads_node1->set_name(read->name);
    reads_node1->add_name_comment(read->comment);
    reads_node1->add_sequence( *read, read->data_type, false, true);

    Node * tmpnode = new Node();
    tmpnode->set_name("(tmp)");

    tmpnode->add_left_child(node);
    tmpnode->add_right_child(reads_node1);

    tmpnode->align_sequences_this_node(mf,true);

    node->set_distance_to_parent(org_dist);

    double score = 0;

    if(tmpnode->node_has_sequence_object)
    {

        // For scoring (below)
        Evol_model model = mf->alignment_model(r_dist+0.001);

        int matching = 0;
        int aligned = 0;

        float subst_score = 0;
        float max_subst_score_l = 0;
        float max_subst_score_r = 0;

        for( int k=1; k < tmpnode->get_sequence()->sites_length()-1; k++ )
        {
            Site *site = tmpnode->get_sequence()->get_site_at(k);

            if(site->get_children()->right_index>=0 && site->get_children()->left_index>=0)
            {

                Site *site1 = tmpnode->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);
                Site *site2 = tmpnode->get_left_child()->get_sequence()->get_site_at(site->get_children()->left_index);

                if(site1->get_state() == site2->get_state())
                    matching++;

                subst_score += model.score(site1->get_state(),site2->get_state());
                max_subst_score_l += model.score(site2->get_state(),site2->get_state());

                aligned++;
            }

            if(site->get_children()->right_index>=0)
            {
                Site *site1 = tmpnode->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);
                max_subst_score_r += model.score(site1->get_state(),site1->get_state());
            }

        }

        double score_s = (double) matching/ (double) reads_node1->get_sequence()->sites_length();
        double score_l = (double) subst_score/ (double) max_subst_score_l;
        double score_r = (double) subst_score/ (double) max_subst_score_r;

        score = score_r;

        if(Settings_handle::st.is("use-identity-score"))
            score = score_s;
        else if(Settings_handle::st.is("use-target-normalised-score"))
            score = score_l;
    }

    tmpnode->has_left_child(false);
    delete tmpnode;

    return score;
}


bool Reads_aligner::correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map)
{

    // correct the sites index at the parent node; insertions corrected later
    //
    vector<int> sites_index;
    int index_delta = 0;

    for(int j=0; j<current_root->get_sequence()->sites_length(); j++)
    {
        if(current_root->has_site_at_alignment_column(j,ref_node_name))
        {
            sites_index.push_back(index_delta);
            index_delta = 0;
        }
        else
            index_delta++;
    }


    Node *current_parent = 0;
    map<string,Node*>::iterator mit = nodes_map->begin();
    bool parent_found = false;

    int is_left_child = true;
    for(;mit != nodes_map->end();mit++)
    {
//        if(!mit->second->is_leaf())
//            cout<<mit->second->get_name()<<" "<<mit->second->get_left_child()->get_name()<<" "<<ref_node_name<<endl;
        if(!mit->second->is_leaf() && mit->second->get_left_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_left_child(current_root);
            parent_found = true;
        }

//        if(!mit->second->is_leaf())
//            cout<<mit->second->get_name()<<" "<<mit->second->get_right_child()->get_name()<<" "<<ref_node_name<<endl;
        if(!mit->second->is_leaf() && mit->second->get_right_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_right_child(current_root);
            is_left_child = false;
            parent_found = true;
        }
    }


    if(parent_found)
    {
        stringstream ss;
        ss<<"Parent of "<<ref_node_name<<" is "<<current_parent->get_name()<<"; "
          <<alignments_done<<" alignments done.";
        Log_output::write_out(ss.str(),3);

        Sequence *parent_sequence = current_parent->get_sequence();

        index_delta = 0;
        int first=0;
        for(int j=1; j<parent_sequence->sites_length(); j++)
        {
            Site *parent_site = parent_sequence->get_site_at(j);

            if(is_left_child && parent_site->get_children()->left_index > 0)
            {
                first = parent_site->get_children()->left_index;
                break;
            }
            else if(!is_left_child && parent_site->get_children()->right_index > 0)
            {
                first = parent_site->get_children()->right_index;
                break;
            }
        }

        for(int j=0;j<first;j++)
            index_delta += sites_index.at(j);

        for(int j=1; j<parent_sequence->sites_length(); j++)
        {
            Site *parent_site = parent_sequence->get_site_at(j);

            if(is_left_child && parent_site->get_children()->left_index >= 0)
            {
                index_delta += sites_index.at(parent_site->get_children()->left_index);
                parent_site->get_children()->left_index += index_delta;
            }
            else if(!is_left_child && parent_site->get_children()->right_index >= 0)
            {
                index_delta += sites_index.at(parent_site->get_children()->right_index);
                parent_site->get_children()->right_index += index_delta;
            }
        }

        if(index_delta>0)
            Log_output::write_out("Site index needs correcting.\n",3);
        else
            Log_output::write_out("Site index not changed.\n",3);


        if(index_delta>0)
        {
            if(is_left_child)
                current_parent->left_needs_correcting_sequence_site_index(true);
            else
                current_parent->right_needs_correcting_sequence_site_index(true);
        }

        return true;

    } // if(parent_found)
    else
    {
        Log_output::write_out("No parent for "+ref_node_name+" found. Assuming that this is root.\n",2);

        return false;
    }

}


/************************************************************************************************/


void Reads_aligner::do_upwards_search(Node *root, Fasta_entry *read, Model_factory *mf)
{

    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    double r_dist = Settings_handle::st.get("query-distance").as<float>();
    Evol_model model = mf->alignment_model(r_dist+0.001);

    Node * current_root = root;
    string previous_hit = "";

    while(true)
    {
        Log_output::write_msg(" aligning read: '"+read->name+"'",0);

        Node * temp_node = new Node();
        temp_node->set_name("temp");

        Node * left_node = current_root;
        float org_distance = left_node->get_distance_to_parent();
        left_node->set_distance_to_parent(0.001);
        temp_node->add_left_child(left_node);

        Node * right_node = new Node();
        this->copy_node_details(right_node,read);
        temp_node->add_right_child(right_node);

        temp_node->align_sequences_this_node(mf,true);


        vector<int> one_score;
        vector<int> all_scores;
        int read_length = 0;
        vector<float> all_dists;
        vector<float> one_dist;
        float read_max_dist = 0;

        vector<Node*> score_nodes;

        bool leaves_only = true;
        if(leaves_only)
        {
            for(int j=0;j<current_root->get_number_of_leaves();j++)
            {
                all_scores.push_back(0);
            }

            for(int j=0;j<current_root->get_number_of_leaves();j++)
            {
                all_dists.push_back(0);
            }
            current_root->get_all_nodes(&score_nodes,true);
        }
        else
        {
            for(int j=0;j<current_root->get_number_of_nodes();j++)
            {
                all_scores.push_back(0);
            }

            for(int j=0;j<current_root->get_number_of_nodes();j++)
            {
                all_dists.push_back(0);
            }
            current_root->get_leaf_nodes(&score_nodes);
        }


        for(int j=1;j<temp_node->get_sequence()->sites_length()-1;j++)
        {
            Site_children *offspring = temp_node->get_sequence()->get_site_at(j)->get_children();

            int lj = offspring->left_index;
            int rj = offspring->right_index;
            if(rj>=0)
            {
                int qs = temp_node->get_right_child()->get_sequence()->get_site_at(rj)->get_state();

                if(lj>=0)
                {
                    one_score.clear();
                    current_root->get_state_identity(lj,qs,&one_score,true,leaves_only);

                    for(int l=0;l<(int)all_scores.size();l++)
                        all_scores.at(l)+=one_score.at(l);

                    //

                    one_dist.clear();
                    current_root->get_subst_distance(lj,qs,&one_dist,&model,true,leaves_only);

                    for(int l=0;l<(int)all_scores.size();l++)
                        all_dists.at(l)+=one_dist.at(l);
                }
                read_length++;

                read_max_dist += model.score(qs,qs);
            }
        }


        double max_score = -1;
        string max_node = current_root->get_name();
        int max_node_index = 0;

        double max_dist = -1;
        string max_dist_node = current_root->get_name();
        int max_dist_node_index = 0;


        for(int l=0;l<(int)all_scores.size();l++)
        {
            float score = (float)all_scores.at(l)/(float)read_length;
            if(score>=max_score)
            {
                max_score = score;
                max_node = score_nodes.at(l)->get_name();
                max_node_index = l;
            }
            float dist = (float)all_dists.at(l)/(float)read_max_dist;
            if(dist>=max_dist)
            {
                max_dist = dist;
                max_dist_node = score_nodes.at(l)->get_name();
                max_dist_node_index = l;
            }
//            cout<<score_nodes.at(l)->get_name()<<" "<<score<<" "<<dist<<endl;
        }

        cout<<endl;
        cout<<"s1 "<<read->name<<": "<<max_node<<" ("<<max_score<<")\n";
        cout<<"s2 "<<read->name<<": "<<max_dist_node<<" ("<<max_dist<<")\n";


//            if(previous_hit=="")
////            if(true)
//            {
//                vector<Fasta_entry> aligned_seqs;
//                temp_node->get_alignment(&aligned_seqs,true);
//                ofstream out(reads->at(i).name.c_str());
//                for(int j=0;j<aligned_seqs.size();j++)
//                    out<<">"<<aligned_seqs.at(j).name<<"\n"<<aligned_seqs.at(j).sequence<<endl;
//            }

        left_node->set_distance_to_parent(org_distance);
        temp_node->has_left_child(false);
        delete temp_node;


        if(max_dist_node == previous_hit)
        {
            read->node_score = max_dist;
            read->node_to_align = max_dist_node;

            break;
        }
        else
        {
            previous_hit = max_dist_node;
            Node *temp = current_root->get_parent_node(max_dist_node);
            if(temp != 0)
            {
                current_root = temp;
//                    cout<<"new root "<<current_root->get_name()<<endl;
            }
            else
            {
                read->node_score = max_dist;
                read->node_to_align = max_dist_node;
                break;
            }
        }

    }


    cout<<"placement "<<read->name<<": "<<read->node_to_align<<" ("<<read->node_score <<")\n";
}

void Reads_aligner::do_upwards_search(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    double r_dist = Settings_handle::st.get("query-distance").as<float>();
    Evol_model model = mf->alignment_model(r_dist+0.001);


    for(int i=0;i<(int)reads->size();i++)
//    for(int i=0;i<1;i++)
    {
        Node * current_root = root;
        string previous_hit = "";

        while(true)
        {
            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: '"+reads->at(i).name+"'",0);

            Node * temp_node = new Node();
            temp_node->set_name("temp");

            Node * left_node = current_root;
            float org_distance = left_node->get_distance_to_parent();
            left_node->set_distance_to_parent(0.001);
            temp_node->add_left_child(left_node);

            Node * right_node = new Node();
            this->copy_node_details(right_node,&reads->at(i));
            temp_node->add_right_child(right_node);

            temp_node->align_sequences_this_node(mf,true);

            vector<int> one_score;
            vector<int> all_scores;
            int read_length = 0;
            vector<float> all_dists;
            vector<float> one_dist;
            float read_max_dist = 0;

            vector<Node*> score_nodes;

            bool leaves_only = true;
            if(leaves_only)
            {
                for(int j=0;j<current_root->get_number_of_leaves();j++)
                {
                    all_scores.push_back(0);
                }

                for(int j=0;j<current_root->get_number_of_leaves();j++)
                {
                    all_dists.push_back(0);
                }
                current_root->get_leaf_nodes(&score_nodes);
            }
            else
            {
                for(int j=0;j<current_root->get_number_of_nodes();j++)
                {
                    all_scores.push_back(0);
                }

                for(int j=0;j<current_root->get_number_of_nodes();j++)
                {
                    all_dists.push_back(0);
                }
                current_root->get_all_nodes(&score_nodes,true);
            }

            vector<bool> scored_sites;
            if(Settings_handle::st.is("score-only-ungapped"))
            {
                float limit = Settings_handle::st.get("score-ungapped-limit").as<float>();
                int num_leaves = temp_node->get_number_of_leaves();
                scored_sites.push_back(false);
                for(int j=1;j<temp_node->get_sequence()->sites_length()-1;j++)
                {
                    int gaps_this_site = temp_node->get_number_of_gaps_at_site(j);
                    if((float)gaps_this_site/float(num_leaves)>=limit)
                        scored_sites.push_back(false);
                    else
                        scored_sites.push_back(true);
                }
                scored_sites.push_back(false);
            }
//            for(int j=1;j<temp_node->get_sequence()->sites_length()-1;j++)
//            {
//                cout<<j<<" "<<scored_sites.at(j)<<endl;
//            }
            for(int j=1;j<temp_node->get_sequence()->sites_length()-1;j++)
            {
                if(not Settings_handle::st.is("score-only-ungapped") || scored_sites.at(j))
                {
                    Site_children *offspring = temp_node->get_sequence()->get_site_at(j)->get_children();

                    int lj = offspring->left_index;
                    int rj = offspring->right_index;
                    if(rj>=0)
                    {
                        int qs = temp_node->get_right_child()->get_sequence()->get_site_at(rj)->get_state();

                        if(lj>=0)
                        {
                            one_score.clear();
                            current_root->get_state_identity(lj,qs,&one_score,true);

                            for(int l=0;l<(int)all_scores.size();l++)
                                all_scores.at(l)+=one_score.at(l);

                            //

                            one_dist.clear();
                            current_root->get_subst_distance(lj,qs,&one_dist,&model,true);

                            for(int l=0;l<(int)all_scores.size();l++)
                                all_dists.at(l)+=one_dist.at(l);
                        }
                        read_length++;

                        read_max_dist += model.score(qs,qs);
                    }
                }
            }

            double max_score = -1;
            string max_node = current_root->get_name();
            int max_node_index = 0;

            double max_dist = -1;
            string max_dist_node = current_root->get_name();
            int max_dist_node_index = 0;

            for(int l=0;l<(int)all_scores.size();l++)
            {
                float score = (float)all_scores.at(l)/(float)read_length;
                if(score>=max_score)
                {
                    max_score = score;
                    max_node = score_nodes.at(l)->get_name();
                    max_node_index = l;
                }
                float dist = (float)all_dists.at(l)/(float)read_max_dist;
                if(dist>=max_dist)
                {
                    max_dist = dist;
                    max_dist_node = score_nodes.at(l)->get_name();
                    max_dist_node_index = l;
                }
            }

            cout<<endl;
            cout<<"s1 "<<reads->at(i).name<<": "<<max_node<<" ("<<max_score<<")\n";
            cout<<"s2 "<<reads->at(i).name<<": "<<max_dist_node<<" ("<<max_dist<<")\n";


//            if(previous_hit=="")
////            if(true)
//            {
//                vector<Fasta_entry> aligned_seqs;
//                temp_node->get_alignment(&aligned_seqs,true);
//                ofstream out(reads->at(i).name.c_str());
//                for(int j=0;j<aligned_seqs.size();j++)
//                    out<<">"<<aligned_seqs.at(j).name<<"\n"<<aligned_seqs.at(j).sequence<<endl;
//            }

            left_node->set_distance_to_parent(org_distance);
            temp_node->has_left_child(false);
            delete temp_node;


            if(max_dist_node == previous_hit)
            {
                reads->at(i).node_score = max_dist;
                reads->at(i).node_to_align = max_dist_node;

                break;
            }
            else
            {
                previous_hit = max_dist_node;
                Node *temp = current_root->get_parent_node(max_dist_node);
                if(temp != 0)
                {
                    current_root = temp;
//                    cout<<"new root "<<current_root->get_name()<<endl;
                }
                else
                {
                    reads->at(i).node_score = max_dist;
                    reads->at(i).node_to_align = max_dist_node;
                    break;
                }
            }

        }


        cout<<"placement "<<reads->at(i).name<<": "<<reads->at(i).node_to_align<<" ("<<reads->at(i).node_score <<")\n";
    }

}

/**********************************************************************/



