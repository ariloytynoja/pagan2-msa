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

#ifndef READS_ALIGNER_H
#define READS_ALIGNER_H

#include <vector>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/model_factory.h"
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "main/node.h"
#include "main/sequence.h"
#include "main/viterbi_alignment.h"

#ifdef NCBI_TOOLKIT
#include "utils/ncbi_blast.h"
#endif

using namespace std;

namespace ppa
{

struct Ncbi_hit {
    string tid;
    string qid;
    //
    //  other stuff
    //
    int score;
    bool plus_strand;
};

class Reads_aligner
{
    Node *global_root;
    map<string,string> codon_to_aa;
    map<string,string> aa_to_codon;

    void pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void translated_pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void query_placement_all(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count, bool is_dna);
    void query_placement_one(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count, bool is_dna);
    void translated_query_placement_all(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void translated_query_placement_one(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);

    void do_upwards_search(Node *root, Fasta_entry *read, Model_factory *mf);
    void do_upwards_search(Node *root, vector<Fasta_entry> *reads, Model_factory *mf);

    void fix_branch_lengths(Node *root,Node *current_root);

    void find_orfs(Fasta_entry *read,vector<Orf> *open_frames);
    void define_translation_tables();
    string reverse_complement(string dna);

    void find_nodes_for_queries(Node *root, vector<Fasta_entry> *reads, Model_factory *mf,bool is_dna);
    void find_nodes_for_query(Node *root, Fasta_entry *read, Model_factory *mf,bool is_dna,bool warnings=false);

    void find_targets_for_queries(Node *root, vector<Fasta_entry> *reads, Model_factory *mf,map<string,string> *target_sequences, bool is_dna);
    void find_targets_for_query(Node *root, Fasta_entry *read, Model_factory *mf,map<string,string> *target_sequences,multimap<string, string> *added_sequences, bool is_dna, bool warnings=false);

    double query_match_score(Node *node, Fasta_entry *read, Model_factory *mf);
    double read_match_score(Node *node, Fasta_entry *read, Model_factory *mf);

    void read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity);
    bool read_alignment_overlaps(Node * node, string read_name, string ref_node_name);

    bool correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map);
    void preselect_target_sequences(Node *root, vector<Fasta_entry> *reads, map<string,string> *target_sequences, bool is_dna);

#ifdef NCBI_TOOLKIT
    void preselect_target_sequences_ncbi(Node *root, Fasta_entry *query, int num_ncbi_targets, mol_type molecyle, Model_factory *mf, string *tid = 0);
    void select_node_for_query(Node *root, vector<Alignment_set> *targets, Fasta_entry *query, Model_factory *mf, mol_type molecyle);
    void query_placement_one_ncbi(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count, mol_type molecyle);
#endif

    static bool better_score(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.node_score>b.node_score);
    }

    void sort_reads_vector(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_aligner::better_score);
    }

    static bool longer_sequence(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.sequence.length()>b.sequence.length());
    }

    void sort_reads_vector_by_length(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_aligner::longer_sequence);
    }

    static bool more_copies(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.num_duplicates>b.num_duplicates);
    }

    void sort_reads_vector_by_duplicate_number(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_aligner::more_copies);
    }

    static bool node_is_smaller(const string& l,const string& r)
    {   char a,b,c,d;
        int vl=0; int vr=0;
        stringstream ls(l);
        stringstream rs(r);
        ls>>a>>vl>>b;
        rs>>c>>vr>>d;

        if(a=='#' && b=='#' && c=='#' && d=='#' && vl>0 && vr>0)
            return (vl<vr);

        if(a=='#' && b=='#'&& vl>0)
            return false;

        if(c=='#' && d=='#'&& vr>0)
            return true;

        return (l<r);

    }

    static bool node_is_smaller_ncbi(const Ncbi_hit& l,const Ncbi_hit& r) { return true; }

    void copy_node_details(Node *reads_node,Fasta_entry *read,bool turn_revcomp=false)
    {
        double r_dist = Settings_handle::st.get("query-distance").as<float>();

        reads_node->set_distance_to_parent(r_dist);
        reads_node->set_name(read->name);
        reads_node->add_name_comment(read->comment);
        reads_node->add_sequence( *read, read->data_type, false, true, turn_revcomp);
        reads_node->get_sequence()->is_read_sequence(true);

        if(read->dna_sequence.length()>0)
        {
            Orf o;
            o.dna_sequence = read->dna_sequence;
            o.translation = read->sequence;

            reads_node->set_Orf(&o);
        }
    }

    void create_temp_node(Node *node,string ss, Node *global_root, Fasta_entry *read,bool is_reverse)
    {
        global_root->set_distance_to_parent(0.001);

        node->set_name(ss);
        node->add_left_child(global_root);

        Node * reads_node = new Node();
        this->copy_node_details(reads_node,read,is_reverse);

        node->add_right_child(reads_node);

        node->set_nhx_tid(node->get_left_child()->get_nhx_tid());
        node->get_right_child()->set_nhx_tid(node->get_left_child()->get_nhx_tid());
    }

    void create_temp_orf_node(Node *node,Node *global_root, Fasta_entry *read, Orf *open_frame)
    {
        node->set_name("#orf#");

        global_root->set_distance_to_parent(0.001);
        node->add_left_child(global_root);

        Node * reads_node = new Node();

        Fasta_entry orf;
        orf.name = read->name;
        orf.comment = read->comment;
        orf.sequence = open_frame->translation;
        orf.dna_sequence = open_frame->dna_sequence;
        orf.quality = "";
        orf.first_read_length = -1;
        orf.tid = read->tid;
        orf.cluster_attempts = 0;
        orf.data_type = Model_factory::protein;

        this->copy_node_details(reads_node,&orf);

        node->add_right_child(reads_node);

    }

    void compute_read_overlap(Node *node,string read_name,string ref_root_name,string global_root_name,float *read_overlap,float *read_identity)
    {
        if(node->node_has_sequence_object)
        {
            if(Settings_handle::st.is("overlap-with-any"))
                this->read_alignment_scores(node, read_name,global_root_name,read_overlap,read_identity);
            else
                this->read_alignment_scores(node, read_name,ref_root_name,read_overlap,read_identity);
        }
    }

    void get_target_node_names(Node *root,multimap<string,string> *tid_nodes, bool *ignore_tid_tags, bool no_preselect = false)
    {
        if(Settings::placement_preselection && !no_preselect)
        {
            root->get_node_names_for_exonerate(tid_nodes);
        }
        else if(Settings_handle::st.is("test-every-internal-node") || Settings_handle::st.is("internal-nodes"))
        {
            root->get_internal_node_names(tid_nodes);
        }
        else if(Settings_handle::st.is("test-every-terminal-node") || Settings_handle::st.is("terminal-nodes"))
        {
            root->get_terminal_node_names(tid_nodes);
        }
        else if(Settings_handle::st.is("test-every-node") || Settings_handle::st.is("all-nodes"))
        {
            root->get_node_names(tid_nodes);
        }
        else
        {
            root->get_node_names_with_tid_tag(tid_nodes);
            if((int)tid_nodes->size()>0)
            {
                *ignore_tid_tags = false;
            }
            else
            {
//                Log_output::write_out("No tagged nodes found. Considering all nodes!\n\n",0);
                tid_nodes->clear();
                root->get_node_names(tid_nodes);
            }
        }
//        cout<<"\nget nodes "<<tid_nodes->size()<<"\n";

    }

public:
    Reads_aligner();
    void align(Node *root, Model_factory *mf,int count);

    Node *get_global_root() { return global_root; }
};
}

#endif // READS_ALIGNER_H
