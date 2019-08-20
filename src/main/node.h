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

#ifndef NODE_H
#define NODE_H

#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <string>
#include <vector>
#include "utils/exceptions.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
#include "main/viterbi_alignment.h"
#include "main/reference_alignment.h"
#include "utils/model_factory.h"
#include "utils/fasta_entry.h"
#include "utils/log_output.h"

#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>

#include <omp.h>

using namespace std;


namespace ppa
{

struct Insertion_at_node
{
    string node_name_wanted;
    int length;
    bool left_child_wanted;
    int start_site;
    int end_site;
};


struct Orf
{
    string dna_sequence;
    string translation;
    int frame;
    int start;
    int end;
};

class Node
{

    //for threading
    static boost::mutex list_mutex;
    static boost::mutex model_mutex;
    static boost::mutex log_mutex;

public:
    
    bool leaf;
    double dist_to_parent;
    string name;
    string name_comment;
    string name_id;
    string nhx_tid;
    string nhx_tag;

    Sequence *sequence;
    
    bool node_has_left_child;
    bool node_has_right_child;
    bool adjust_left_node_site_index;
    bool adjust_right_node_site_index;

    static int number_of_nodes;
    static int alignment_number;

    bool node_has_sequence_object;
    bool node_has_sequence;

    Orf orf;
    bool has_orf;
    
    Node *left_child;
    Node *right_child;

    bool visited;
    bool use_for_exonarate;

    Node() : leaf(true), dist_to_parent(0), name("undefined"), nhx_tid(""),
             node_has_left_child(false), node_has_right_child(false),
             adjust_left_node_site_index(false), adjust_right_node_site_index(false),
             node_has_sequence_object(false), node_has_sequence(false),
             has_orf(false), visited(false), use_for_exonarate(false) {}
    ~Node();


    /**************************************/

    void set_distance_to_parent(double d, bool correct = true)
    {
        dist_to_parent = d;

        if(dist_to_parent<=0 && correct)
        {
            if( Settings_handle::st.is("min-branch-length") && Settings_handle::st.get("min-branch-length").as<float>() > 0 )
            {
                dist_to_parent = Settings_handle::st.get("min-branch-length").as<float>();
            }
            else
            {
                dist_to_parent = 0.001;
            }
        }
        if(correct && !Settings_handle::st.is("real-branches") &&
                ( Settings_handle::st.is("scale-branches") || Settings_handle::st.is("truncate-branches") )  )
        {
            if( Settings_handle::st.is("scale-branches") &&
                  Settings_handle::st.get("scale-branches").as<float>() > 0 )
            {
                dist_to_parent *= Settings_handle::st.get("scale-branches").as<float>();
            }

            if( Settings_handle::st.is("real-branches") )
                ;
            else if( Settings_handle::st.is("truncate-branches") &&
                  Settings_handle::st.get("truncate-branches").as<float>() > 0 &&
                    dist_to_parent > Settings_handle::st.get("truncate-branches").as<float>() )
            {
                dist_to_parent = Settings_handle::st.get("truncate-branches").as<float>();
            }
        }
        if(correct &&  Settings_handle::st.is("fixed-branches") )
        {
            dist_to_parent = Settings_handle::st.get("fixed-branches").as<float>();
        }
    }

    double get_distance_to_parent() { return dist_to_parent; }

    /**************************************/

    void reset_alignment_number() { alignment_number = 1; }

    void set_name(string nm)
    {
        name = nm;
    }

    void set_name_ids(int *count)
    {
        if(!leaf)
        {
            left_child->set_name_ids(count);
            right_child->set_name_ids(count);
        } else {
            stringstream ss;
            ss << "seq" << (++*count);
            this->name_id = ss.str();
        }
    }

    void set_name_id(int count)
    {
        stringstream ss;
        ss << "seq" << (++count);
        this->name_id = ss.str();
    }

    void set_nhx_tid(string nm)
    {
        nhx_tid = nm;
    }

    void set_nhx_tag(string nm)
    {
        nhx_tag = nm;
    }


    void set_Orf(Orf *o) { orf = *o; has_orf = true; ;}
    Orf get_Orf() { return orf; }
    bool has_Orf() { return has_orf; }

    string get_name() const { return name; }

    string get_name_id() const { return name_id; }

    string get_nhx_tid() const { return nhx_tid; }

    string get_nhx_tag() const { return nhx_tag; }

    string get_id_for_name(string query) const
    {
        if(this->name == query)
        {
            return this->name_id;
        }
        else if(!leaf)
        {
            string tmp;
            tmp = left_child->get_id_for_name(query);
            if(tmp!="")
                return tmp;

            tmp = right_child->get_id_for_name(query);
            if(tmp!="")
                return tmp;

            return "";
        }

        return "";
    }

    /**************************************/

    void add_name_comment( string comment ) { name_comment = comment; }

    string get_name_comment() { return name_comment; }

    bool is_leaf() { return leaf; }
    void is_leaf(bool i) { leaf = i; }

    bool has_left_child() { return node_has_left_child; }
    void has_left_child(bool h) { node_has_left_child = h; }

    bool has_right_child() { return node_has_right_child; }
    void has_right_child(bool h) { node_has_right_child = h; }

    bool left_needs_correcting_sequence_site_index() { return adjust_left_node_site_index; }
    void left_needs_correcting_sequence_site_index(bool h) { adjust_left_node_site_index = h; }

    bool right_needs_correcting_sequence_site_index() { return adjust_right_node_site_index; }
    void right_needs_correcting_sequence_site_index(bool h) { adjust_right_node_site_index = h; }

    /**************************************/

    void add_left_child(Node *child)
    {
        left_child = child; leaf = false; this->has_left_child(true);
    }

    void add_right_child(Node *child)
    {
        right_child = child; leaf = false; this->has_right_child(true);
    }
    Node * get_left_child() { return left_child; }

    Node * get_right_child() { return right_child; }

    /**************************************/


    void get_all_nodes(vector<Node*> *nodes,bool visit_leaves_first=false)
    {
        if(!leaf)
            left_child->get_all_nodes(nodes,visit_leaves_first);

        if(not visit_leaves_first)
            nodes->push_back(this);

        if(!leaf)
            right_child->get_all_nodes(nodes,visit_leaves_first);

        if(visit_leaves_first)
            nodes->push_back(this);
    }

    void get_leaf_nodes(vector<Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_leaf_nodes(nodes);
            right_child->get_leaf_nodes(nodes);
        }
        else
            nodes->push_back(this);
    }

    void get_leaf_node_names(set<string> *names)
    {
        if(!this->is_leaf())
        {
            left_child->get_leaf_node_names(names);
            right_child->get_leaf_node_names(names);
        }
        else
            names->insert(this->get_name());
    }

    void get_leaf_nodes(map<string,Node*> *nodes)
    {
        if(this->is_leaf())
        {
            nodes->insert(pair<string,Node*>(this->get_name(),this));
        }
        else
        {
            left_child->get_leaf_nodes(nodes);
            right_child->get_leaf_nodes(nodes);
        }
    }

    void get_internal_nodes(vector<Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_nodes(nodes);
            nodes->push_back(this);
            right_child->get_internal_nodes(nodes);
        }
    }

    void get_read_nodes_below(vector<Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_read_nodes_below(nodes);
            right_child->get_read_nodes_below(nodes);
        }
        else if(this->get_sequence()->is_read_sequence() )
        {
            nodes->push_back(this);
        }
        return;
    }

    void get_all_nodes(map<string,Node*> *nodes)
    {
        if(!this->is_leaf())
            left_child->get_all_nodes(nodes);

        nodes->insert(pair<string,Node*>(this->get_name(),this));

        if(!this->is_leaf())
            right_child->get_all_nodes(nodes);
    }

    void get_internal_nodes(map<string,Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_nodes(nodes);
            nodes->insert(pair<string,Node*>(this->get_name(),this));
            right_child->get_internal_nodes(nodes);
        }
    }

    void get_internal_node_names(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_node_names(list);
            right_child->get_internal_node_names(list);

            list->insert(pair<string,string>(this->get_name(),this->get_name()));
//            list->insert(pair<string,string>(this->get_nhx_tid(),this->get_name()));
        }
    }

    void get_internal_node_names(set<string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_node_names(list);
            right_child->get_internal_node_names(list);

            list->insert(this->get_name());
        }
    }


    void get_node_names(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_node_names(list);
            right_child->get_node_names(list);
        }
        list->insert(pair<string,string>(this->get_name(),this->get_name()));
    }

    void get_node_names(set<string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_node_names(list);
            right_child->get_node_names(list);
        }
        list->insert(this->get_name());
    }

    void get_terminal_node_names(multimap<string,string> *list)
    {
        if(this->is_leaf())
            list->insert(pair<string,string>(this->get_name(),this->get_name()));
        else
        {
            left_child->get_terminal_node_names(list);
            right_child->get_terminal_node_names(list);
        }
    }

    void get_internal_node_names_with_tid_tag(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_node_names_with_tid_tag(list);
            right_child->get_internal_node_names_with_tid_tag(list);

            if(this->get_nhx_tid()!="")
            {
                list->insert(pair<string,string>(this->get_nhx_tid(),this->get_name()));
            }
        }
    }

    void get_node_names_with_tid_tag(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_node_names_with_tid_tag(list);
            right_child->get_node_names_with_tid_tag(list);
        }
        if(this->get_nhx_tid()!="")
        {
            list->insert(pair<string,string>(this->get_nhx_tid(),this->get_name()));
        }
    }

    void get_node_names_with_tid_tag(set<string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_node_names_with_tid_tag(list);
            right_child->get_node_names_with_tid_tag(list);
        }
        if(this->get_nhx_tid()!="")
        {
            list->insert(this->get_name());
        }
    }
    void get_node_names_with_tid_exact(set<string> *list, string *tid)
    {
        if(!this->is_leaf())
        {
            left_child->get_node_names_with_tid_exact(list, tid);
            right_child->get_node_names_with_tid_exact(list, tid);
        }
        if(tid != 0 && this->get_nhx_tid() == *tid)
        {
            list->insert(this->get_name());
        }
    }

    void name_internal_nodes(int *count)
    {
        if(leaf)
            return;
        else
        {
            left_child->name_internal_nodes(count);
            right_child->name_internal_nodes(count);

            stringstream ss;
            ss<<"#"<<(*count)<<"#";

            this->set_name(ss.str());
            (*count)++;
        }

    }

    void get_all_terminal_node_names(set<string> *node_names)
    {
        if(this->is_leaf())
        {
            node_names->insert(this->get_name());
        }
        else
        {
            left_child->get_all_terminal_node_names(node_names);
            right_child->get_all_terminal_node_names(node_names);
        }
    }

    void get_read_node_names(set<string> *read_names)
    {
        if(this->is_leaf())
        {
            if(this->get_sequence()->is_read_sequence())
            {
                read_names->insert(this->get_name());
            }
        }
        else
        {
            left_child->get_read_node_names(read_names);
            right_child->get_read_node_names(read_names);
        }
    }

    void get_read_dna_sequences(map<string,string> *dna_seqs)
    {
        if(this->is_leaf())
        {
            if(this->get_sequence()->is_read_sequence())
            {
                dna_seqs->insert(pair<string,string>(this->get_name(),this->get_Orf().dna_sequence));
            }
        }
        else
        {
            left_child->get_read_dna_sequences(dna_seqs);
            right_child->get_read_dna_sequences(dna_seqs);
        }
    }

    void set_node_names_for_exonerate(set<string> *names)
    {
        if(!leaf)
        {
            this->get_left_child()->set_node_names_for_exonerate(names);
            this->get_right_child()->set_node_names_for_exonerate(names);
        }

        if(names->find(this->get_name()) != names->end() )
            use_for_exonarate = true;
        else
            use_for_exonarate = false;
    }

    void get_node_names_for_exonerate(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            this->get_left_child()->get_node_names_for_exonerate(list);
            this->get_right_child()->get_node_names_for_exonerate(list);
        }

        if(use_for_exonarate)
            list->insert(pair<string,string>(this->get_name(),this->get_name()));
    }

    void set_use_for_exonerate(bool v) { use_for_exonarate = v; }
    bool get_use_for_exonerate() { return use_for_exonarate; }

    bool sequence_site_index_needs_correcting()
    {
        if(this->is_leaf())
            return false;
        else if(adjust_left_node_site_index || adjust_right_node_site_index)
            return true;
        else if(left_child->sequence_site_index_needs_correcting())
            return true;
        else if(right_child->sequence_site_index_needs_correcting())
            return true;
        else
            return false;
    }

    /************************************/

    void get_state_identity(int pos,int query_state,vector<int> *scores,bool visit_leaves_first=false,bool visit_leaves_only=false)
    {
        if(leaf)
        {
            if(this->get_sequence()->get_site_at(pos)->get_state() == query_state)
                scores->push_back(1);
            else
                scores->push_back(0);
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(pos)->get_children();

            int lj = offspring->left_index;
            if(lj>=0)
            {
                left_child->get_state_identity(lj,query_state,scores,visit_leaves_first);
            }
            else
            {
                for(int i=0;i<left_child->get_number_of_nodes();i++)
                     scores->push_back(0);
            }

            if(not visit_leaves_first && not visit_leaves_only)
            {
                if(this->get_sequence()->get_site_at(pos)->get_state() == query_state)
                    scores->push_back(1);
                else
                    scores->push_back(0);
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                right_child->get_state_identity(rj,query_state,scores,visit_leaves_first);
            }
            else
            {
                for(int i=0;i<right_child->get_number_of_nodes();i++)
                     scores->push_back(0);
            }

            if(visit_leaves_first && not visit_leaves_only)
            {
                if(this->get_sequence()->get_site_at(pos)->get_state() == query_state)
                    scores->push_back(1);
                else
                    scores->push_back(0);
            }
        }
    }


    void get_subst_distance(int pos,int query_state,vector<float> *scores,Evol_model *model,bool visit_leaves_first=false,bool visit_leaves_only=false)
    {
        if(leaf)
        {
            scores->push_back( model->score(this->get_sequence()->get_site_at(pos)->get_state(), query_state) );
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(pos)->get_children();

            int lj = offspring->left_index;
            if(lj>=0)
            {
                left_child->get_subst_distance(lj,query_state,scores,model,visit_leaves_first);
            }
            else
            {
                for(int i=0;i<left_child->get_number_of_nodes();i++)
                     scores->push_back(0);
            }

            if(not visit_leaves_first && not visit_leaves_only)
            {
                scores->push_back( model->score(this->get_sequence()->get_site_at(pos)->get_state(), query_state) );
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                right_child->get_subst_distance(rj,query_state,scores,model,visit_leaves_first);
            }
            else
            {
                for(int i=0;i<right_child->get_number_of_nodes();i++)
                     scores->push_back(0);
            }

            if(visit_leaves_first && not visit_leaves_only)
            {
                scores->push_back( model->score(this->get_sequence()->get_site_at(pos)->get_state(), query_state) );
            }
        }
    }

    Node * get_parent_node(string node_name)
    {
        if(leaf)
            return 0;

        if(this->get_left_child()->get_name()==node_name)
            return this;
        if(this->get_right_child()->get_name()==node_name)
            return this;

        Node *tmp;

        tmp = this->get_left_child()->get_parent_node(node_name);
        if(tmp!=0)
            return tmp;
        tmp = this->get_right_child()->get_parent_node(node_name);
        if(tmp!=0)
            return tmp;

        return 0;
    }


    /************************************/

    void has_sequence(bool s) { node_has_sequence = s; }
    bool has_sequence() { return node_has_sequence; }

    void prune_tree() { this->prune_down(); this->prune_up(); }
    void prune_up();
    void prune_down();

    void unset_has_sequence()
    {
        if(this->is_leaf())
        {
            this->has_sequence(false);
        }
        else
        {
            this->get_left_child()->unset_has_sequence();
            this->get_right_child()->unset_has_sequence();
        }
    }

    void unset_has_sequence(set<string> *remove)
    {
        if(this->is_leaf())
        {
            if(remove->find(this->get_name())!=remove->end())
            {
                this->has_sequence(false);
            }
        }
        else
        {
            this->get_left_child()->unset_has_sequence(remove);
            this->get_right_child()->unset_has_sequence(remove);
        }
    }

    void set_has_sequence()
    {
        if(this->is_leaf())
        {
            this->has_sequence(true);
        }
        else
        {
            this->get_left_child()->set_has_sequence();
            this->get_right_child()->set_has_sequence();
        }
    }

    void set_has_sequence(set<string> *keep)
    {
        if(this->is_leaf())
        {
            if(keep->find(this->get_name())!=keep->end())
            {
                this->has_sequence(true);
            }
        }
        else
        {
            this->get_left_child()->set_has_sequence(keep);
            this->get_right_child()->set_has_sequence(keep);
        }
    }


    void set_reads_has_sequence()
    {
        if(this->is_leaf() )
        {
            if(this->get_sequence()->is_read_sequence())
            {
                this->has_sequence(true);
            }
        }
        else
        {
            this->get_left_child()->set_reads_has_sequence();
            this->get_right_child()->set_reads_has_sequence();
        }
    }

    bool all_descendants_reads()
    {
        if(this->is_leaf())
        {
            if(this->get_sequence()->is_read_sequence())
                return true;
            else
                return false;
        }
        else
        {
            return this->get_left_child()->all_descendants_reads()
                && this->get_right_child()->all_descendants_reads();
        }
    }

    void get_closest_reference_leaf(float *dist,string *name)
    {
        if(this->is_leaf())
        {
            if(this->get_sequence()->is_read_sequence())
            {
                *dist = -1;
                *name = "";
            }
            else
            {
                *dist = this->get_distance_to_parent();
                *name = this->get_name();
            }
        }
        else
        {
            float ld = -1;
            string ln = "";
            this->get_left_child()->get_closest_reference_leaf(&ld,&ln);
            float rd = -1;
            string rn = "";
            this->get_right_child()->get_closest_reference_leaf(&rd,&rn);

            if(ld>=0 && ld<rd)
            {
                *dist = ld+this->get_distance_to_parent();
                *name = ln;
            }
            else if(rd>=0)
            {
                *dist = rd+this->get_distance_to_parent();
                *name = rn;
            }
        }
    }

    void get_closest_reference_leaves(set<string> *names)
    {
        if(!this->is_leaf())
        {
            this->get_left_child()->get_closest_reference_leaves(names);
            this->get_right_child()->get_closest_reference_leaves(names);

            bool lr = this->get_left_child()->all_descendants_reads();
            bool rr = this->get_right_child()->all_descendants_reads();

            if(not lr && rr)
            {
                float dist = -1;
                string name = "";
                this->get_left_child()->get_closest_reference_leaf(&dist,&name);
                names->insert(name);
            }
            else if(lr && not rr)
            {
                float dist = -1;
                string name = "";
                this->get_right_child()->get_closest_reference_leaf(&dist,&name);
                names->insert(name);
            }
        }
    }


    /************************************/

    /*
     * Note on graphical output: posterior plots have to written *during* the alignment
     * as the dynamic programming matrices aren't stored for later use.
     * In contrast, sequence graphs exist till the very end.
     */
    void start_alignment(Model_factory *mf)
    {
        
       // Log_output::write_header("Performing multiple alignment ",0);

        if(Settings_handle::st.is("mpost-posterior-plot-file"))
        {
            this->start_mpost_plot_file();
        }

        this->number_of_nodes = this->get_number_of_leaves()-1;
        this->alignment_number = 1;

        this->align_sequences(mf);

        if(Settings_handle::st.is("mpost-posterior-plot-file"))
        {
            this->finish_mpost_plot_file();
        }

        if(Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors"))
        {
//                show_seqs();
            this->reconstruct_parsimony_ancestor(mf);
        }
      
    }


    void read_reference_alignment(Model_factory *mf, bool quick=false)
    {
        //Log_output::write_header("Reading reference alignment ",0);

        if(quick || Settings_handle::st.is("quick"))
        {
            Evol_model model = mf->alignment_model(.5);
            this->read_alignment(&model);
        }
        else
        {
            this->read_alignment(mf);
        }

        this->reconstruct_parsimony_ancestor(mf);
    }



    void align_sequences(Model_factory *mf)
    {
        if(leaf)
            return;
             
        left_child->align_sequences(mf);
        right_child->align_sequences(mf);

        this->align_sequences_this_node(mf);
        
    }

    /*******************************************************************************/

    void align_sequences_this_node(Model_factory *mf, bool is_reads_sequence=false);

    /*******************************************************************************/

    void start_threaded_alignment(Model_factory *mf, int n_threads);
    
    /*******************************************************************************/

    void start_openmp_alignment(Model_factory *mf, int n_threads);

    /*******************************************************************************/

    void build_queues(vector<Node*>& wait_nodes, vector<Node*>& run_nodes);

    /*******************************************************************************/

    void threaded_function(Model_factory *mf, vector<Node*>& wait_nodes, vector<Node*>& run_nodes);

    /*******************************************************************************/

    void align_sequences_this_node_threaded(Model_factory *mf);

    /*******************************************************************************/

    void align_sequences_this_node_openmp(Model_factory *mf);

    /*******************************************************************************/

    void print_alignment()
    {
        vector<Fasta_entry> alignment;
        this->get_alignment(&alignment);
        for(int i=0;i<(int)alignment.size();i++)
            Log_output::write_out(">"+alignment.at(i).name+"\n"+alignment.at(i).sequence+"\n",1);
    }

    void read_alignment(Model_factory *mf)
    {
        if(leaf)
            return;

        left_child->read_alignment(mf);
        right_child->read_alignment(mf);

        this->read_alignment_this_node(mf);
    }

    void read_alignment_this_node(Model_factory *mf)
    {

        Log_output::write_msg("reading alignment node "+this->get_name()+": "+left_child->get_name()+" - "+right_child->get_name()+".",0);

        double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();
        Evol_model model = mf->alignment_model(dist);

        Reference_alignment ra;
        ra.read_alignment(left_child->get_sequence(),right_child->get_sequence(),&model,
                 left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

        this->add_ancestral_sequence( ra.get_simple_sequence() );
    }

    void read_alignment(Evol_model *model)
    {
        if(leaf)
            return;

        left_child->read_alignment(model);
        right_child->read_alignment(model);

        this->read_alignment_this_node(model);
    }

    void read_alignment_this_node(Evol_model *model)
    {

        Log_output::write_msg("reading alignment node "+this->get_name()+": "+left_child->get_name()+" - "+right_child->get_name()+".",0);

        Reference_alignment ra;
        ra.read_alignment(left_child->get_sequence(),right_child->get_sequence(),model,
                 left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

        this->add_ancestral_sequence( ra.get_simple_sequence() );
    }

    void reconstruct_parsimony_ancestor(Model_factory *mf)
    {
        for(int i=1;i<this->get_sequence()->sites_length()-1;i++)
        {
            int state = this->get_sequence()->get_site_at(i)->get_state();
            this->reconstruct_parsimony_ancestor_at_site(mf,i,state,false);
        }
    }

    void reconstruct_parsimony_ancestor_at_site(Model_factory *mf,int pos,int parent_state,bool is_matched)
    {

        if(leaf)
            return;

        Site *site = this->get_sequence()->get_site_at(pos);

        int pstate = site->get_path_state();

        if( pstate == Site::matched )
        {
            int new_state = mf->get_child_parsimony_state(parent_state,site->get_state());

            site->set_state(new_state);

            is_matched = true;
        }
        else
        {
            site->set_state(parent_state);
        }

        if( !is_matched )
        {
            site->set_site_type(Site::non_real);
        }

        if(site->children.left_index >= 0)
        {
            this->left_child->reconstruct_parsimony_ancestor_at_site(mf,site->children.left_index,site->get_state(),is_matched);
        }
        if(site->children.right_index >= 0)
        {
            this->right_child->reconstruct_parsimony_ancestor_at_site(mf,site->children.right_index,site->get_state(),is_matched);
        }

    }

    void reconstruct_one_parsimony_ancestor(Model_factory *mf,bool do_left)
    {
        for(int i=1;i<this->get_sequence()->sites_length()-1;i++)
        {
            Site *site = this->get_sequence()->get_site_at(i);
            int state = site->get_state();
            if(do_left)
            {
                if(site->children.left_index >= 0)
                    this->get_left_child()->reconstruct_one_parsimony_ancestor_at_site(mf,site->children.left_index,state);
            }
            else
            {
                if(site->children.right_index >= 0)
                    this->get_right_child()->reconstruct_one_parsimony_ancestor_at_site(mf,site->children.right_index,state);
            }
        }
    }

    void reconstruct_one_parsimony_ancestor_at_site(Model_factory *mf,int pos,int parent_state)
    {

        if(leaf)
            return;

        Site *site = this->get_sequence()->get_site_at(pos);

        int new_state = mf->get_child_parsimony_state(parent_state,site->get_state());
        site->set_state(new_state);

    }

    bool has_site_at_alignment_column(int j,string node_name)
    {

        if(this->get_name() == node_name)
        {
            return true;
        }
        else if(leaf)
        {
            return false;
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(j)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
            {
                bool l = left_child->has_site_at_alignment_column(lj,node_name);
                if(l)
                    return true;
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                bool r = right_child->has_site_at_alignment_column(rj,node_name);
                if(r)
                    return true;
            }
        }
        return false;
    }

    bool any_other_has_site_at_alignment_column(int j,string node_name)
    {
        if(leaf && this->get_name() != node_name)
        {
            return true;
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(j)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
            {
                bool l = left_child->any_other_has_site_at_alignment_column(lj,node_name);
                if(l)
                    return true;
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                bool r = right_child->any_other_has_site_at_alignment_column(rj,node_name);
                if(r)
                    return true;
            }
        }
        return false;
    }


    int get_state_at_alignment_column(int j,string node_name)
    {

        if(this->get_name() == node_name)
        {
            return this->get_sequence()->get_site_at(j)->get_state();
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(j)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
            {
                int l = left_child->get_state_at_alignment_column(lj,node_name);
                if(l>=0)
                    return l;
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                int r = right_child->get_state_at_alignment_column(rj,node_name);
                if(r>=0)
                    return r;
            }
        }
        return -1;
    }

    void find_closest_reference(string *ref_name,string query_name)
    {

        if(leaf)
            *ref_name = this->get_name();
        else
        {
            string r_name;
            if(this->get_left_child()->has_reference_child(&r_name) &&
                    this->get_right_child()->has_child_named(query_name))
            {
                *ref_name = r_name;
            }
            if(!leaf)
            {
                this->get_left_child()->find_closest_reference(ref_name,query_name);
                this->get_right_child()->find_closest_reference(ref_name,query_name);
            }
        }
    }

    bool has_reference_child(string *ref_name)
    {
        if(leaf && !this->get_sequence()->is_read_sequence())
        {
            *ref_name = this->get_name();
            return true;
        }
        if(!leaf)
        {
            bool has = this->get_left_child()->has_reference_child(ref_name);
            if(has)
                return true;
            return this->get_right_child()->has_reference_child(ref_name);
        }
    }

    bool has_child_named(string search_name)
    {
        if(this->get_name()==name)
            return true;
        if(!leaf)
        {
            bool has = this->get_left_child()->has_child_named(search_name);
            if(has)
                return true;
            return this->get_right_child()->has_child_named(search_name);
        }
    }

    bool has_additional_sites_before_alignment_column(int j)
    {

        Site_children *offspring = sequence->get_site_at(j)->get_children();

        int lj = offspring->left_index;
        if(lj>=0)
        {
            bool l = left_child->has_additional_sites_before_alignment_column(lj);
            if(l)
                return true;
        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            bool r = right_child->has_additional_sites_before_alignment_column(rj);
            if(r)
                return true;
        }

        if(!(this->adjust_left_node_site_index && this->adjust_right_node_site_index))
        {
            return false;
        }
        else
        {
            int prev_lj = -1;
            int prev_rj = -1;

            if(j>0)
            {
                Site_children *prev_offspring = sequence->get_site_at(j-1)->get_children();
                prev_lj = prev_offspring->left_index;
                prev_rj = prev_offspring->right_index;
            }

            if(this->adjust_left_node_site_index && lj-prev_lj!=1)
                return true;
            else if(this->adjust_right_node_site_index && rj-prev_rj!=1)
                return true;
            else
                return false;
        }
    }

    void additional_sites_before_alignment_column(int j,vector<Insertion_at_node> *addition);

    void start_mpost_plot_file()
    {
        string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), (ios::out) );
        if (! output) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        ofstream output2(path2.c_str(), (ios::out) );
        if (! output2) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }


        output<<"vardef circle(expr a,c,col) =\n path p;\n h := 0.6mm;\n"
                " p = a+(-h/2,0){up}..{right}a+(0,h/2){right}..{down}a+(h/2,0){down}..{left}a+(0,-h/2){left}..{up}cycle;\n"
                " fill p withcolor col;\n draw(p);\n pair x;\n x = ((point(0) of p)+(point(2) of p))/2;\n"
                " label(c,x);\n p\nenddef;\n\n";

        output<<"vardef edge(expr px,ax,py,ay,s) =\n pair x,y;\n x = ((point(0) of px)+(point(2) of px))/2;\n"
                " y = ((point(0) of py)+(point(2) of py))/2;\n path p; p = (x){dir(ax)}..{dir(ay)}(y);\n"
                " draw (p) cutbefore px cutafter py withpen pencircle scaled (0.3*s);\n point .5*length p of p\n"
                "enddef;\n\n";

        output<<"def edgetop(expr px,py,a,c,s) =\n label.top(c,edge(px,a,py,-1*a,s));\nenddef;\n"
                "def edgebot(expr px,py,a,c,s) =\n label.bot(c,edge(px,360-a,py,360+1*a,s));\nenddef;\n\n";
        output<<"def edgelft(expr px,py,a,c,s) =\n label.lft(c,edge(px,a,py,180-1*a,s));\nenddef;\n"
                "def edgergt(expr px,py,a,c,s) =\n label.rt(c,edge(px,180-a,py,a,s));\nenddef;\n\n";

        output2<<"\\documentclass{article}\n\\usepackage{geometry,graphicx}\n"
                 "\\geometry{a4paper,tmargin=0.5in,bmargin=0.5in,lmargin=.5in,rmargin=0.5in}\n"
                 "\\begin{document}\n";


        output.close();
        output2.close();
    }

    void finish_mpost_plot_file()
    {
        string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();
        Log_output::write_out("Plot file: "+file+"\n",1);

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), (ios::out|ios::app) );
        if (! output) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        ofstream output2(path2.c_str(), (ios::out|ios::app) );
        if (! output2) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        output<<"end;\n";
        output2<<"\\end{document}\n";

        output.close();
        output2.close();

        Log_output::write_out("\nThe posterior probability plot files can generated using following commands:\n"
              "  mpost "+file+".mp\n"
              "  latex "+file+".tex\n"
              "  dvipdf "+file+".dvi\n\n",1);

    }

    /************************************/

    void show_seqs()
    {
        stringstream ss;
        ss<<"node "<<get_name();
        for(int i=0;i<sequence->sites_length();i++)
            ss<<"; ["<<i<<"] "<<sequence->get_site_at(i)->get_state();
        ss<<endl<<endl;
        Log_output::write_out(ss.str(),1);

        if(!is_leaf())
            get_left_child()->show_seqs();

        if(!is_leaf())
            get_right_child()->show_seqs();
    }

    void get_node_sequence(Fasta_entry *seq);
    void get_alignment(vector<Fasta_entry> *aligned_sequences, bool include_internal_nodes=false);
    void get_alignment_for_nodes(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes);
    void add_root_consensus(vector<Fasta_entry> *aligned_sequences);

    void get_dna_sequences(map<string,string*> *dna_sequences)
    {
        if(leaf)
        {
            dna_sequences->insert(make_pair( this->get_name(),this->get_sequence()->get_dna_sequence() ) );
        }
        else
        {
            this->get_left_child()->get_dna_sequences(dna_sequences);
            this->get_right_child()->get_dna_sequences(dna_sequences);
        }
    }

    void get_dna_sequences(map<string,string> *dna_sequences)
    {
        if(leaf)
        {
            dna_sequences->insert(make_pair( this->get_name(),*this->get_sequence()->get_dna_sequence() ) );
        }
        else
        {
            this->get_left_child()->get_dna_sequences(dna_sequences);
            this->get_right_child()->get_dna_sequences(dna_sequences);
        }
    }

    void get_unaligned_sequences(map<string,string*> *unaligned_sequences)
    {
        if(leaf)
        {
            unaligned_sequences->insert(make_pair( this->get_name(),this->get_sequence()->get_unaligned_sequence() ) );
        }
        else
        {
            this->get_left_child()->get_unaligned_sequences(unaligned_sequences);
            this->get_right_child()->get_unaligned_sequences(unaligned_sequences);
        }
    }

    void get_unaligned_sequences(map<string,string> *unaligned_sequences)
    {
        if(leaf)
        {
            unaligned_sequences->insert(make_pair( this->get_name(),*this->get_sequence()->get_unaligned_sequence() ) );
        }
        else
        {
            this->get_left_child()->get_unaligned_sequences(unaligned_sequences);
            this->get_right_child()->get_unaligned_sequences(unaligned_sequences);
        }
    }

    string *get_dna_sequence_for_node(string node_name)
    {
        if(is_leaf() && this->get_name() == node_name)
        {
            return this->get_sequence()->get_dna_sequence();
        }
        else if(is_leaf())
        {
            return 0;
        }
        else
        {
            string *rv = this->get_left_child()->get_dna_sequence_for_node(node_name);
            if(rv != 0)
                return rv;
            rv = this->get_right_child()->get_dna_sequence_for_node(node_name);
            if(rv != 0)
                return rv;
        }
        return 0;
    }

    void get_alignment_column_at(int j,vector<string> *column, bool include_internal_nodes);

    void get_multiple_alignment_columns_before(Insertion_at_node ins,vector<string> *columns,bool include_internal_nodes);

    void get_multiple_alignment_columns_before(int j,vector< vector<string> > *columns, string node_name_wanted,
                                                bool left_child_wanted, bool include_internal_nodes);

    int get_number_of_leaves()
    {
        if(leaf)
             return 1;
        else
            return left_child->get_number_of_leaves()+right_child->get_number_of_leaves();
    }

    int get_weighted_number_of_leaves()
    {
        if(leaf)
             return this->get_sequence()->get_num_duplicates();
        else
            return left_child->get_weighted_number_of_leaves()+right_child->get_weighted_number_of_leaves();
    }

    int get_number_of_read_leaves_for_read_node()
    {
        if(!this->get_sequence()->is_read_sequence())
            return 0;

        if(leaf)
             return 1;
        else
            return left_child->get_number_of_read_leaves()+right_child->get_number_of_read_leaves();
    }

    int get_number_of_read_leaves()
    {
        if(leaf)
             return 1;
        else
            return left_child->get_number_of_read_leaves()+right_child->get_number_of_read_leaves();
    }

    int get_number_of_query_leaves()
    {
        if(leaf)
        {
            if(this->get_sequence()->is_read_sequence())
                return 1;
            else
                return 0;
        }
        else
            return left_child->get_number_of_query_leaves()+right_child->get_number_of_query_leaves();
    }

    int get_number_of_nodes()
    {
        if(leaf)
             return 1;
        else
            return 1+left_child->get_number_of_nodes()+right_child->get_number_of_nodes();
    }

    void set_number_of_nodes(int i) { number_of_nodes = i; }

    /************************************/

    int get_number_of_gaps_at_site(int i)
    {
        int ng = 0;

        if(!leaf)
        {
            Site_children *offspring = sequence->get_site_at(i)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
                ng += left_child->get_number_of_gaps_at_site(lj);
            else
                ng += left_child->get_number_of_leaves();

            int rj = offspring->right_index;
            if(rj>=0)
                ng += right_child->get_number_of_gaps_at_site(rj);
            else
                ng += right_child->get_number_of_leaves();
        }

        return ng;
    }


    /************************************/

    void get_parent_child_pairs(vector<pair<string,string> > *pairs)
    {
        if(!leaf)
        {
            left_child->get_parent_child_pairs(pairs);
            right_child->get_parent_child_pairs(pairs);
            pairs->push_back(pair<string,string>(this->get_name(),left_child->get_name()));
            pairs->push_back(pair<string,string>(this->get_name(),right_child->get_name()));
        }
        else
            return;
    }

    string print_tree(bool int_names=false) {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree(int_names)<<","<<right_child->print_subtree(int_names)<<")";
            if(int_names)
                ss<<this->get_name()<<":0";
            ss<<";";
            return ss.str();
        } else {
            return "";
        }
    }


    string print_subtree(bool int_names=false) {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree(int_names)<<","<<right_child->print_subtree(int_names)<<")";
            if(int_names)
                ss<<this->get_name();
            ss<<":"<<dist_to_parent;
            return ss.str();
        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent;
            return ss.str();
        }
    }

    /************************************/

    void write_nhx_tree(string path, string suffix, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+"."+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app));

        if (! output) { throw IOException ("Node::write_nhx_tree. Failed to open file"); }

        output<<this->print_nhx_tree();
        output.close();
    }

    string print_nhx_tree() const {
        if(!leaf)
        {
            string tid = this->get_nhx_tag();
            if(this->get_nhx_tid() != "")
                tid += ":TID="+this->get_nhx_tid();

            if(tid != "")
                tid = "["+tid+"]";

            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree()<<","<<right_child->print_nhx_subtree()<<"):"<<dist_to_parent<<tid<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    string print_nhx_subtree() const {

        string tid = this->get_nhx_tag();
        if(this->get_nhx_tid() != "")
            tid += ":TID="+this->get_nhx_tid();

        if(tid != "")
            tid = "["+tid+"]";

        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree()<<","<<right_child->print_nhx_subtree()<<"):"<<dist_to_parent<<tid;
            return ss.str();

        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent<<tid;
            return ss.str();
        }
    }

    /************************************/

    string print_nhx_tree_with_intIDs() const {
        if(!leaf)
        {
            string tid = this->get_nhx_tag();
            if(this->get_nhx_tid() != "")
                tid += ":TID="+this->get_nhx_tid();

            if(tid != "")
                tid = "["+tid+"]";

            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree_with_intIDs()<<","<<right_child->print_nhx_subtree_with_intIDs()<<")"<<this->get_name()<<":"<<dist_to_parent<<tid<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    string print_nhx_subtree_with_intIDs() const {

        string tid = this->get_nhx_tag();
        if(this->get_nhx_tid() != "")
            tid += ":TID="+this->get_nhx_tid();

        if(tid != "")
            tid = "["+tid+"]";

        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree_with_intIDs()<<","<<right_child->print_nhx_subtree_with_intIDs()<<")"<<this->get_name()<<":"<<dist_to_parent<<tid;
            return ss.str();

        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent<<tid;
            return ss.str();
        }
    }

    /************************************/

    string print_bppa_tree(int *add) const {
        if(!leaf)
        {
            stringstream ss(name);
            char b,e; int n;
            ss >> b >> n >> e;
            ss.clear();
            ss.str("");
            ss<<"("<<left_child->print_bppa_subtree(add)<<","<<right_child->print_bppa_subtree(add)<<")"<<":0[&&NHX:ND="<<n<<"];";
            return ss.str();
        } else {
            return "";
        }
    }

    string print_bppa_subtree(int *add) const {
        if(!leaf)
        {
            stringstream ss(name);
            char b,e; int n;
            ss >> b >> n >> e;
            ss.clear();
            ss.str("");
            ss<<"("<<left_child->print_bppa_subtree(add)<<","<<right_child->print_bppa_subtree(add)<<"):"<<dist_to_parent<<"[&&NHX:ND="<<n<<"]";
            return ss.str();

        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent<<"[&&NHX:ND="<<*add<<"]";
            (*add)++;
            return ss.str();
        }
    }


    /************************************/

    string print_xml_tree() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_xml_subtree()<<","<<right_child->print_xml_subtree()<<")"<<name<<":0";

            if(Settings_handle::st.is("xml-nhx"))
            {
                string tid = this->get_nhx_tag();
                if(this->get_nhx_tid() != "")
                    tid += ":TID="+this->get_nhx_tid();

                if(tid != "")
                    tid = "["+tid+"]";
                ss<<tid;
            }

            ss<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    string print_xml_subtree() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_xml_subtree()<<","<<right_child->print_xml_subtree()<<")"<<name<<":"<<dist_to_parent;

            if(Settings_handle::st.is("xml-nhx"))
            {
                string tid = this->get_nhx_tag();
                if(this->get_nhx_tid() != "")
                    tid += ":TID="+this->get_nhx_tid();

                if(tid != "")
                    tid = "["+tid+"]";
                ss<<tid;
            }

            return ss.str();
        } else {
            stringstream ss;
            ss<<name_id<<":"<<dist_to_parent;

            if(Settings_handle::st.is("xml-nhx"))
            {
                string tid = this->get_nhx_tag();
                if(this->get_nhx_tid() != "")
                    tid += ":TID="+this->get_nhx_tid();

                if(tid != "")
                    tid = "["+tid+"]";
                ss<<tid;
            }

            return ss.str();
        }
    }

    /************************************/

    string print_tagged_topology() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_tagged_subtopology()<<","<<right_child->print_tagged_subtopology()<<")"<<name<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    string print_tagged_subtopology() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_tagged_subtopology()<<","<<right_child->print_tagged_subtopology()<<")"<<name;
            return ss.str();
        } else {
            return name;
        }
    }

    /************************************/


    void get_alignment_for_reads(vector<Fasta_entry> *aligned_sequences, bool show_ref_insertions);
    void get_alignment_for_read_nodes(vector<Fasta_entry> *aligned_sequences, bool show_ref_insertions);
    void get_alignment_column_for_reads_at(int j,vector<string> *column, bool *has_characters);

    bool site_in_reference(int i)
    {
        if(!this->get_sequence()->is_read_sequence())
        {
            return true;
        }
        else
        {
            if(!is_leaf())
            {
                Site_children *offspring = sequence->get_site_at(i)->get_children();
                int li = offspring->left_index;
                if(li>=0)
                {
                    if(left_child->site_in_reference(li))
                        return true;
                }
                int ri = offspring->right_index;
                if(ri>=0)
                {
                    if(right_child->site_in_reference(ri))
                        return true;
                }
            }
        }
        return false;
    }

    bool site_in_read(int i)
    {
        if(this->get_sequence()->is_read_sequence())
        {
            return true;
        }
        else
        {
            if(!is_leaf())
            {
                Site_children *offspring = sequence->get_site_at(i)->get_children();
                int li = offspring->left_index;
                if(li>=0)
                {
                    if(left_child->site_in_read(li))
                        return true;
                }
                int ri = offspring->right_index;
                if(ri>=0)
                {
                    if(right_child->site_in_read(ri))
                        return true;
                }
            }
        }
        return false;
    }


    string find_first_nonread_left_parent()
    {
        string name = "";
        if(this->get_sequence()->is_read_sequence() && !this->is_leaf())
            name = left_child->find_first_nonread_left_parent();
        else
            name = this->get_name();

        return name;
    }

    void reconstruct_contigs(vector<Fasta_entry> *contigs,bool parent_is_read_sequence, bool consensus_only = false)
    {

        bool this_is_read_sequence = this->get_sequence()->is_read_sequence();
        bool show_ref_insertions = false;
        if(!parent_is_read_sequence && this_is_read_sequence)
        {

            if(Settings_handle::st.is("inlude-parent-in-contig"))
            {
                if(this->get_number_of_leaves() == this->get_number_of_read_leaves_for_read_node()+1)
                {
                    int seq_length = this->get_sequence()->sites_length();
                    string parent_name = this->find_first_nonread_left_parent();

                    Fasta_entry ref;
                    ref.name = parent_name;

                    string alpha = Model_factory::get_dna_full_char_alphabet();
                    if(this->get_sequence()->get_data_type() == Model_factory::protein)
                        alpha = Model_factory::get_protein_full_char_alphabet();

                    for(int j=1;j<seq_length-1;j++)
                    {
                        int st = this->get_state_at_alignment_column(j,parent_name);
                        if(st>=0)
                            ref.sequence.append(alpha.substr(st,1));
                        else
                            ref.sequence.append("-");
                    }

                    contigs->push_back(ref);
                    show_ref_insertions = true;
                }
            }

            Sequence *seq = this->get_sequence();
            int seq_length = seq->sites_length();
            Fasta_entry entry;
            entry.name = "consensus_"+this->find_first_nonread_left_parent();
            entry.comment = this->find_first_nonread_left_parent();

            int min_num_seqs = int( this->get_weighted_number_of_leaves() * (float)Settings_handle::st.get("consensus-minimum-proportion").as<float>() );

            if(min_num_seqs < Settings_handle::st.get("consensus-minimum").as<int>())
                min_num_seqs = Settings_handle::st.get("consensus-minimum").as<int>();

            if(seq->get_data_type() == Model_factory::dna)
            {

                string alpha = Model_factory::get_dna_full_char_alphabet();

                for(int j=1;j<seq_length-1;j++)
                {
                    Site *site = seq->get_site_at(j);
                    int sA = site->get_sumA();
                    int sC = site->get_sumC();
                    int sG = site->get_sumG();
                    int sT = site->get_sumT();

                    bool included_in_reference = this->site_in_reference(j);

                    if(included_in_reference && sA+sC+sG+sT==0)
                    {
                        int state = site->get_state();
                        int path_state = site->get_path_state();

                        if(path_state != Site::xskipped && path_state != Site::yskipped)
                        {
                            if(Settings_handle::st.is("show-contig-ancestor") && state>=0 && state < (int)alpha.length())
                            {
                                char c = alpha.at(site->get_state());
                                entry.sequence.append(1,tolower(c));
                            }
                            else
                            {
                                entry.sequence.append("n");
                            }
                        }
                        else if(show_ref_insertions)
                        {
                            entry.sequence.append("-");
                        }
                    }
                    else if(!included_in_reference && sA+sC+sG+sT < min_num_seqs)
                    {
                        entry.sequence.append("-");
                    }
                    else
                    {
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
            }
            else if(seq->get_data_type() == Model_factory::protein)
            {
/*
                string alpha = Model_factory::get_protein_full_char_alphabet();

                for(int j=1;j<seq_length-1;j++)
                {
                    Site *site = seq->get_site_at(j);
                    int sAmino = site->get_sumAmino();

                    bool included_in_reference = this->site_in_reference(j);

                    if(included_in_reference && sAmino==0)
                    {
                        int state = site->get_state();
                        int path_state = site->get_path_state();

                        if(path_state != Site::xskipped && path_state != Site::yskipped)
                        {
                            if(Settings_handle::st.is("show-contig-ancestor") && state>=0 && state < (int)alpha.length())
                            {
                                char c = alpha.at(site->get_state());
                                entry.sequence.append(1,tolower(c));
                            }
                            else
                            {
                                entry.sequence.append("x");
                            }
                        }
                        else if(show_ref_insertions)
                        {
                            entry.sequence.append("-");
                        }
                    }
                    else if(!included_in_reference && sAmino < min_num_seqs)
                    {
                        entry.sequence.append("-");
                    }
                    else
                    {
                        char c = alpha.at(site->get_state());
                        entry.sequence.append(1,c);
                    }
                }
*/
            }

            contigs->push_back(entry);

            if(!consensus_only)
            {
                vector<Fasta_entry> reads;
                this->get_alignment_for_reads(&reads, show_ref_insertions);

                for(int i=0;i<(int)reads.size();i++)
                    contigs->push_back(reads.at(i));
            }

        }

        if(!left_child->is_leaf())
            left_child->reconstruct_contigs(contigs,this_is_read_sequence);
        if(!right_child->is_leaf())
            right_child->reconstruct_contigs(contigs,this_is_read_sequence);
    }

    /************************************/

    void write_sequence_graphs(bool overwrite=true) const throw (Exception)
    {

        string file = Settings_handle::st.get("mpost-graph-file").as<string>();
        Log_output::write_out("Graph file: "+file+"\n",1);

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        if (! output) { throw IOException ("Node::write_sequence_graphs. Failed to open file"); }

        output<<"vardef circle(expr a,c,col) =\n path p;\n h := 0.2cm;\n"
                " p = a+(-h/2,0){up}..{right}a+(0,h/2){right}..{down}a+(h/2,0){down}..{left}a+(0,-h/2){left}..{up}cycle;\n"
                " fill p withcolor col;\n draw(p);\n pair x;\n x = ((point(0) of p)+(point(2) of p))/2;\n"
                " label(c,x);\n p\nenddef;\n\n";

        output<<"vardef edge(expr px,ax,py,ay,s) =\n pair x,y;\n x = ((point(0) of px)+(point(2) of px))/2;\n"
                " y = ((point(0) of py)+(point(2) of py))/2;\n path p; p = (x){dir(ax)}..{dir(ay)}(y);\n"
                " draw (p) cutbefore px cutafter py withpen pencircle scaled s;\n point .5*length p of p\n"
                "enddef;\n\n";

        output<<"def edgetop(expr px,py,a,c,s) =\n label.top(c,edge(px,a,py,-1*a,s));\nenddef;"
                "def edgebot(expr px,py,a,c,s) =\n label.bot(c,edge(px,360-a,py,360+1*a,s));\nenddef;\n\n";

        ofstream output2(path2.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        if (! output2) { throw IOException ("Node::write_sequence_graphs. Failed to open file"); }

        output2<<"\\documentclass{article}\n\\usepackage{geometry,graphicx}\n"
                "\\geometry{a4paper,tmargin=0.5in,bmargin=0.5in,lmargin=.5in,rmargin=0.5in}\n"
                "\\begin{document}\n";

        int fig_count = 0;
        int root_length = sequence->get_sites()->size();

        this->write_metapost_graphs(&output,&output2,&fig_count,root_length);

        output<<"end;\n";
        output.close();

        output2<<"\\end{document}\n";
        output2.close();


        Log_output::write_out("\nThe sequence graph files can generated using following commands:\n"
              "  mpost "+file+".mp\n"
              "  latex "+file+".tex\n"
              "  dvipdf "+file+".dvi\n\n",1);

    }

    void write_metapost_graphs(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
    {
        if(leaf)
        {
            if(Settings_handle::st.is("output-leaf-graphs"))
                this->write_metapost_sequence_graph(output, output2, count, root_length);
        }
        else
        {
            left_child->write_metapost_graphs(output, output2, count, root_length);

            this->write_metapost_sequence_graph(output, output2, count, root_length);

            if(Settings_handle::st.is("output-alignment-graphs"))
                this->write_metapost_alignment_graph(output, output2, count, root_length);

            right_child->write_metapost_graphs(output, output2, count, root_length);
        }
    }

    void write_metapost_sequence_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception);

    void write_metapost_alignment_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception);

    static string get_node_fill_color(char c)
    {
        string color = "white";
        switch(c)
        {
            case 'A': color = "0.5blue"; break;
            case 'C': color = "0.5green"; break;
            case 'G': color = "(0.85,0.85,0)"; break;
            case 'T': color = "0.5red"; break;
            case 'R': color = "(0,0.5,0.5)"; break;
            case 'Y': color = "(0,0.5,0.5)"; break;
            case 'M': color = "(0,0.5,0.5)"; break;
            case 'K': color = "(0,0.5,0.5)"; break;
            case 'W': color = "(0,0.5,0.5)"; break;
            case 'S': color = "(0,0.5,0.5)"; break;
            case 'B': color = "(0.5,0,0.5)"; break;
            case 'D': color = "(0.5,0,0.5)"; break;
            case 'H': color = "(0.5,0,0.5)"; break;
            case 'V': color = "(0.5,0,0.5)"; break;
        }
        return color;
    }

    void add_sequence( Fasta_entry seq_entry, int data_type, bool gapped = false, bool no_trimming = false, bool turn_revcomp = false);

    void add_ancestral_sequence( Sequence* s ) { sequence = s;  node_has_sequence_object = true;}

    Sequence *get_sequence() { return sequence; }

    std::string get_sequence_string(bool with_gaps=false);

    void fix_ambiguous_states(Evol_model *model);
    void get_ambiguous_states(int pos, set<int>* states);
    bool set_ambiguous_state(int pos, int state);

    void check_valid_graph() const;

};

}
#endif // NODE_H

