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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "utils/fasta_entry.h"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

namespace ppa{

struct Edge
{
private:
    int index;

    int start_site_index;
    int end_site_index;

    float posterior_weight;
    float log_posterior_weight;

    int next_fwd_edge_index;
    int next_bwd_edge_index;

    bool used_in_alignment;

    int branch_count_since_last_used;
    float branch_distance_since_last_used;
    int branch_count_as_skipped_edge;


public:
    Edge(int s, int e): index(-1), start_site_index(s), end_site_index(e), posterior_weight(1.0), log_posterior_weight(0),
                            next_fwd_edge_index(-1), next_bwd_edge_index(-1), used_in_alignment(false),
                            branch_count_since_last_used(0), branch_distance_since_last_used(0),
                            branch_count_as_skipped_edge(0){}

    Edge(int s, int e, float w): index(-1), start_site_index(s), end_site_index(e),
                            posterior_weight(w), log_posterior_weight(log(w)),
                            next_fwd_edge_index(-1), next_bwd_edge_index(-1), used_in_alignment(false),
                            branch_count_since_last_used(0), branch_distance_since_last_used(0),
                            branch_count_as_skipped_edge(0){}

    bool is_used() { return used_in_alignment; }
    void is_used(bool set) { used_in_alignment = set; }

    int get_branch_count_since_last_used() { return branch_count_since_last_used; }
    float get_branch_distance_since_last_used() { return branch_distance_since_last_used; }
    int get_branch_count_as_skipped_edge() { return branch_count_as_skipped_edge; }

    void set_branch_count_since_last_used(int s) { branch_count_since_last_used = s; }
    void set_branch_distance_since_last_used(float s) { branch_distance_since_last_used = s; }
    void set_branch_count_as_skipped_edge(int s) { branch_count_as_skipped_edge = s; }

    void increase_branch_count_as_skipped_edge() { ++branch_count_as_skipped_edge; }

    int get_start_site_index() { return start_site_index; }
    int get_end_site_index() { return end_site_index; }

    void set_start_site_index(int s) { start_site_index = s; }
    void set_end_site_index(int s) { end_site_index = s; }

    int get_index() { return index; }
    void set_index(int i) { index = i; }

    int get_next_fwd_edge_index() { return next_fwd_edge_index; }
    void set_next_fwd_edge_index(int i) { next_fwd_edge_index = i; }

    int get_next_bwd_edge_index() { return next_bwd_edge_index; }
    void set_next_bwd_edge_index(int i) { next_bwd_edge_index = i; }

    double get_posterior_weight() { return posterior_weight; }
    double get_log_posterior_weight() { return log_posterior_weight; }

    void set_weight(float w) { posterior_weight = w; log_posterior_weight = log(w); }
    void multiply_weight(float w) { posterior_weight *= w; log_posterior_weight = log(posterior_weight); }

    bool operator==(const Edge& b)
    {
        return start_site_index==b.start_site_index && end_site_index==b.end_site_index;
    }

    bool operator!=(const Edge& b)
    {
        return !(start_site_index==b.start_site_index && end_site_index==b.end_site_index);
    }

    bool operator<(const Edge& b)
    {
        return start_site_index<b.start_site_index && end_site_index<b.end_site_index;
    }

    bool operator>(const Edge& b)
    {
        return start_site_index>b.start_site_index && end_site_index>b.end_site_index;
    }

    friend ostream& operator<< (ostream &out, Edge& b)
    {
        out<<"index: "<<b.index<<"; start_site_index: "<<b.start_site_index<<"; end_site_index: "<<b.end_site_index;
        return out;
    }

};

/**************************************/

struct Site_children
{
    int left_index;
    int right_index;

public:

    bool operator==(const Site_children& b)
    {
        return left_index==b.left_index && right_index==b.right_index;
    }

    bool operator!=(const Site_children& b)
    {
        return !(left_index==b.left_index && right_index==b.right_index);
    }

    bool operator<(const Site_children& b)
    {
        return left_index<b.left_index && right_index<b.right_index;
    }

    bool operator>(const Site_children& b)
    {
        return left_index>b.left_index && right_index>b.right_index;
    }
};

/**************************************/

struct Unique_index
{
    int left_index;
    int right_index;
    int match_state;
    int site_index;
    enum Match_state {match,xgap,ygap};

public:
    Unique_index(): left_index(-1), right_index(-1), match_state(-1), site_index(-1) {}

    Unique_index(int l, int r, int s): left_index(l), right_index(r), match_state(s) {}
    Unique_index(int l, int r, int s, int i): left_index(l), right_index(r), match_state(s), site_index(i) {}

    bool operator==(const Unique_index& b)
    {
        return left_index==b.left_index && right_index==b.right_index && match_state == b.match_state;
    }

    bool operator!=(const Unique_index& b)
    {
        return !(left_index==b.left_index && right_index==b.right_index && match_state == b.match_state);
    }

    bool operator<(const Unique_index& b)
    {
        if(b.match_state == Unique_index::match)
            return (left_index<=b.left_index && right_index<b.right_index ) || (left_index<b.left_index && right_index<=b.right_index);
        else if(b.match_state == Unique_index::xgap)
            return left_index<b.left_index && right_index<=b.right_index;
        else if(b.match_state == Unique_index::xgap)
            return left_index<=b.left_index && right_index<b.right_index;
        return false;
    }

    bool operator>(const Unique_index& b)
    {
        if(b.match_state == Unique_index::match)
            return (left_index>=b.left_index && right_index>b.right_index ) || (left_index>b.left_index && right_index>=b.right_index);
        else if(b.match_state == Unique_index::xgap)
            return left_index>b.left_index && right_index>=b.right_index;
        else if(b.match_state == Unique_index::xgap)
            return left_index>=b.left_index && right_index>b.right_index;
        return false;
    }

    friend ostream& operator<< (ostream &out, Unique_index& b)
    {
        out<<"left_index: "<<b.left_index<<"; right_index: "<<b.right_index<<"; match_state: "<<b.match_state;
        return out;
    }
};

/**************************************/

struct Site
{
    int index;
    Site_children children;
    Unique_index unique_index;

    int character_state;
    string character_symbol;

    int site_type;
    enum Site_type {start_site,real_site,stop_site,break_start_site,break_stop_site,non_real};

    int path_state;
    enum Path_state {ends_site,terminal,matched,xgapped,ygapped,xskipped,yskipped};

    vector<Edge> *edges;

    int first_fwd_edge_index;
    int current_fwd_edge_index;

    int first_bwd_edge_index;
    int current_bwd_edge_index;

    float posterior_support;

    int branch_count_since_last_used;
    float branch_distance_since_last_used;

    int sumA, sumC, sumG, sumT, sumAmino;
    bool ambiguous = false;
//    int sumAA[20];
public:
    Site(vector<Edge> *e,int type=Site::real_site,int p_state=Site::terminal):index(-1),character_state(-1),character_symbol("0"),
            site_type(type),path_state(p_state),edges(e),first_fwd_edge_index(-1),current_fwd_edge_index(-1),
            first_bwd_edge_index(-1),current_bwd_edge_index(-1),posterior_support(1),
            branch_count_since_last_used(0),branch_distance_since_last_used(0),sumA(0),sumC(0),sumG(0),sumT(0),sumAmino(0) {}

    int get_branch_count_since_last_used() { return branch_count_since_last_used; }
    float get_branch_distance_since_last_used() { return branch_distance_since_last_used; }

    void set_branch_count_since_last_used(int s) { branch_count_since_last_used = s; }
    void set_branch_distance_since_last_used(float s) { branch_distance_since_last_used = s; }

    void set_edge_vector(vector<Edge> *e) { edges = e; }

    void set_state(int c) { character_state = c; }
    int  get_state() { return character_state; }

    void set_symbol(char c) {  character_symbol = string(1,c); }
    void set_symbol(string c) {  character_symbol = c; }
    string get_symbol() {  return character_symbol; }

    void set_path_state(int c) { path_state = c; }
    int  get_path_state() { return path_state; }

    void set_index(int i) { index = i; }
    int  get_index() { return index; }

    void set_site_type(int i) { site_type = i; }
    int get_site_type() { return site_type; }

    void set_posterior_support(float s) { posterior_support = s; }
    float get_posterior_support() { return posterior_support; }

    /**************************************/

    int get_sumA() { return sumA; }
    int get_sumC() { return sumC; }
    int get_sumG() { return sumG; }
    int get_sumT() { return sumT; }

    void add_sumA(int i) { sumA += i; }
    void add_sumC(int i) { sumC += i; }
    void add_sumG(int i) { sumG += i; }
    void add_sumT(int i) { sumT += i; }

//    int get_sumAmino() { return sumAmino; }
//    void add_sumAmino(int i) { sumAmino += i; }

//    int get_sumAA(int j) { return sumAA[j]; }
//    void set_sumAA(int j,int i) { sumAA[j] = i; }
//    void add_sumAA(int j,int i) { sumAA[j] += i; }

    bool is_ambiguous() { return ambiguous; }
    void is_ambiguous(bool a) { ambiguous = a; }
    /**************************************/

    void set_empty_children()
    {
        children.left_index = -1; children.right_index = -1;
    }

    void set_children(int li, int ri)
    {
        children.left_index = li; children.right_index = ri;
    }

    Site_children * get_children() { return &children; }

    /**************************************/

    void set_unique_index(int li, int ri, int st)
    {
        unique_index.left_index = li; unique_index.right_index = ri; unique_index.match_state = st;
    }

    void set_unique_index(int li, int ri, int st, int si)
    {
        unique_index.left_index = li; unique_index.right_index = ri;
        unique_index.match_state = st; unique_index.site_index = si;
    }

    Unique_index * get_unique_index() { return &unique_index; }

    /**************************************/
    void set_first_fwd_edge_index(int i)
    {
        first_fwd_edge_index = current_fwd_edge_index = i;
    }

    void set_first_bwd_edge_index(int i)
    {
        first_bwd_edge_index = current_bwd_edge_index = i;
    }

    /**************************************/

    void add_new_fwd_edge_index(int i)
    {
        if(first_fwd_edge_index<0)
        {
            this->set_first_fwd_edge_index(i);
            return;
        }
        int prev_fwd_edge_index = current_fwd_edge_index;
        current_fwd_edge_index = i;
        edges->at(prev_fwd_edge_index).set_next_fwd_edge_index(current_fwd_edge_index);
    }

    void add_new_bwd_edge_index(int i) {
        if(first_bwd_edge_index<0)
        {
            this->set_first_bwd_edge_index(i);
            return;
        }
        int prev_bwd_edge_index = current_bwd_edge_index;
        current_bwd_edge_index = i;
        edges->at(prev_bwd_edge_index).set_next_bwd_edge_index(current_bwd_edge_index);
    }

    /**************************************/

    bool has_fwd_edge()
    {
        return first_fwd_edge_index >= 0;
    }

    Edge * get_first_fwd_edge()
    {
        current_fwd_edge_index = first_fwd_edge_index;
        return &edges->at(current_fwd_edge_index);
    }

    bool has_next_fwd_edge()
    {
        return edges->at(current_fwd_edge_index).get_next_fwd_edge_index() >= 0;
    }

    Edge * get_next_fwd_edge()
    {
        if(edges->at(current_fwd_edge_index).get_next_fwd_edge_index() < 0)
            return 0;
        current_fwd_edge_index = edges->at(current_fwd_edge_index).get_next_fwd_edge_index();

        return &edges->at(current_fwd_edge_index);
    }

    /**************************************/

    bool has_bwd_edge()
    {
        return first_bwd_edge_index >= 0;
    }

    Edge * get_first_bwd_edge()
    {
        current_bwd_edge_index = first_bwd_edge_index;
        return &edges->at(current_bwd_edge_index);
    }

    bool has_next_bwd_edge()
    {
        return edges->at(current_bwd_edge_index).get_next_bwd_edge_index() >= 0;
    }

    Edge * get_next_bwd_edge()
    {
        if(edges->at(current_bwd_edge_index).get_next_bwd_edge_index() < 0)
            return 0;
        current_bwd_edge_index = edges->at(current_bwd_edge_index).get_next_bwd_edge_index();
        return &edges->at(current_bwd_edge_index);
    }

    bool contains_bwd_edge(Edge *copy, bool thorough=false)
    {
        if( this->has_bwd_edge() )
        {
            Edge *edge = this->get_first_bwd_edge();
            if(thorough)
            {
                if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                    return true;
            }
            else
            {
                if(*copy == *edge)
                    return true;
            }
            while(this->has_next_bwd_edge())
            {
                edge = this->get_next_bwd_edge();
                if(thorough)
                {
                    if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                        return true;
                }
                else
                {
                    if(*copy == *edge)
                        return true;
                }
            }
        }
        return false;
    }

    void update_bwd_edge_details(Edge *copy, bool thorough=false)
    {
        if( this->has_bwd_edge() )
        {
            Edge *edge = this->get_first_bwd_edge();
            if(thorough)
            {
                if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                {
                    edge->set_branch_count_as_skipped_edge( copy->get_branch_count_as_skipped_edge() );
                    edge->set_branch_count_since_last_used( copy->get_branch_count_since_last_used() );
                    edge->set_branch_distance_since_last_used( copy->get_branch_distance_since_last_used() );
                    edge->set_weight( copy->get_posterior_weight() );
                }
            }
            else
            {
                if(*copy == *edge)
                {
                    edge->set_branch_count_as_skipped_edge( copy->get_branch_count_as_skipped_edge() );
                    edge->set_branch_count_since_last_used( copy->get_branch_count_since_last_used() );
                    edge->set_branch_distance_since_last_used( copy->get_branch_distance_since_last_used() );
                    edge->set_weight( copy->get_posterior_weight() );
                }
            }
            while(this->has_next_bwd_edge())
            {
                edge = this->get_next_bwd_edge();
                if(thorough)
                {
                    if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                    {
                        edge->set_branch_count_as_skipped_edge( copy->get_branch_count_as_skipped_edge() );
                        edge->set_branch_count_since_last_used( copy->get_branch_count_since_last_used() );
                        edge->set_branch_distance_since_last_used( copy->get_branch_distance_since_last_used() );
                        edge->set_weight( copy->get_posterior_weight() );
                    }
                }
                else
                {
                    if(*copy == *edge)
                    {
                        edge->set_branch_count_as_skipped_edge( copy->get_branch_count_as_skipped_edge() );
                        edge->set_branch_count_since_last_used( copy->get_branch_count_since_last_used() );
                        edge->set_branch_distance_since_last_used( copy->get_branch_distance_since_last_used() );
                        edge->set_weight( copy->get_posterior_weight() );
                    }
                }
            }
        }
    }

    bool contains_fwd_edge(Edge *copy, bool thorough=false)
    {
        if( this->has_fwd_edge() )
        {
            Edge *edge = this->get_first_fwd_edge();
            if(thorough)
            {
                if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                    return true;
            }
            else
            {
                if(*copy == *edge)
                    return true;
            }
            while(this->has_next_fwd_edge())
            {
                edge = this->get_next_fwd_edge();
                if(thorough)
                {
                    if(*copy == *edge && copy->get_posterior_weight()==edge->get_posterior_weight())
                        return true;
                }
                else
                {
                    if(*copy == *edge)
                        return true;
                }
            }
        }
        return false;
    }

    void delete_bwd_edge(int edge_ind)
    {
        if(this->has_bwd_edge())
        {
            Edge *edge = this->get_first_bwd_edge();

            if(edge->get_index() == edge_ind)
            {
                if(this->has_next_bwd_edge())
                {
                    Edge *edge2 = this->get_next_bwd_edge();
                    int next_ind = edge2->get_index();
                    current_bwd_edge_index = first_bwd_edge_index = next_ind;
                }
                else
                {
                    current_bwd_edge_index = first_bwd_edge_index = -1;
                }

                return;
            }


            while(this->has_next_bwd_edge())
            {
                int prev_ind = edge->get_index();

                edge = this->get_next_bwd_edge();
                if(edge->get_index() == edge_ind)
                {
                    if(this->has_next_bwd_edge())
                    {
                        edge = this->get_next_bwd_edge();
                        int next_ind = edge->get_index();
                        edges->at(prev_ind).set_next_bwd_edge_index(next_ind);
                    }
                    else
                    {
                        edges->at(prev_ind).set_next_bwd_edge_index(-1);
                    }

                }
            }
        }
    }

    void delete_fwd_edge(int edge_ind)
    {
        if(this->has_fwd_edge())
        {
            Edge *edge = this->get_first_fwd_edge();

            if(edge->get_index() == edge_ind)
            {
                if(this->has_next_fwd_edge())
                {
                    edge = this->get_next_fwd_edge();
                    int next_ind = edge->get_index();
                    current_fwd_edge_index = first_fwd_edge_index = next_ind;
                }
                else
                {
                    current_fwd_edge_index = first_fwd_edge_index = -1;
                }

                return;
            }

            while(this->has_next_fwd_edge())
            {
                int prev_ind = edge->get_index();

                edge = this->get_next_fwd_edge();
                if(edge->get_index() == edge_ind)
                {
                    if(this->has_next_fwd_edge())
                    {
                        edge = this->get_next_fwd_edge();
                        int next_ind = edge->get_index();
                        edges->at(prev_ind).set_next_fwd_edge_index(next_ind);
                    }
                    else
                    {
                        edges->at(prev_ind).set_next_fwd_edge_index(-1);
                    }
                }
            }
        }
    }

    friend ostream& operator<< (ostream &out, Site& b)
    {
        out<<"index: "<<b.index<<"; character_state: "<<b.character_state<<"; site_type: "<<b.site_type<<"; path_state: "<<b.path_state;
        return out;
    }

    bool operator==(const Site& b)
    {
        return children==b.children;
    }

    bool operator!=(const Site& b)
    {
        return !(children==b.children);
    }

    bool operator<(const Site& b)
    {
        return children<b.children;
    }

    bool operator>(const Site& b)
    {
        return children>b.children;
    }

    static bool comesBefore(const Site& a,const Site& b)
    {
        return  (a.unique_index.left_index<b.unique_index.left_index && a.unique_index.right_index<=b.unique_index.right_index) ||
                (a.unique_index.left_index<=b.unique_index.left_index && a.unique_index.right_index<b.unique_index.right_index);
    }
};

/**************************************/


class Sequence
{
    int curr_site_index;
    int prev_site_index;
    int curr_edge_index;

    vector<Site> sites;
    vector<Edge> edges;
    string full_char_alphabet;
    int data_type;
    string gap_symbol;

    vector<Unique_index> unique_index;

    bool read_sequence;
    bool has_read_descendants;
    bool terminal_sequence;
    string gapped_seq;
    string unaligned_seq;
    string dna_seq;

    int num_duplicates;
public:

    Sequence(Fasta_entry &seq_entry,const int data_type,bool gapped = false, bool no_trimming=false, bool turn_revcomp=false);
    Sequence(const int length,const int data_type, string gapped_s="");

    bool is_terminal_sequence() { return terminal_sequence; }
    void is_terminal_sequence(bool t) { terminal_sequence = t; }

    bool is_read_sequence() { return read_sequence; }
    void is_read_sequence(bool t) { read_sequence = t; }

    bool is_read_descendants() { return has_read_descendants; }
    void is_read_descendants(bool t) { has_read_descendants = t; }

    int get_num_duplicates() { return num_duplicates; }

    void initialise_indeces() {
        curr_site_index = prev_site_index = curr_edge_index = 0;
    }

    void set_gap_symbol(char c) {  gap_symbol = string(1,c); }
    void set_gap_symbol(string c) {  gap_symbol = c; }
    string get_gap_symbol() {  return gap_symbol; }

    int get_data_type() { return data_type; }

    int get_new_site_index() {curr_site_index++; return curr_site_index; }
    int get_new_edge_index() {curr_edge_index++; return curr_edge_index; }

    int get_current_site_index() { return curr_site_index; }
    int get_previous_site_index() { return prev_site_index; }
    int get_current_edge_index() { return curr_edge_index; }

    void push_back_site(Site site)
    {
        site.set_index(sites.size());
        sites.push_back(site);

        prev_site_index = sites.size()-2;
        curr_site_index = sites.size()-1;
    }

    void push_back_edge(Edge edge)
    {
        edge.set_index(edges.size());
        edges.push_back(edge);

        curr_edge_index = edges.size()-1;
    }

    bool contains_edge(Edge edge)
    {
        vector<Edge>::iterator it = edges.begin();
        for(;it!=edges.end();it++)
        {
            if(*it==edge)
                return true;
        }
        return false;
    }

    Site *get_previous_site()
    {
        return &sites.at(prev_site_index);
    }

    Site *get_current_site()
    {
        return &sites.at(curr_site_index);
    }

    int get_bwd_edge_index_at_site(int site,Edge *copy)
    {
        if(!sites.at(site).has_bwd_edge())
            return -1;

        Edge *edge = sites.at(site).get_first_bwd_edge();
        if(*edge == *copy)
            return edge->get_index();

        while(sites.at(site).has_next_bwd_edge())
        {
            edge = sites.at(site).get_next_bwd_edge();
            if(*edge == *copy)
                return edge->get_index();
        }
        return -1;
    }

    bool contains_this_bwd_edge_at_site(int site, Edge *copy)
    {
        if(this->get_bwd_edge_index_at_site(site,copy)<0)
            return false;

        return true;
    }

    int get_fwd_edge_index_at_site(int site,Edge *copy)
    {
        if(!sites.at(site).has_fwd_edge())
            return -1;

        Edge *edge = sites.at(site).get_first_fwd_edge();
        if(*edge == *copy)
            return edge->get_index();

        while(sites.at(site).has_next_fwd_edge())
        {
            edge = sites.at(site).get_next_fwd_edge();
            if(*edge == *copy)
                return edge->get_index();
        }
        return -1;
    }

    bool contains_this_fwd_edge_at_site(int site, Edge *copy)
    {
        if(this->get_fwd_edge_index_at_site(site,copy)<0)
            return false;

        return true;
    }

    string *get_unaligned_sequence();
    string get_sequence_string(bool with_gaps);

    string print_sequence(vector<Site> *sites);
    string print_sequence() { return this->print_sequence(this->get_sites()); }
    string print_path(vector<Site> *sites);
    string print_path(){ return this->print_path(this->get_sites()); }

    string *get_gapped_sequence() { return &gapped_seq; }
    string *get_dna_sequence() { return &dna_seq; }


    void create_default_sequence(Fasta_entry &seq_entry);
    void create_codon_sequence(Fasta_entry &seq_entry);
    void create_fastq_sequence(Fasta_entry &seq_entry, bool no_trimming=false);
    void create_graph_sequence(Fasta_entry &seq_entry);

    vector<Site> *get_sites() { return &sites; }
    vector<Edge> *get_edges() { return &edges; }

    int sites_length() { return sites.size(); }
    int edges_length() { return edges.size(); }

    Site *get_site_at(int i) { return &sites.at(i); }
    Edge *get_first_bwd_edge_at(int i) { return sites.at(i).get_first_bwd_edge(); }

    string get_full_alphabet() { return full_char_alphabet; }

    void delete_all_bwd_edges_at_site(int index)
    {
        Site *site = this->get_site_at(index);

        if(site->has_bwd_edge())
        {
            Edge *edge = site->get_first_bwd_edge();
            this->get_site_at(edge->get_start_site_index())->delete_fwd_edge(edge->get_index());

            while(site->has_next_bwd_edge())
            {
                edge = site->get_next_bwd_edge();
                this->get_site_at(edge->get_start_site_index())->delete_fwd_edge(edge->get_index());
            }
        }
        site->set_first_bwd_edge_index(-1);
    }

    void delete_all_fwd_edges_at_site(int index)
    {
        Site *site = this->get_site_at(index);

        if(site->has_fwd_edge())
        {
            Edge *edge = site->get_first_fwd_edge();
            this->get_site_at(edge->get_end_site_index())->delete_bwd_edge(edge->get_index());

            while(site->has_next_fwd_edge())
            {
                edge = site->get_next_fwd_edge();
                this->get_site_at(edge->get_end_site_index())->delete_bwd_edge(edge->get_index());
            }
        }
        site->set_first_fwd_edge_index(-1);
    }

    void initialise_unique_index()
    {
        unique_index.clear();
        int prev_left = 0;
        int prev_right = 0;
        for(int i=0;i<this->sites_length();i++)
        {
            Site *ts = this->get_site_at(i);
            int this_left = ts->get_children()->left_index;
            int this_right = ts->get_children()->right_index;
            if(this_left>0 && this_right>0)
            {
                unique_index.push_back( Unique_index(this_left,this_right, Unique_index::match,i) );
                ts->set_unique_index(this_left,this_right, Unique_index::match,i);

                prev_left = this_left;
                prev_right = this_right;
            }
            else if(this_left>0)
            {
                unique_index.push_back( Unique_index(this_left,prev_right, Unique_index::xgap,i) );
                ts->set_unique_index(this_left,prev_right, Unique_index::xgap,i);

                prev_left = this_left;
            }
            else if(this_right>0)
            {
                unique_index.push_back( Unique_index(prev_left,this_right, Unique_index::match,i) );
                ts->set_unique_index(prev_left,this_right, Unique_index::match,i);

                prev_right = this_right;
            }
            else if(ts->get_site_type()==Site::start_site)
            {
                unique_index.push_back( Unique_index(0,0, Unique_index::match,i) );
                ts->set_unique_index(0,0, Unique_index::match,i);
            }
            else if(ts->get_site_type()==Site::stop_site)
            {
                unique_index.push_back( Unique_index(prev_left+1,prev_right+1, Unique_index::match,i) );
                ts->set_unique_index(prev_left+1,prev_right+1, Unique_index::match,i);
            }
        }
    }

    void add_term_in_unique_index(Unique_index ind)
    {
        unique_index.push_back( ind );
    }

    vector<Unique_index> *get_unique_index() { return &unique_index; }

    Unique_index *get_unique_index_at(int i) { return &unique_index.at(i); }

    bool is_unique_index_ordered()
    {
        for(int i=0;i<(int) this->unique_index.size()-1;i++)
        {
            if( ! (this->unique_index.at(i) < this->unique_index.at(i+1) ) )
                return false;
        }
        return true;
    }

    bool unique_index_contains_term(Unique_index *s)
    {
        for(int i=0;i<(int) this->unique_index.size();i++)
        {
            if(this->unique_index.at(i) == (*s))
                return true;
        }

        return false;
    }

    int unique_index_of_term(Unique_index *s)
    {
        for(int i=0;i<(int) this->unique_index.size();i++)
        {
            if(this->unique_index.at(i) == (*s))
                return unique_index.at(i).site_index;
        }

        return -1;
    }

    void copy_site_details(Site* original, Site *copy)
    {
        copy->set_children( original->get_children()->left_index, original->get_children()->right_index );
        copy->set_unique_index( original->get_unique_index()->left_index, original->get_unique_index()->right_index,
                                original->get_unique_index()->match_state, original->get_unique_index()->site_index );
        copy->set_state( original->get_state() );
        copy->set_path_state( original->get_path_state() );
        copy->set_branch_count_since_last_used( original->get_branch_count_since_last_used() );
        copy->set_branch_distance_since_last_used( original->get_branch_distance_since_last_used() );
    }

    void copy_edge_details(Edge* original, Edge *copy)
    {
        copy->set_branch_count_as_skipped_edge( original->get_branch_count_as_skipped_edge() );
        copy->set_branch_count_since_last_used( original->get_branch_count_since_last_used() );
        copy->set_branch_distance_since_last_used( original->get_branch_distance_since_last_used() );
        copy->set_weight( original->get_posterior_weight() );
    }

    void sort_sites_vector()
    {
        sort(sites.begin(),sites.end(),Site::comesBefore);
    }
    
    void remap_edges_vector()
    {
        vector<int> new_site_index;
        new_site_index.resize( this->sites_length() );

        for(int i=0;i<this->sites_length();i++)
            new_site_index.at( this->get_site_at(i)->get_index() ) = i;

        for(int i=0;i<(int) edges.size();i++)
        {
             int s = edges.at(i).get_start_site_index();
             int e = edges.at(i).get_end_site_index();

             edges.at(i).set_start_site_index( new_site_index.at(s) );
             edges.at(i).set_end_site_index( new_site_index.at(e) );
        }

        for(int i=0;i<this->sites_length();i++)
            this->get_site_at(i)->set_index(i);

    }

};
}

#endif // SEQUENCE_H
