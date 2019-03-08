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

#include <iostream>
#include <map>
#include "utils/settings.h"
#include "utils/model_factory.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"
#include "utils/text_utils.h"
#include "main/sequence.h"
#include <iomanip>
#include <sstream>

using namespace std;
using namespace ppa;

Sequence::Sequence(Fasta_entry &seq_entry,const int data_t,bool gapped, bool no_trimming, bool turn_revcomp)
{
    data_type = data_t;

    this->set_gap_symbol("-");
    if( ( data_type == Model_factory::dna && Settings_handle::st.is("codons") ) || data_type == Model_factory::codon)
        this->set_gap_symbol("---");

    if(gapped)
    {
        gapped_seq = seq_entry.sequence;

        /*
        for (string::iterator si = seq_entry.sequence.begin();si != seq_entry.sequence.end();)
            if(*si == '-')
                seq_entry.sequence.erase(si);
            else
                si++;
        */
        stringstream seqss;
        for(int i = 0; i < (int) seq_entry.sequence.size(); i++)
        {
            if(seq_entry.sequence[i] != '-')
            {
                seqss << seq_entry.sequence[i];
            }
        }
        seq_entry.sequence = seqss.str();
    }
    else
    {
        gapped_seq = "";
    }

    unaligned_seq = seq_entry.sequence;

    if(seq_entry.dna_sequence.length()>0 && ( Settings_handle::st.is("translate") || Settings_handle::st.is("mt-translate") ||
            Settings_handle::st.is("find-orfs") ) )
    {

        stringstream seqss;
        for(int i = 0; i < (int) seq_entry.dna_sequence.size(); i++)
        {
            if(seq_entry.dna_sequence[i] != '-')
            {
                seqss << seq_entry.dna_sequence[i];
            }
        }

        dna_seq = seqss.str();
    }

    if(turn_revcomp)
    {
        string revcomp = seq_entry.sequence;

        reverse(revcomp.begin(), revcomp.end());

        Text_utils tu;

        tu.replace_all(revcomp,'A','Z');
        tu.replace_all(revcomp,'T','A');
        tu.replace_all(revcomp,'Z','T');
        tu.replace_all(revcomp,'C','Z');
        tu.replace_all(revcomp,'G','C');
        tu.replace_all(revcomp,'Z','G');
        tu.replace_all(revcomp,'R','Z');
        tu.replace_all(revcomp,'Y','R');
        tu.replace_all(revcomp,'Z','Y');
        tu.replace_all(revcomp,'K','Z');
        tu.replace_all(revcomp,'M','K');
        tu.replace_all(revcomp,'Z','M');
        tu.replace_all(revcomp,'B','Z');
        tu.replace_all(revcomp,'V','B');
        tu.replace_all(revcomp,'Z','V');
        tu.replace_all(revcomp,'D','Z');
        tu.replace_all(revcomp,'H','D');
        tu.replace_all(revcomp,'Z','H');

        seq_entry.sequence = revcomp;

    }

    this->initialise_indeces();

    full_char_alphabet = Model_factory::get_dna_full_char_alphabet();
    if(data_type == Model_factory::protein)
        full_char_alphabet = Model_factory::get_protein_full_char_alphabet();


    sites.reserve(seq_entry.sequence.size()+2);
    edges.reserve(seq_entry.sequence.size()+3);

    if( seq_entry.quality != "" && !Settings_handle::st.is("no-fastq") )
        this->create_fastq_sequence(seq_entry, no_trimming);

    else if( seq_entry.edges.size()>0 )
        this->create_graph_sequence(seq_entry);

    else
        if( ( data_type == Model_factory::dna && Settings_handle::st.is("codons") ) || data_type == Model_factory::codon)
            this->create_codon_sequence(seq_entry);
        else
            this->create_default_sequence(seq_entry);

    terminal_sequence = true;
    read_sequence = false;
    has_read_descendants = false;
    num_duplicates = seq_entry.num_duplicates;

    if(Settings::noise>5)
    {
        this->print_sequence(&sites);
    }

}

void Sequence::create_default_sequence(Fasta_entry &seq_entry)
{

    Site first_site( &edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_empty_children();
    this->push_back_site(first_site);

    int in_row = 1;
    int prev_row = 1;
    int prev_state = -1;

    Edge first_edge( -1,this->get_current_site_index() );
    this->push_back_edge(first_edge);

    string::iterator si = seq_entry.sequence.begin();
    string::iterator qi = seq_entry.quality.begin();

    for(;si!=seq_entry.sequence.end();si++,qi++)
    {

        if(*si=='0')
        {
            continue;
        }

        Site site( &edges );
        site.set_state( full_char_alphabet.find( *si ) );
        site.set_symbol( *si );
        site.set_empty_children();
        this->push_back_site(site);

//        Edge edge( this->get_previous_site_index(),this->get_current_site_index() );
//        this->push_back_edge(edge);

//        this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
//        this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

        // Check for homopolymers
        if( site.get_state() == prev_state)
        {
            in_row++;
            prev_row = 1;
        }
        else
        {
            prev_row = in_row;
            in_row = 1;
            prev_state = site.get_state();
        }

        // If 454 data, correct for homopolymer error
        //
        if(Settings_handle::st.is("454") &&  prev_row > 2  )
        {
            // first edge
            float weight = 1;

            Edge edge( this->get_previous_site_index(),this->get_current_site_index(), weight );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

            if( prev_row < 5 )
            {
                // second edge
                int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();
//                Edge edge_2( prev_ind ,this->get_current_site_index(), 1.0-weight );
                Edge edge_2( prev_ind ,this->get_current_site_index(), 0.9 );
                this->push_back_edge(edge_2);

                this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

            }
            else
            {
                // second edge
                int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();
//                Edge edge_2( prev_ind ,this->get_current_site_index(), 1.0-weight-0.02 );
                Edge edge_2( prev_ind ,this->get_current_site_index(), 0.9 );
                this->push_back_edge(edge_2);

                this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

                // third edge
                int prev_prev_ind = get_site_at(prev_ind)->get_first_bwd_edge()->get_start_site_index();
//                Edge edge_3( prev_prev_ind ,this->get_current_site_index(), 0.02 );
                Edge edge_3( prev_prev_ind ,this->get_current_site_index(), 0.9 );
                this->push_back_edge(edge_3);

                this->get_site_at(prev_prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

            }
        }

        // If 454 data, correct for homopolymer error
        //
        else if(Settings_handle::st.is("homopolymer") &&  prev_row >= 2  )
        {
            // first edge
            float weight = 0.25;

            Edge edge( this->get_previous_site_index(),this->get_current_site_index(), 1 );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

            int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();

            while(prev_row>=2)
            {
                Edge edge_2( prev_ind ,this->get_current_site_index(), weight );
                this->push_back_edge(edge_2);

                this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

                prev_ind = get_site_at(prev_ind)->get_first_bwd_edge()->get_start_site_index();

                prev_row--;
            }
        }

        else
        {
            Edge edge( this->get_previous_site_index(),this->get_current_site_index() );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );
        }


    }

    Site last_site( &edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_empty_children();
    this->push_back_site(last_site);

    Edge last_edge( this->get_previous_site_index(),this->get_current_site_index() );
    this->push_back_edge(last_edge);

    this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
    this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

}


void Sequence::create_codon_sequence(Fasta_entry &seq_entry)
{
    this->set_gap_symbol("---");

    Site first_site( &edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_empty_children();
    this->push_back_site(first_site);

    Edge first_edge( -1,this->get_current_site_index() );
    this->push_back_edge(first_edge);

    vector<string> *fca = Model_factory::get_codon_character_alphabet();
    map<string,int> codons;
    int count = 0;
    for(vector<string>::iterator it = fca->begin();it != fca->end();it++)
        codons.insert(make_pair(*it,count++));

    for(int i=0;i<(int)seq_entry.sequence.length();i+=3)
    {
        int state = 61;

        string codon = seq_entry.sequence.substr(i,3);
        map<string,int>::iterator fi = codons.find(codon);
        if(fi!=codons.end())
            state = fi->second;
        else
            codon = "NNN";

        Site site( &edges );
        site.set_state( state );
        site.set_symbol( codon );
        site.set_empty_children();
        this->push_back_site(site);

        Edge edge( this->get_previous_site_index(),this->get_current_site_index() );
        this->push_back_edge(edge);

        this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
        this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );
    }

    Site last_site( &edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_empty_children();
    this->push_back_site(last_site);

    Edge last_edge( this->get_previous_site_index(),this->get_current_site_index() );
    this->push_back_edge(last_edge);

    this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
    this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

}
void Sequence::create_fastq_sequence(Fasta_entry &seq_entry, bool no_trimming)
{

    Site first_site( &edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_empty_children();
    this->push_back_site(first_site);

    Edge first_edge( -1,this->get_current_site_index() );
    this->push_back_edge(first_edge);

    int quality_threshold = Settings_handle::st.get("qscore-minimum").as<int>();

    if(no_trimming)
        quality_threshold = 0;

    int in_row = 1;
    int prev_row = 1;
    int prev_state = -1;

    string::iterator si = seq_entry.sequence.begin();
    string::iterator tsi;
    string::iterator qi = seq_entry.quality.begin();

    int site_qscore = quality_threshold;

    for(;si!=seq_entry.sequence.end();si++,qi++)
    {
        // this site is a paired-read break point
        //
        if(*si=='0')
        {
            continue;
        }


        int prev_site_qscore = site_qscore;

        Site site( &edges );
        site.set_empty_children();


        // check if the previous site was a paired-read break point
        //
        tsi = si;
        if(tsi!=seq_entry.sequence.begin())
        {
            tsi--;
            if(*tsi=='0')
            {
                site.set_site_type(Site::break_stop_site);
            }
        }

        // check if the next site will be a paired-read break point
        //
        tsi = si;
        tsi++;
        if(tsi!=seq_entry.sequence.end())
        {
            if(*tsi=='0')
            {
                site.set_site_type(Site::break_start_site);
            }
        }

        site_qscore = static_cast<int>(*qi)-33;

        if(site_qscore < quality_threshold)
        {
            site.set_state( full_char_alphabet.find( 'N' ) );
            site.set_symbol( tolower(*si) );
        }
        else
        {
            site.set_state( full_char_alphabet.find( *si ) );
            site.set_symbol( *si );
        }

        this->push_back_site(site);

        // Check for homopolymers
        if( site.get_state() == prev_state)
        {
            in_row++;
            prev_row = 1;
        }
        else
        {
            prev_row = in_row;
            in_row = 1;
            prev_state = site.get_state();
        }

        // If 454 data, correct for homopolymer error
        //
        if(Settings_handle::st.is("454") && ( prev_row > 2 ||  prev_site_qscore < quality_threshold ) )
        {
            // first edge
            float weight = 0.9;
            if(prev_site_qscore < quality_threshold)
                weight = 0.6;

            Edge edge( this->get_previous_site_index(),this->get_current_site_index(), weight );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

            if( prev_row < 5 )
            {
                // second edge
                int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();
                Edge edge_2( prev_ind ,this->get_current_site_index(), 1.0-weight );
                this->push_back_edge(edge_2);

                this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

            }
            else
            {
                // second edge
                int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();
                Edge edge_2( prev_ind ,this->get_current_site_index(), 1.0-weight-0.02 );
                this->push_back_edge(edge_2);

                this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

                // third edge
                int prev_prev_ind = get_site_at(prev_ind)->get_first_bwd_edge()->get_start_site_index();
                Edge edge_3( prev_prev_ind ,this->get_current_site_index(), 0.02 );
                this->push_back_edge(edge_3);

                this->get_site_at(prev_prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
                this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );

            }

        }

        else if( Settings_handle::st.is("allow-skip-low-qscore") && prev_site_qscore < quality_threshold )
        {
            Edge edge( this->get_previous_site_index(),this->get_current_site_index(), 0.6 );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

            // second edge
            int prev_ind = this->get_previous_site()->get_first_bwd_edge()->get_start_site_index();
            Edge edge_2( prev_ind ,this->get_current_site_index(), 0.4 );
            this->push_back_edge(edge_2);

            this->get_site_at(prev_ind)->add_new_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->add_new_bwd_edge_index( this->get_current_edge_index() );
        }

        // All other data
        else
        {
            Edge edge( this->get_previous_site_index(),this->get_current_site_index() );
            this->push_back_edge(edge);

            this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
            this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );
        }
    }

    Site last_site( &edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_empty_children();
    this->push_back_site(last_site);

    Edge last_edge( this->get_previous_site_index(),this->get_current_site_index() );
    this->push_back_edge(last_edge);

    this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
    this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );
}

void Sequence::create_graph_sequence(Fasta_entry &seq_entry)
{
    Site first_site( &edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_empty_children();
    this->push_back_site(first_site);

    Edge first_edge( -1,this->get_current_site_index() );
    this->push_back_edge(first_edge);

    string::iterator si = seq_entry.sequence.begin();
    string::iterator qi = seq_entry.quality.begin();

    for(;si!=seq_entry.sequence.end();si++,qi++)
    {

        Site site( &edges );
        site.set_state( full_char_alphabet.find( *si ) );
        site.set_symbol( *si );
        site.set_empty_children();
        this->push_back_site(site);

    }

    Site last_site( &edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_empty_children();
    this->push_back_site(last_site);

    for(int i=0; i<(int)seq_entry.edges.size(); i++)
    {
        int s = seq_entry.edges.at(i).start_site;
        int e = seq_entry.edges.at(i).end_site;
        double w = seq_entry.edges.at(i).weight;
        Edge edge( s, e, w );
        this->push_back_edge(edge);

        if( this->get_site_at(s)->has_fwd_edge() )
            this->get_site_at(s)->add_new_fwd_edge_index( this->get_current_edge_index() );
        else
            this->get_site_at(s)->set_first_fwd_edge_index( this->get_current_edge_index() );

        if( this->get_site_at(e)->has_bwd_edge() )
            this->get_site_at(e)->add_new_bwd_edge_index( this->get_current_edge_index() );
        else
            this->get_site_at(e)->set_first_bwd_edge_index( this->get_current_edge_index() );
    }
}


Sequence::Sequence(const int length,const int data_t, string gapped_s)
{
    data_type = data_t;

    this->set_gap_symbol("-");
    if( ( data_type == Model_factory::dna && Settings_handle::st.is("codons") ) || data_type == Model_factory::codon)
        this->set_gap_symbol("---");

    gapped_seq = gapped_s;
    unaligned_seq = "";

    this->initialise_indeces();

    sites.reserve(length+2);
    edges.reserve(length+3);

    full_char_alphabet = Model_factory::get_dna_full_char_alphabet();
    if(data_type == Model_factory::protein)
        full_char_alphabet = Model_factory::get_protein_full_char_alphabet();

    read_sequence = false;
    terminal_sequence = false;
    has_read_descendants = false;

}



string Sequence::print_sequence(vector<Site> *sites)
{

    vector<string> *alphabet = 0;
    if( ( data_type == Model_factory::dna && Settings_handle::st.is("codons") ) || data_type == Model_factory::codon)
        alphabet = Model_factory::get_codon_full_character_alphabet();
    else if(data_type == Model_factory::dna)
        alphabet = Model_factory::get_dna_full_character_alphabet();
    else if(data_type == Model_factory::protein)
        alphabet = Model_factory::get_protein_full_character_alphabet();

    stringstream ss;
    ss<<endl;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);
        ss<<i<<": ";
        if(tsite->get_site_type()==Site::real_site)
            ss<<tsite->get_index()<<" "<<setw(2)<<alphabet->at(tsite->get_state());
        else
            ss<<tsite->get_index()<<" +";

        ss<<"("<<tsite->get_unique_index()->left_index<<","<<tsite->get_unique_index()->right_index<<") ";
        ss<<"["<<tsite->get_path_state()<<"] ["<<tsite->get_site_type()<<"] ";
        ss<<"\t";
        ss << setprecision (2);

        if(true)
        {
            if(tsite->has_fwd_edge())
            {
                Edge *tedge = tsite->get_first_fwd_edge();
                ss<<" F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                        <<" "<<scientific<<tedge->get_posterior_weight()<<fixed<<" "<<tedge->get_branch_count_since_last_used()<<" "
                        <<tedge->get_branch_distance_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<"]";
                while(tsite->has_next_fwd_edge())
                {
                    tedge = tsite->get_next_fwd_edge();
                    ss<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                        <<" "<<scientific<<tedge->get_posterior_weight()<<fixed<<" "<<tedge->get_branch_count_since_last_used()<<" "
                        <<tedge->get_branch_distance_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<"]";
                }
            }
            ss<<"; \t";
        }
        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            ss<<"B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                    <<" "<<scientific<<tedge->get_posterior_weight()<<fixed<<" "<<tedge->get_branch_count_since_last_used()<<" "
                    <<tedge->get_branch_distance_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<"]";
            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                ss<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                    <<" "<<scientific<<tedge->get_posterior_weight()<<fixed<<" "<<tedge->get_branch_count_since_last_used()<<" "
                    <<tedge->get_branch_distance_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<"]";
            }
        }
        ss << setprecision (4);
        
        ss<<"\n";
    }

    return ss.str();
}

/***************************************************************************/

string *Sequence::get_unaligned_sequence()
{
    if(unaligned_seq.length()>0) { return &unaligned_seq; }

    unaligned_seq = gapped_seq;

    for (string::iterator si = unaligned_seq.begin();si != unaligned_seq.end();)
        if(*si == '-')
            unaligned_seq.erase(si);
        else
            si++;

    return &unaligned_seq;
}

string Sequence::get_sequence_string(bool with_gaps)
{

    int seq_length = this->sites_length();
    string seq;

    bool codonstr = ( ( this->get_data_type()==Model_factory::dna && Settings_handle::st.is("codons") ) || this->get_data_type() == Model_factory::codon);

    if(this->is_terminal_sequence())
    {
        for(int j=1;j<seq_length-1;j++)
            seq += this->get_site_at(j)->get_symbol();
    }
    else
    {
        for(int j=1;j<seq_length-1;j++)
        {
            string c = Model_factory::get_ancestral_character_alphabet_at( this->get_site_at(j)->get_state() );

            int pstate = this->get_site_at(j)->get_path_state();
            int ptype  = this->get_site_at(j)->get_site_type();

            if( pstate != Site::xskipped && pstate != Site::yskipped && ptype != Site::non_real)
                seq += c;
            else
                if(with_gaps)
                {
                    if(codonstr)
                        seq +=  "---";
                    else
                        seq +=  '-';
                }
        }
    }

    return seq;
}

/***************************************************************************/


string Sequence::print_path(vector<Site> *sites)
{
    stringstream ss;
    ss<<endl;

   for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

        ss<<i<<" "<<tsite->get_state()<<" "<<endl;

        int ps = tsite->path_state;
        switch(ps)
        {
            case Site::matched:
                ss<<"M";
                continue;
            case Site::xgapped:
                ss<<"X";
                continue;
            case Site::ygapped:
                ss<<"Y";
                continue;
            case Site::xskipped:
                ss<<"x";
                continue;
            case Site::yskipped:
                ss<<"y";
                continue;
            default:
                ss<<"o";
                continue;
        }
        ss<<": ";

        if(tsite->get_site_type()==Site::real_site)
            ss<<tsite->get_index()<<" "<<full_char_alphabet.at(tsite->get_state());
        else
            ss<<tsite->get_index()<<" +";

    }
    return ss.str();
}
