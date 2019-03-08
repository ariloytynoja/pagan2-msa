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

#include "exonerate_queries.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/regex.hpp>
#include <sys/stat.h>

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;
using namespace ppa;

Exonerate_queries::Exonerate_queries()
{

}

bool Exonerate_queries::test_executable()
{

    int status = -1;

    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    exoneratepath = epath;
    epath = epath+"exonerate.exe > /dev/null 2>/dev/null";
    status = system(epath.c_str());

    #else

    char path[200];
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';
    epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #endif

    exoneratepath = epath;
    epath = epath+"exonerate >/dev/null 2>/dev/null";
    status = system(epath.c_str());

    #endif

    if(WEXITSTATUS(status) == 1)
        return true;

    exoneratepath = "";
    status = system("`exonerate  >/dev/null 2>/dev/null`");

    return WEXITSTATUS(status) == 1;
}

bool Exonerate_queries::split_sugar_string(const string& row,hit *h)
{

    const boost::regex pattern("sugar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );
    }

    return valid;
}

bool Exonerate_queries::split_vulgar_string(const string& row,hit *h)
{

    const boost::regex pattern("vulgar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\\s+(.+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );

    }



    return valid;
}


void Exonerate_queries::write_exonerate_input(string *str1, string *str2, int *r)
{
    // create exonerate input
    //
    ofstream q_output;
    ofstream t_output;

    while(true)
    {

        stringstream q_name;
        stringstream t_name;

        string tmp_dir = this->get_temp_dir();

        q_name <<tmp_dir<<"q"<<*r<<".fas";
        t_name <<tmp_dir<<"t"<<*r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        *r = rand();
    }

    t_output<<">t\n"<<*str2<<endl;
    t_output.close();

    q_output<<">q\n"<<*str1<<endl;
    q_output.close();

    return;
}

void Exonerate_queries::write_exonerate_input(map<string,string> *target_sequences, vector<Fasta_entry> *reads, int *r)
{

    // create exonerate input
    //
    ofstream q_output;
    ofstream t_output;

    while(true)
    {

        stringstream q_name;
        stringstream t_name;

        string tmp_dir = this->get_temp_dir();

        q_name <<tmp_dir<<"q"<<*r<<".fas";
        t_name <<tmp_dir<<"t"<<*r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        *r = rand();
    }

    vector<Fasta_entry>::iterator it = reads->begin();
    for(;it!=reads->end();it++)
    {
        if(Settings_handle::st.is("score-as-dna") && it->dna_sequence.length()>0)
            q_output<<">"<<it->name<<endl<<it->dna_sequence<<endl;
        else
            q_output<<">"<<it->name<<endl<<it->sequence<<endl;
    }
    q_output.close();


    map<string,string>::iterator it2 = target_sequences->begin();
    for(;it2!=target_sequences->end();it2++)
        t_output<<">"<<it2->first<<endl<<it2->second<<endl;

    t_output.close();

    return;
}

void Exonerate_queries::write_exonerate_input(map<string,string> *target_sequences, Fasta_entry *read, int *r)
{

    // create exonerate input
    //
    ofstream q_output;
    ofstream t_output;

    while(true)
    {

        stringstream q_name;
        stringstream t_name;

        string tmp_dir = this->get_temp_dir();

        q_name <<tmp_dir<<"q"<<*r<<".fas";
        t_name <<tmp_dir<<"t"<<*r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        *r = rand();
    }

    q_output<<">"<<read->name<<endl<<read->sequence<<endl;

    q_output.close();


    map<string,string>::iterator it2 = target_sequences->begin();
    for(;it2!=target_sequences->end();it2++)
        t_output<<">"<<it2->first<<endl<<it2->second<<endl;

    t_output.close();

    return;
}

void Exonerate_queries::write_exonerate_input(Node *root, vector<Fasta_entry> *reads, map<string,string> *names, int *r)
{

    // create exonerate input
    //
    ofstream q_output;
    ofstream t_output;

    while(true)
    {

        stringstream q_name;
        stringstream t_name;

        string tmp_dir = this->get_temp_dir();

        q_name <<tmp_dir<<"q"<<*r<<".fas";
        t_name <<tmp_dir<<"t"<<*r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        *r = rand();
    }

    vector<Fasta_entry>::iterator it = reads->begin();
    for(;it!=reads->end();it++)
    {
        if(Settings_handle::st.is("score-as-dna") && it->dna_sequence.length()>0)
            q_output<<">"<<it->name<<endl<<it->dna_sequence<<endl;
        else
            q_output<<">"<<it->name<<endl<<it->sequence<<endl;
    }
    q_output.close();


    map<string,string*> unaligned_seqs;

    if(Settings_handle::st.is("score-as-dna"))
        root->get_dna_sequences(&unaligned_seqs);
    else
        root->get_unaligned_sequences(&unaligned_seqs);

    map<string,string*>::iterator it2 = unaligned_seqs.begin();
    for(;it2!=unaligned_seqs.end();it2++)
    {
        if(names->find(it2->first) != names->end())
        {
            t_output<<">"<<it2->first<<endl<<*(it2->second)<<endl;
        }
    }
    t_output.close();

    return;
}

void Exonerate_queries::write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int *r)
{

    // create exonerate input
    //
    ofstream q_output;
    ofstream t_output;

    while(true)
    {

        stringstream q_name;
        stringstream t_name;

        string tmp_dir = this->get_temp_dir();

        q_name <<tmp_dir<<"q"<<*r<<".fas";
        t_name <<tmp_dir<<"t"<<*r<<".fas";

        ifstream q_file(q_name.str().c_str());
        ifstream t_file(t_name.str().c_str());

        if(!q_file && !t_file)
        {
            q_output.open( q_name.str().c_str(), (ios::out) );
            t_output.open( t_name.str().c_str(), (ios::out) );

            break;
        }
        *r = rand();
    }

    if(Settings_handle::st.is("score-as-dna") && read->dna_sequence.length()>0)
        q_output<<">"<<read->name<<endl<<read->dna_sequence<<endl;
    else
        q_output<<">"<<read->name<<endl<<read->sequence<<endl;
    q_output.close();


    if(Settings_handle::st.is("score-as-dna"))
    {
        map<string,string*> dna_seqs;
        root->get_dna_sequences(&dna_seqs);

        map<string,string*>::iterator it = dna_seqs.begin();
        for(;it!=dna_seqs.end();it++)
        {
            if(names->find(it->first) != names->end())
            {
                string seq = *(it->second);

                for (string::iterator si = seq.begin();si != seq.end();)
                    if(*si == '-')
                        seq.erase(si);
                    else
                        si++;
                t_output<<">"<<it->first<<endl<<seq<<endl;
            }
        }

    }
    else
    {
        vector<Fasta_entry> aligned_sequences;
        root->get_alignment(&aligned_sequences,true);

        vector<Fasta_entry>::iterator it = aligned_sequences.begin();
        for(;it!=aligned_sequences.end();it++)
        {
            if(names->find(it->name) != names->end())
            {
                string seq = it->sequence;

                for (string::iterator si = seq.begin();si != seq.end();)
                    if(*si == '-')
                        seq.erase(si);
                    else
                        si++;
                t_output<<">"<<it->name<<endl<<seq<<endl;
            }
        }
    }
    t_output.close();

    return;
}


void Exonerate_queries::preselect_targets(map<string,string> *target_sequences, vector<Fasta_entry> *reads, map<string, string> *selected_sequences,map<string, multimap<string,hit> > *best_hits, bool is_dna)
{
    Log_output::append_msg(" preselecting targets with Exonerate.",0);

    int r = rand();
    string tmp_dir = this->get_temp_dir();

    this->write_exonerate_input(target_sequences,reads,&r);


    string data_type = "-T protein -Q protein";
    if(is_dna)
        data_type = "-T dna -Q dna";

    // exonerate command for local alignment

    stringstream command;
    command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no "<<data_type<<" 2>&1";

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }

    Log_output::write_out("Exonerate_queries: command: "+command.str()+"\n",2);

    // read exonerate output, summing the multiple hit scores

    char line[256];
    map<string,multimap<string,hit> > all_hits;

    while ( fgets( line, sizeof line, fpipe))
    {
        this->read_output_line(&all_hits,line);
    }
    pclose(fpipe);


    this->find_hits_for_queries(&all_hits, reads, best_hits);

//    Log_output::write_out("Exonerate_reads: preselection finds "+Log_output::itos(all_hits.size())+" hits\n",2);


    for(vector<Fasta_entry>::iterator ri = reads->begin() ;ri != reads->end(); ri++ )
    {

        map<string,multimap<string,hit> >::iterator iter = best_hits->find(ri->name);
        if( iter != best_hits->end() )
        {
            multimap<string,hit>::iterator iter2 = iter->second.begin();
            for( ;iter2 != iter->second.end(); iter2++ )
            {
                if(selected_sequences->find(iter2->second.node)==selected_sequences->end())
                {
                    string name = iter2->second.node;
                    string seq = target_sequences->find(name)->second;
                    selected_sequences->insert(make_pair(name,seq));
                }
            }
        }
    }

}


void Exonerate_queries::read_output_line(map<string,multimap<string,hit> > *all_hits, string line)
{
    hit h;
    bool valid = this->split_sugar_string(string(line),&h);

    if(valid)
    {

        map<string,multimap<string,hit> >::iterator iter = all_hits->find(h.query);

        if( iter != all_hits->end() )
        {

            multimap<string,hit>::iterator iter2 = iter->second.find(h.node);

            if( iter2 != iter->second.end() )
            {
                if(iter2->second.t_strand == h.t_strand && iter2->second.q_strand == h.q_strand)
                {
                    iter2->second.score += h.score;

                    if(iter2->second.q_start > h.q_start)
                        iter2->second.q_start = h.q_start;
                    if(iter2->second.q_end < h.q_end)
                        iter2->second.q_end = h.q_end;
                    if(iter2->second.t_start > h.t_start)
                        iter2->second.t_start = h.t_start;
                    if(iter2->second.t_end < h.t_end)
                        iter2->second.t_end = h.t_end;
                }
                else if(iter2->second.score < h.score)
                {
                    iter2->second = h;
                }
            }
            else
            {
                iter->second.insert( make_pair(h.node, h) );
            }
        }
        else
        {
            multimap<string,hit> new_hit;
            new_hit.insert( make_pair(h.node, h) );
            all_hits->insert( make_pair(h.query, new_hit ) );
        }
    }
}


void Exonerate_queries::find_hits_for_queries(map<string,multimap<string,hit> > *all_hits, vector<Fasta_entry> *reads, map<string,multimap<string,hit> > *hits, bool is_local)
{
    vector<Fasta_entry>::iterator ri = reads->begin();

    for( ;ri != reads->end(); ri++ )
    {

        map<string,multimap<string,hit> >::iterator iter = all_hits->find(ri->name);
        if( iter != all_hits->end() )
        {

            vector<hit> hits_for_this;

            multimap<string,hit>::iterator iter2 = iter->second.begin();
            for( ;iter2 != iter->second.end(); iter2++ )
            {
                hits_for_this.push_back(iter2->second);
            }

            sort (hits_for_this.begin(), hits_for_this.end(), Exonerate_queries::better);

            multimap<string,hit> best_hits_for_this;

            if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                    Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

                ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                    Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
            {
                int lim = hits_for_this.at(0).score;
                if(is_local)
                    lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
                else
                    lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );

                for(int i=0; i<(int)hits_for_this.size(); i++)
                {
                    if(hits_for_this.at(i).score > lim)
                    {
                        best_hits_for_this.insert(pair<string,hit>(hits_for_this.at(i).node,hits_for_this.at(i)));

                        Log_output::write_out("Exonerate_reads: adding "+hits_for_this.at(i).node+" "+Log_output::itos(hits_for_this.at(i).score)+"\n",3);
                    }
                }
            }

            // keep a fixed number of hits

            else
            {
                int lim = hits_for_this.size();
                if( is_local )
                    lim = Settings::exonerate_local_keep_best;

                if( !is_local )
                    lim = Settings::exonerate_gapped_keep_best;

                for(int i=0; i<lim && i<(int)hits_for_this.size(); i++)
                {
                    best_hits_for_this.insert(pair<string,hit>(hits_for_this.at(i).node,hits_for_this.at(i)));

                    Log_output::write_out("Exonerate_reads: adding "+hits_for_this.at(i).node+" "+Log_output::itos(hits_for_this.at(i).score)+"\n",3);
                }
            }

            hits->insert( pair<string,multimap<string,hit> >(ri->name,best_hits_for_this) );
        }

        else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
        {
            ri->node_to_align = "discarded_read";
        }

    }

    if(Settings::noise>1)
    {
        ri = reads->begin();
        for( ;ri != reads->end(); ri++ )
        {

            map<string,multimap<string,hit> >::iterator iter = hits->find(ri->name);
            if( iter != hits->end() )
            {
                Log_output::write_out("Exonerate_reads: "+iter->first+" has "+Log_output::itos(iter->second.size())+" hits\n",2);

                multimap<string,hit>::iterator iter2 = iter->second.begin();

                for( ;iter2 != iter->second.end(); iter2++ )
                {
                    Log_output::write_out("  "+ri->name+" matches "+iter2->first+" with score "+Log_output::itos(iter2->second.score)+"\n",2);

                }
            }
        }
    }

}



void Exonerate_queries::local_alignment(map<string,string> *target_sequences, Fasta_entry *read, map<string,hit> *hits, bool is_local, bool is_dna)
{
    if(is_local)
        Log_output::append_msg(" running Exonerate with one query sequence (ungapped).",0);
    else
        Log_output::append_msg(" running Exonerate with one query sequence (gapped).",0);

    int r = rand();
    string tmp_dir = this->get_temp_dir();


    this->write_exonerate_input(target_sequences,read,&r);

    // exonerate command for local alignment

    string data_type = "-T protein -Q protein";
    if(is_dna)
        data_type = "-T dna -Q dna";

    stringstream command;
    if(is_local)
        command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no "<<data_type<<" 2>&1";
    else
        command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -m affine:local "<<data_type<<" 2>&1";

    Log_output::write_out("Exonerate_local: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }

    // read exonerate output, summing the multiple hit scores

    char line[256];
    map<string,hit> all_hits;
    vector<string> hit_names;

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {
            map<string,hit>::iterator iter = all_hits.find(h.node);
            if( iter != all_hits.end() )
            {
                if(iter->second.t_strand == h.t_strand && iter->second.q_strand == h.q_strand)
                {
                    iter->second.score += h.score;

                    if(iter->second.q_start > h.q_start)
                        iter->second.q_start = h.q_start;
                    if(iter->second.q_end < h.q_end)
                        iter->second.q_end = h.q_end;
                    if(iter->second.t_start > h.t_start)
                        iter->second.t_start = h.t_start;
                    if(iter->second.t_end < h.t_end)
                        iter->second.t_end = h.t_end;
                }
                else if(iter->second.score < h.score)
                {
                    iter->second = h;
                }
            }
            else
            {
                all_hits.insert( make_pair(h.node, h) );
                hit_names.push_back(h.node);
            }
        }
    }
    pclose(fpipe);


    Log_output::write_out("Exonerate_reads: "+read->name+" has "+Log_output::itos(hit_names.size())+" hits\n",2);

//    Log_output::write_out("Exonerate_reads: "+read->name+" has "+Log_output::itos(hit_names.size())+" hits\n",2);

    if(hit_names.size()>0)
    {
        hits->clear();

        vector<hit> best_hits;
        vector<string>::iterator iter = hit_names.begin();

        for(;iter!=hit_names.end();iter++)
        {
             map<string,hit>::iterator iter2 = all_hits.find(*iter);
             if( iter2 != all_hits.end() )
             {
               best_hits.push_back(iter2->second);
             }
        }

        sort (best_hits.begin(), best_hits.end(), Exonerate_queries::better);


        // keep hits that are above a relative threshold

        if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

            ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
        {
            int lim = best_hits.at(0).score;
            if(is_local)
                lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
            else
                lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );


            for(int i=0; i<(int)hit_names.size(); i++)
            {
                if(best_hits.at(i).score > lim)
                {
                    hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                    Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
                }
            }
        }


        // keep a fixed number of hits

        else
        {
            int lim = hit_names.size();
            if( is_local )
                lim = Settings::exonerate_local_keep_best;

            if( !is_local )
                lim = Settings::exonerate_gapped_keep_best;

            for(int i=0; i<lim && i<(int)hit_names.size(); i++)
            {
                hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
            }
        }
    }
    else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
    {
        read->node_to_align = "discarded_read";
    }

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);

}

void Exonerate_queries::local_alignment(Node *root, Fasta_entry *read, multimap<string,string> *tid_nodes, map<string,hit> *hits, bool is_local, bool is_dna, bool all_nodes)
{
    if(is_local)
        Log_output::append_msg(" running Exonerate with one query sequence (ungapped).",0);
    else
        Log_output::append_msg(" running Exonerate with one query sequence (gapped).",0);

    int r = rand();
    string tmp_dir = this->get_temp_dir();


    map<string,string> names;
    multimap<string,string>::iterator it = tid_nodes->begin();
    for(;it!=tid_nodes->end();it++)
    {
        if(all_nodes)
        {
            names.insert(pair<string,string>(it->second,"empty"));
        }
        else
        {
            if(it->first==read->tid)
            {
                names.insert(pair<string,string>(it->second,it->first));
            }
        }
    }

    this->write_exonerate_input(root,read,&names,&r);

    string data_type = "-T protein -Q protein";
    if(is_dna)
        data_type = "-T dna -Q dna";

   // exonerate command for local alignment

    stringstream command;
    if(is_local)
        command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no "<<data_type<<" 2>&1";
    else
        command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -m affine:local "<<data_type<<" 2>&1";

    Log_output::write_out("Exonerate_local: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }

    // read exonerate output, summing the multiple hit scores

    char line[256];
    map<string,hit> all_hits;
    vector<string> hit_names;

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {
            map<string,hit>::iterator iter = all_hits.find(h.node);
            if( iter != all_hits.end() )
            {
                if(iter->second.t_strand == h.t_strand && iter->second.q_strand == h.q_strand)
                {
                    iter->second.score += h.score;

                    if(iter->second.q_start > h.q_start)
                        iter->second.q_start = h.q_start;
                    if(iter->second.q_end < h.q_end)
                        iter->second.q_end = h.q_end;
                    if(iter->second.t_start > h.t_start)
                        iter->second.t_start = h.t_start;
                    if(iter->second.t_end < h.t_end)
                        iter->second.t_end = h.t_end;
                }
                else if(iter->second.score < h.score)
                {
                    iter->second = h;
                }
            }
            else
            {
                all_hits.insert( make_pair(h.node, h) );
                hit_names.push_back(h.node);
            }
        }
    }
    pclose(fpipe);


    Log_output::write_out("Exonerate_reads: "+read->name+" has "+Log_output::itos(hit_names.size())+" hits\n",2);

    if(hit_names.size()>0)
    {
        tid_nodes->clear();
        hits->clear();

        vector<hit> best_hits;
        vector<string>::iterator iter = hit_names.begin();

        for(;iter!=hit_names.end();iter++)
        {
             map<string,hit>::iterator iter2 = all_hits.find(*iter);
             if( iter2 != all_hits.end() )
             {
               best_hits.push_back(iter2->second);
             }
        }

        sort (best_hits.begin(), best_hits.end(), Exonerate_queries::better);


        // keep hits that are above a relative threshold

        if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

            ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
        {
            int lim = best_hits.at(0).score;
            if(is_local)
                lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
            else
                lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );


            for(int i=0; i<(int)hit_names.size(); i++)
            {
                if(best_hits.at(i).score > lim)
                {
                    string tid = names.find(best_hits.at(i).node)->second;
                    tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                    hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                    Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
                }
            }
        }


        // keep a fixed number of hits

        else
        {
            int lim = hit_names.size();
            if( is_local )
                lim = Settings::exonerate_local_keep_best;

            if( !is_local )
                lim = Settings::exonerate_gapped_keep_best;

            for(int i=0; i<lim && i<(int)hit_names.size(); i++)
            {
                string tid = names.find(best_hits.at(i).node)->second;
                tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
            }
        }
    }
    else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
    {
        tid_nodes->clear();
        read->node_to_align = "discarded_read";
    }

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);
}

/****************************************************************************/

void Exonerate_queries::local_pairwise_alignment(string *str1,string *str2,vector<Substring_hit> *hits,int *best_reverse_hit)
{
    int r = rand();
    string tmp_dir = this->get_temp_dir();

    this->write_exonerate_input(str1,str2,&r);

    // exonerate command for local alignment

    stringstream command;
    command <<exoneratepath << "exonerate -q "+tmp_dir+"q"<<r<<".fas -t "+tmp_dir+"t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";

    Log_output::write_out("Exonerate_pairwise: command: "+command.str()+"\n",2);


    FILE *fpipe;

//    #pragma omp critical
//    {
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }
//    }

    // read exonerate output, summing the multiple hit scores

    char line[256];
    vector<hit> best_hits;

    while ( fgets( line, sizeof line, fpipe))
    {
//        cout<<line;
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
            best_hits.push_back( h);
    }
    pclose(fpipe);


    if(best_hits.size()>0)
    {

        sort (best_hits.begin(), best_hits.end(), Exonerate_queries::better);


        // keep hits above a threshold
        int hit_score = -1;
        if( Settings_handle::st.is("exonerate-hit-score") )
            hit_score = Settings_handle::st.get("exonerate-hit-score").as<int>();

        int hit_length = Settings_handle::st.get("exonerate-hit-length").as<int>();

        for(int i=0; i<(int)best_hits.size(); i++)
        {
            if( best_hits.at(i).q_strand == '-' || best_hits.at(i).t_strand == '-' ||
                    best_hits.at(i).q_start > best_hits.at(i).q_end || best_hits.at(i).t_start > best_hits.at(i).t_end )
            {
                if(best_reverse_hit != 0)
                    *best_reverse_hit = best_hits.at(i).score;
                if(i==0)
                    Log_output::write_out("Warning: Anchoring suggests reverse strand hit!\n",0);
                continue;
            }

            if(hit_score>0 && best_hits.at(i).score > hit_score)
            {
                Substring_hit hit;
                hit.start_site_1 = best_hits.at(i).q_start;
                hit.start_site_2 = best_hits.at(i).t_start;
                hit.length = best_hits.at(i).q_end - best_hits.at(i).q_start;
                hit.score = best_hits.at(i).score;

                hits->push_back(hit);
//                cout<<"hit: "<<hit.start_site_1<<" "<<hit.start_site_2<<" "<<hit.length<<"\n";
                Log_output::write_out("Exonerate_queries: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
            }
            else if(hit_length>0 && best_hits.at(i).q_end-best_hits.at(i).q_start > hit_length)
            {
                Substring_hit hit;
                hit.start_site_1 = best_hits.at(i).q_start;
                hit.start_site_2 = best_hits.at(i).t_start;
                hit.length = best_hits.at(i).q_end - best_hits.at(i).q_start;
                hit.score = best_hits.at(i).score;

                hits->push_back(hit);
//                cout<<"hit: "<<hit.start_site_1<<" "<<hit.start_site_2<<" "<<hit.length<<"\n";
                Log_output::write_out("Exonerate_queries: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
            }

        }
    }

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);

}

/****************************************************************************/


void Exonerate_queries::delete_files(int r)
{
    string tmp_dir = this->get_temp_dir();

    stringstream q_name;
    q_name <<tmp_dir<<"q"<<r<<".fas";

    stringstream t_name;
    t_name <<tmp_dir<<"t"<<r<<".fas";

    if( remove( q_name.str().c_str() ) != 0 )
       Log_output::write_out( "Error deleting file", 1);
    if( remove( t_name.str().c_str() ) != 0 )
       Log_output::write_out( "Error deleting file", 1 );
}
