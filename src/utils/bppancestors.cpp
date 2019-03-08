/***************************************************************************
 *   Copyright (C) 2013 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
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

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>
#include "utils/bppancestors.h"
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "utils/settings_handle.h"
#include "utils/model_factory.h"
#include "utils/codon_translation.h"

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;
using namespace ppa;

BppAncestors::BppAncestors() {}

bool BppAncestors::test_executable()
{

    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    # else

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

    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 0)
        return true;

    bppdistpath = "";
    status = system("bppancestor >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 0;

    #endif
}

bool BppAncestors::infer_ancestors(Node *root,vector<Fasta_entry> *aligned_sequences,bool isCodon)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    stringstream t_name;
    stringstream o_name;

    int r = rand();
    while(true)
    {

        f_name.str("");
        t_name.str("");
        o_name.str("");

        f_name <<tmp_dir<<"f"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"t"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        o_name <<tmp_dir<<"o"<<r<<".fas";
        ifstream o_file(t_name.str().c_str());

        if(!f_file && !t_file && !o_file)
        {
            ofstream f_tmp;
            f_tmp.open(f_name.str().c_str(), (ios::out) );
            ofstream t_tmp;
            t_tmp.open(t_name.str().c_str(), (ios::out) );
            ofstream o_tmp;
            o_tmp.open(o_name.str().c_str(), (ios::out) );

            break;
        }
        r = rand();
    }


    ////////////


    bool isDna = root->get_sequence()->get_data_type() == Model_factory::dna;


    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );

    int count = 0;
    map<string,string> tmp_names;
    vector<Fasta_entry>::iterator si = aligned_sequences->begin();
    for(;si!=aligned_sequences->end();si++)
    {
        stringstream name;
        name<<"seq"<<count++;
        tmp_names.insert(make_pair(si->name,name.str()));

        if(isCodon)
        {
            string seq = si->sequence;
            string codon = seq.substr(seq.length()-3,3);
            if(codon == "TGA" || codon == "TAG" || codon == "TAA")
                seq.replace(seq.length()-3,3,"NNN");
            f_output<<">"<<name.str()<<"\n"<<seq<<"\n";
        }
        else
            f_output<<">"<<name.str()<<"\n"<<si->sequence<<"\n";
    }
    f_output.close();


    int add = root->get_number_of_leaves();
    string tree = root->print_bppa_tree(&add);

    for(map<string,string>::iterator it = tmp_names.begin();it != tmp_names.end(); it++)
    {
        size_t pos = 0;
        if((pos = tree.find(it->first)) != std::string::npos)
            tree.replace(pos, it->first.length(), it->second);
    }

    ofstream t_output;
    t_output.open( t_name.str().c_str(), (ios::out) );
    t_output<<tree<<endl;
    t_output.close();


    ////////////


    stringstream command;
    command << bppdistpath<<"bppancestor input.sequence.file="<<f_name.str()<<" input.sequence.format=Fasta input.sequence.sites_to_use=all input.tree.file="<<t_name.str()<<
            " input.tree.format=NHX input.sequence.max_gap_allowed=100% initFreqs=observed output.sequence.file="<<o_name.str()<<" output.sequence.format=Phylip";

    if(isCodon)
    {
        command << " alphabet=Codon\\(letter=DNA,type=Standard\\) model=YN98\\(kappa=2,omega=0.5\\)";
    }
    else
    {
        if(!isDna)
            command << " alphabet=Protein model=WAG01";
        else
        {
            if(Settings_handle::st.is("codons") )
                command << " alphabet=Codon\\(letter=DNA,type=Standard\\) model=YN98\\(kappa=2,omega=0.5\\)";
            else
                command << " alphabet=DNA model=HKY85";
        }
    }

    command << " 2>&1";

    Log_output::write_out("BppAncestors: command: "+command.str()+"\n",2);


    try {

        FILE *fpipe;
        char line[256];

        if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
        {
            Log_output::write_out("Problems with bppancestor pipe.\nExiting.\n",0);
            exit(1);
        }
        while ( fgets( line, sizeof line, fpipe))
        {
            Log_output::write_out("BppAncestors: "+string(line),3);
        }
        pclose(fpipe);


        ////////////


        map<string,string> bppa_sequences;

        Fasta_reader fr;

        fr.read_bpp_phylip(o_name.str().c_str(),&bppa_sequences);

        si = aligned_sequences->begin();
        for(;si!=aligned_sequences->end();si++)
        {
            map<string,string>::iterator fi = bppa_sequences.find(si->name);
            if(fi != bppa_sequences.end())
            {
                for(int i=0;i<(int)si->sequence.length();i++)
                {
                    if( si->sequence.at(i)!='-' && si->sequence.at(i)!='.' )
                        si->sequence.at(i) = fi->second.at(i);
                }
            }
        }

    } catch(Exception e)
    {
        Log_output::write_out("\nReconstructing ML ancestral sequences failed. Outputting parsimony ancestors.\n\n",0);

        if(!Settings_handle::st.is("keep-temp-files"))
            this->delete_files(r);

        return false;
    }

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);


    ////////////


    bool copy_left =
        root->get_left_child()->get_distance_to_parent() < root->get_right_child()->get_distance_to_parent();

    string left_seq,right_seq,root_seq;

    si = aligned_sequences->begin();
    for(;si!=aligned_sequences->end();si++)
    {
        if(si->name == root->get_left_child()->get_name())
            left_seq = si->sequence;
        if(si->name == root->get_right_child()->get_name())
            right_seq = si->sequence;
        if(si->name == root->get_name())
            root_seq = si->sequence;
    }

    for(int i=0;i<(int)root_seq.length();i++)
    {
        bool left_gap  = left_seq.at(i)=='-'  || left_seq.at(i)=='.';
        bool right_gap = right_seq.at(i)=='-' || right_seq.at(i)=='.';

        if(left_gap && right_gap)
            ;
        else if(left_gap && !right_gap)
            root_seq.at(i) = right_seq.at(i);
        else if(!left_gap && right_gap)
            root_seq.at(i) = left_seq.at(i);
        else if(!left_gap && !right_gap)
        {
            if(copy_left)
                root_seq.at(i) = left_seq.at(i);
            else
                root_seq.at(i) = right_seq.at(i);
        }
    }

    si = aligned_sequences->begin();
    for(;si!=aligned_sequences->end();si++)
    {
        if(si->name == root->get_name())
            si->sequence = root_seq;
    }

    return true;
}

void BppAncestors::count_events(Node *root,vector<Fasta_entry> *aligned_sequences,string outfile,bool isCodon)
{

    bool isDna = root->get_sequence()->get_data_type() == Model_factory::dna;

    if( isDna && Settings_handle::st.is("codons") )
        isCodon = true;

    // make map of sequences
    map<string,string> sequences;
    for(int j=0;j<(int)aligned_sequences->size();j++)
    {
        sequences.insert( sequences.begin(),pair<string,string>(aligned_sequences->at(j).name,aligned_sequences->at(j).sequence) );
    }


    // list states
    map<string,int> alphabet;

    int wordsize = 1;
    if( isCodon )
    {
        std::vector<std::string> *a = Model_factory::get_codon_full_character_alphabet();

        for(int i=0;i<(int)a->size();i++)
            alphabet.insert(alphabet.begin(),pair<string,int>(a->at(i),i));

        alphabet.insert(alphabet.begin(),pair<string,int>("---",-1));
        alphabet.insert(alphabet.begin(),pair<string,int>("...",-1));
        wordsize = 3;
    }
    else
    {
        std::vector<std::string> *a;
        if(isDna)
            a = Model_factory::get_dna_full_character_alphabet();
        else
            a = Model_factory::get_protein_full_character_alphabet();

        for(int i=0;i<(int)a->size();i++)
            alphabet.insert(alphabet.begin(),pair<string,int>(a->at(i),i));

        alphabet.insert(alphabet.begin(),pair<string,int>("-",-1));
        alphabet.insert(alphabet.begin(),pair<string,int>(".",-1));
    }


    // list parent-child pairs
    vector<pair<string,string> > pairs;
    root->get_parent_child_pairs(&pairs);

    Codon_translation ct;
    ct.define_translation_tables();


    ofstream output;
    if(Settings_handle::st.is("events"))
    {
        output.open( (outfile+".events").c_str() );

        output<<"Alignment topology with node labels:"<<endl<<endl;
        output<<root->print_tagged_topology()<<endl<<endl;

        output<<"Inferred evolutionary events per branch:"<<endl;
    }

    // go through all pairs, count/list events
    for(int j=0;j<(int)pairs.size();j++)
    {
        string p = pairs.at(j).first;
        string c = pairs.at(j).second;

        string ps = sequences.find(p)->second;
        string cs = sequences.find(c)->second;

        stringstream subs;
        stringstream ins;
        stringstream dels;

        bool pg = false;
        bool cg = false;

        if(Settings_handle::st.is("events"))
            output<<endl<<"branch "<<c<<endl;

        for(int i=0;i<(int)ps.length();i+=wordsize)
        {
            int site = i+1;
            if(isCodon)
                site = i/3+1;

            string pc = ps.substr(i,wordsize);
            string cc = cs.substr(i,wordsize);

            int pi = alphabet.find(pc)->second;
            int ci = alphabet.find(cc)->second;

            if(pi>=0 && ci<0)
            {
                if(not cg)
                {
                    cg = true;
                    if(isCodon)
                        dels << " "<<i/3+1;
                    else
                        dels << " "<<i+1;
                }
            }
            if(pi<0 && ci>=0)
            {
                if(not pg)
                {
                    pg = true;
                    if(isCodon)
                        ins << " "<<i/3+1;
                    else
                        ins << " "<<i+1;
                }
            }
            if(pi>=0 && pg)
            {
                pg = false;
                ins << ".."<<site-1<<" insertion"<<endl;
            }
            if(ci>=0 && cg)
            {
                cg = false;
                dels << ".."<<site-1<<" deletion"<<endl;
            }
            if(pi>=0 && ci>=0 && pi != ci)
            {
                subs<<" "<<site<<" "<<pc<<" -> "<<cc;
                if(isCodon)
                {
                    if(ct.codon_to_amino(pc) == ct.codon_to_amino(cc))
                        subs<<" ("<<ct.codon_to_amino(pc)<<")";
                    else
                        subs<<" ("<<ct.codon_to_amino(pc) << " -> " << ct.codon_to_amino(cc)<<")";
                }
                subs<<endl;
            }

        }
        if(Settings_handle::st.is("events"))
        {
            output<<subs.str();
            output<<ins.str();
            output<<dels.str();
        }
    }

}

void BppAncestors::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream t_name;
    t_name <<tmp_dir<<"t"<<r<<".tre";

    stringstream f_name;
    f_name <<tmp_dir<<"f"<<r<<".fas";

    stringstream o_name;
    o_name <<tmp_dir<<"o"<<r<<".fas";

    if ( remove( t_name.str().c_str() ) != 0 )
        perror( "Error deleting file" );
    if ( remove( f_name.str().c_str() ) != 0 )
        perror( "Error deleting file");
    if ( remove( o_name.str().c_str() ) != 0 )
        perror( "Error deleting file");

}
