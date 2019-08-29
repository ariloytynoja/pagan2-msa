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

#include "utils/mafft_alignment.h"
#include "utils/log_output.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;
using namespace ppa;

Mafft_alignment::Mafft_alignment()
{
}

bool Mafft_alignment::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    mafftpath = epath;
    epath = epath+"sh.exe "+epath+"mafft -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());
    return WEXITSTATUS(status) == 1;

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

    mafftpath = epath;
    epath = epath+"mafft -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 1)
        return true;

    mafftpath = "";
    status = system("mafft -h >/dev/null 2>/dev/null");

    if(Settings_handle::st.is("docker"))
        return true;

    return WEXITSTATUS(status) == 1;

    #endif
}

void Mafft_alignment::align_sequences(vector<Fasta_entry> *sequences)
{
    ofstream m_output;
    string tmp_dir = this->get_temp_dir();

    int input_sequences = (int) sequences->size();

    stringstream m_name;

    int r = rand();
    while(true)
    {
        m_name <<tmp_dir<<"m"<<r<<".fas";
        ifstream m_file(m_name.str().c_str());

        if(!m_file)
        {
            m_output.open( m_name.str().c_str(), (ios::out) );
            break;
        }
        r = rand();
    }


    map<string,string> dna_seqs;
    bool has_dna_seqs = (sequences->at(0).sequence.length() > 0);

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        string seq = si->sequence;
        if(seq.find('.') != string::npos)
            seq.erase( remove( seq.begin(), seq.end(), '.' ), seq.end() ) ;

        if(seq.find('-') != string::npos)
            seq.erase( remove( seq.begin(), seq.end(), '-' ), seq.end() ) ;

        m_output<<">"<<si->name<<endl<<seq<<endl;
        if(has_dna_seqs)
            dna_seqs.insert(pair<string,string>(si->name,si->dna_sequence));
    }
    m_output.close();
    sequences->clear();


    stringstream command;
    command << mafftpath<<"mafft "<<m_name.str()<<" 2>/dev/null";

    Log_output::write_out("Mafft: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with mafft pipe.\nExiting.\n",0);
        exit(1);
    }

    // read mafft output
    string name, sequence = "";  // Initialization
    char temp[256];

    while ( fgets( temp, sizeof temp, fpipe))
    {
    	string line(temp);

        if (line[0] == '>')
        {
            line = this->remove_last_whitespaces(line);

            // If a name and a sequence were found
            if ((name != "") && (sequence != ""))
            {
                Fasta_entry s;
                s.name = name;
                sequence = this->remove_whitespaces(sequence);
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
                s.sequence = sequence;

                if(has_dna_seqs)
                {
                    map<string,string>::iterator it = dna_seqs.find(name);
                    if(it!=dna_seqs.end())
                        s.dna_sequence = it->second;
                }

                sequences->push_back(s);
                name = "";
                sequence = "";
            }
            name = line;
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            sequence += temp;  // Sequence isolation
        }
    }

    // Addition of the last sequence in file
    if ((name != "") && (sequence != ""))
    {
        Fasta_entry s;
        s.name = name;
        sequence = this->remove_whitespaces(sequence);
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
        s.sequence = sequence;

        if(has_dna_seqs)
        {
            map<string,string>::iterator it = dna_seqs.find(name);
            if(it!=dna_seqs.end())
                s.dna_sequence = it->second;
        }

        sequences->push_back(s);
    }

    pclose(fpipe);

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);

    if( input_sequences != (int) sequences->size() )
    {
        Log_output::write_out("Problems with mafft.\nExiting.\n",0);
        exit(1);
    }


}

void Mafft_alignment::align_sequences_fifo(vector<Fasta_entry> *sequences)
{
//    ofstream m_output;
    string tmp_dir = this->get_temp_dir();

    stringstream m_name;

    int r = rand();
    while(true)
    {
//        stringstream m_name;
        m_name <<tmp_dir<<"m"<<r<<".fas";
        ifstream m_file(m_name.str().c_str());

        if(!m_file)
        {
//            m_output.open( m_name.str().c_str(), (ios::out) );
            break;
        }
        r = rand();
    }


    map<string,string> dna_seqs;
    bool has_dna_seqs = (sequences->at(0).sequence.length() > 0);

//    vector<Fasta_entry>::iterator si = sequences->begin();
//    for(;si!=sequences->end();si++)
//    {
//        m_output<<">"<<si->name<<endl<<si->sequence<<endl;
//        if(has_dna_seqs)
//            dna_seqs.insert(pair<string,string>(si->name,si->dna_sequence));
//    }
//    m_output.close();
//    sequences->clear();


    stringstream command;
    command << mafftpath<<"mafft "+tmp_dir+"m"<<r<<".fas  2>/dev/null";

    Log_output::write_out("Mafft: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with mafft pipe.\nExiting.\n",0);
        exit(1);
    }

    int rv = mkfifo (m_name.str().c_str(), S_IRUSR| S_IWUSR);
    if (rv < 0) {
        Log_output::write_out("Problems with mafft fifo.\nExiting.\n",0);
        exit(1);
    }

    FILE *ffifo;

    if((ffifo = fopen(m_name.str().c_str(), "w")) == NULL)
    {
        Log_output::write_out("Problems with mafft fifo.\nExiting.\n",0);
        exit(1);
    }

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        fprintf(ffifo,">%s\n%s\n",si->name.c_str(),si->sequence.c_str());

        if(has_dna_seqs)
            dna_seqs.insert(pair<string,string>(si->name,si->dna_sequence));
    }
    fclose(ffifo);
    sequences->clear();
    //

    // read mafft output
    string name, sequence = "";  // Initialization
    char temp[256];

    while ( fgets( temp, sizeof temp, fpipe))
    {
        string line(temp);

        if (line[0] == '>')
        {
            line = this->remove_last_whitespaces(line);

            // If a name and a sequence were found
            if ((name != "") && (sequence != ""))
            {
                Fasta_entry s;
                s.name = name;
                sequence = this->remove_whitespaces(sequence);
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
                s.sequence = sequence;

                if(has_dna_seqs)
                {
                    map<string,string>::iterator it = dna_seqs.find(name);
                    if(it!=dna_seqs.end())
                        s.dna_sequence = it->second;
                }

                sequences->push_back(s);
                name = "";
                sequence = "";
            }
            name = line;
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            sequence += temp;  // Sequence isolation
        }
    }

    // Addition of the last sequence in file
    if ((name != "") && (sequence != ""))
    {
        Fasta_entry s;
        s.name = name;
        sequence = this->remove_whitespaces(sequence);
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
        s.sequence = sequence;

        if(has_dna_seqs)
        {
            map<string,string>::iterator it = dna_seqs.find(name);
            if(it!=dna_seqs.end())
                s.dna_sequence = it->second;
        }

        sequences->push_back(s);
    }

    pclose(fpipe);

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);
}

void Mafft_alignment::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream m_name;
    m_name <<tmp_dir<<"m"<<r<<".fas";

    if ( remove( m_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);
}
