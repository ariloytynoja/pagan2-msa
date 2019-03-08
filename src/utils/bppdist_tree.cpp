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

#include "utils/bppdist_tree.h"
#include "utils/model_factory.h"

using namespace std;
using namespace ppa;

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

BppDist_tree::BppDist_tree()
{
}

bool BppDist_tree::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"bppdist >/dev/null 2>/dev/null";
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
    epath = epath+"bppdist >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 0)
        return true;

    bppdistpath = "";
    status = system("bppdist >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 0;

    #endif
}


string BppDist_tree::infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein,int n_threads)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    stringstream t_name;

    int r = rand();
    while(true)
    {

        f_name <<tmp_dir<<"d"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"d"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        if(!f_file && !t_file)
        {
            break;
        }
        r = rand();
    }


    stringstream command;
    command << bppdistpath<<"bppdist method=bionj output.tree.file="<<t_name.str()<<" output.matrix.file=none input.sequence.file="<<f_name.str()<<" input.sequence.format=Fasta "
               "input.sequence.sites_to_use=all optimization.method=init optimization.verbose=0 input.sequence.max_gap_allowed=100 initFreqs=observed ";
    if(is_protein)
        command << "alphabet=Protein model=WAG01";
    else
        command << "alphabet=DNA model=HKY85";

//    command << " > bppdist.log";

    Log_output::write_out("BppDist_tree: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with bppdist pipe.\nExiting.\n",0);
        exit(1);
    }

    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        f_output<<">"<<si->name<<"\n"<<si->sequence<<"\n";
    }
    f_output.close();

    char line[256];

    while ( fgets( line, sizeof line, fpipe))
    {
        Log_output::write_out("BppDist: "+string(line),2);
    }
    pclose(fpipe);

    ifstream input(t_name.str().c_str(), ios::in);
    string temp,tree = "";

    while(!input.eof())
    {
        getline(input, temp, '\n');
        tree += temp;
    }
    input.close();


    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);

    return tree;
}


string BppDist_tree::infer_phylogeny_fifo(std::vector<Fasta_entry> *sequences,bool is_protein,int n_threads)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    stringstream t_name;

    int r = rand();
    while(true)
    {
        f_name <<tmp_dir<<"d"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"d"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        if(!f_file && !t_file)
        {
            break;
        }
        r = rand();
    }


    stringstream command;
    command << bppdistpath<<"bppdist method=bionj output.tree.file="<<t_name.str()<<" output.matrix.file=none input.sequence.file="<<f_name.str()<<" input.sequence.format=Fasta "
               "input.sequence.sites_to_use=all optimization.method=init optimization.verbose=0 input.sequence.max_gap_allowed=100 initFreqs=observed ";
    if(is_protein)
        command << "alphabet=Protein model=WAG01";
    else
        command << "alphabet=DNA model=HKY85";

//    command << " > bppdist.log";

    Log_output::write_out("BppDist_tree: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with bppdist pipe.\nExiting.\n",0);
        exit(1);
    }

    int rv = mkfifo (f_name.str().c_str(), S_IRUSR| S_IWUSR);
    if (rv < 0) {
        Log_output::write_out("Problems with biodist fifo.\nExiting.\n",0);
        exit(1);
    }

    rv = mkfifo (t_name.str().c_str(), S_IRUSR| S_IWUSR);
    if (rv < 0) {
        Log_output::write_out("Problems with biodist fifo.\nExiting.\n",0);
        exit(1);
    }

    FILE *ffifo;
    if((ffifo = fopen(f_name.str().c_str(), "w")) == NULL)
    {
        Log_output::write_out("Problems with biodist fifo.\nExiting.\n",0);
        exit(1);
    }

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        fprintf(ffifo,">%s\n%s\n",si->name.c_str(),si->sequence.c_str());
    }


    fclose(ffifo);

    cout<<"fasta closed\n";

    FILE *tfifo;
    if((tfifo = fopen(t_name.str().c_str(), "r")) == NULL)
    {
        Log_output::write_out("Problems with biodist fifo.\nExiting.\n",0);
        exit(1);
    }

    cout<<"tree opened\n";

    char line[5000];
    string tree = "";

    while ( fgets( line, sizeof line, tfifo))
    {
        cout<<line<<endl;
        tree += string(line);
    }

//    ifstream input(t_name.str().c_str(), ios::in);
//    string temp,tree = "";

//    while(!input.eof())
//    {
//        getline(input, temp, '\n');
//        tree += temp;
//    }
//    input.close();


    fclose(tfifo);

    cout<<"tree closed\n";

    while ( fgets( line, sizeof line, fpipe))
    {
        Log_output::write_out("BppDist: "+string(line),2);
    }
    pclose(fpipe);


    this->delete_files(r);

    return tree;
}

void BppDist_tree::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    f_name <<tmp_dir<<"d"<<r<<".fas";

    if ( remove( f_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);

    stringstream t_name;
    t_name <<tmp_dir<<"d"<<r<<".tre";

    if ( remove( t_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);

}
