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

#include "utils/raxml_tree.h"
#include "utils/model_factory.h"

using namespace std;
using namespace ppa;

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

RAxML_tree::RAxML_tree()
{
}

bool RAxML_tree::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));
    raxmlpath = epath;
    epath = epath+"raxml -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    # else
    int status = system("raxml -h >/dev/null 2>/dev/null");

    if(WEXITSTATUS(status) == 0)
        return true;

    char path[200];
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));
    //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';
    epath = string(path).substr(0,length);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));

    #endif

    raxmlpath = epath;
    epath = epath+"raxml -h >/dev/null 2>/dev/null";
    status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    #endif
}


string RAxML_tree::infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein,int n_threads)
{

    ofstream m_output;
    string tmp_dir = this->get_temp_dir();

    int r = rand();
    while(true)
    {
        stringstream m_name;
        m_name <<tmp_dir<<"RAxML_input.r"<<r;
        ifstream m_file(m_name.str().c_str());

        stringstream r_name;
        r_name <<tmp_dir<<"RAxML_info.r"<<r;
        ifstream r_file(r_name.str().c_str());

        if(!m_file && !r_file)
        {
            m_output.open( m_name.str().c_str(), (ios::out) );
            break;
        }
        r = rand();
    }
    vector<Fasta_entry>::iterator si = sequences->begin();
    m_output <<sequences->size()<<" "<<sequences->at(0).sequence.length()<<endl;
    for(;si!=sequences->end();si++)
    {
        m_output<<si->name<<endl<<si->sequence<<endl;
    }
    m_output.close();

    if(n_threads == 1)
        n_threads = 2;

    stringstream command;
    if(is_protein)
        command << raxmlpath<<"raxml -s "+tmp_dir+"RAxML_input.r"<<r<<" -c 4 -f d -m PROTCATJTT -w "+tmp_dir+" -n r"<<r<<" -p "<<r<<" -T "<<n_threads<<" 2>/dev/null";
    else
        command << raxmlpath<<"raxml -s "+tmp_dir+"RAxML_input.r"<<r<<" -c 4 -f d -m GTRCAT -w "+tmp_dir+" -n r"<<r<<" -p "<<r<<" -T "<<n_threads<<" 2>/dev/null";

    Log_output::write_out("Raxml: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with raxml pipe.\nExiting.\n",0);
        exit(1);
    }


    char line[256];

    while ( fgets( line, sizeof line, fpipe))
    {
        Log_output::write_out("RAxML: "+string(line),2);
    }
    pclose(fpipe);

    stringstream t_name;
    t_name <<tmp_dir<<"RAxML_bestTree.r"<<r;

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


void RAxML_tree::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    vector<string> f_names;
    f_names.push_back("RAxML_input.r");
    f_names.push_back("RAxML_bestTree.r");
    f_names.push_back("RAxML_info.r");
    f_names.push_back("RAxML_log.r");
    f_names.push_back("RAxML_parsimonyTree.r");
    f_names.push_back("RAxML_result.r");

    for(int i=0;i<(int)f_names.size();i++)
    {
        stringstream m_name;
        m_name <<tmp_dir<<f_names.at(i)<<r;

        if ( remove( m_name.str().c_str() ) != 0 )
            Log_output::write_out( "Error deleting file", 1);
    }
}
