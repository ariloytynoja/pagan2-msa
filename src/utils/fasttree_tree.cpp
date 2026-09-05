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

#include "utils/fasttree_tree.h"
#include "utils/helper_probe.h"
#include "utils/model_factory.h"
#include <netdb.h>

using namespace std;
using namespace ppa;

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif


// FastTree given no input file reads the ALIGNMENT FROM STANDARD INPUT -- it
// announces this itself, "Alignment: standard input". So a bare `fasttree`
// never exits 0, and both ways it can end are wrong for a presence probe:
//
//   stdin closed  ->  reads EOF, has no alignment to build a tree from,
//                     exits 1
//   stdin open    ->  BLOCKS FOREVER waiting for one
//
// The second case is not hypothetical. Any parent that did not redirect its
// own stdin -- an interactive shell, a pipeline, a job launcher, a test
// harness -- hangs pagan2 inside this probe, before it has read a single
// sequence. The `</dev/null` in the probes below is therefore load-bearing,
// not tidiness.
//
// The first case is why the probe used to answer "FastTree is not available"
// on machines where FastTree is installed and works: it demanded exit 0, a
// status this invocation cannot produce. The previous code compensated with a
// hard-coded `gethostname() == "wasabi2"` special case that accepted exit 1 on
// that one host -- the same observation, confined to one machine.
//
// helper_was_found() replaces the exit-status comparison; see helper_probe.h
// for why the shell's 127/126 is the part that does not depend on FastTree.

// The FastTree executable to probe for and to run.
//
// Upstream distributes FastTree, FastTreeMP and FastTreeUPGMA; none of them is
// called "fasttree", so the hard-coded name below finds a binary only where a
// distribution has added a lowercase name or symlink. Naming it is also the
// only way to CHOOSE: FastTreeMP is the OpenMP build, and FastTreeUPGMA
// computes UPGMA rather than approximately-ML trees, which is a modelling
// decision and not something to guess at by sniffing the filesystem.
//
// The default keeps today's behaviour exactly.
static string fasttree_exec()
{
    if(Settings_handle::st.is("fasttree-exec"))
        return Settings_handle::st.get("fasttree-exec").as<string>();

    return "fasttree";
}

FastTree_tree::FastTree_tree()
{
}

bool FastTree_tree::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    // readlink() returns -1 on failure (no /proc, or not readable) and never
    // NUL-terminates; it also silently truncates a path that does not fit.
    // path[-1] would be a write outside the buffer, so treat both as "no
    // directory known here" and fall through to looking the program up on PATH.
    if(length < 0 || length >= 200-1)
        length = 0;
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    progpath = epath;
    epath = epath+fasttree_exec()+" </dev/null >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return helper_was_found(status);

    # else

    char path[200];
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    // readlink() returns -1 on failure (no /proc, or not readable) and never
    // NUL-terminates; it also silently truncates a path that does not fit.
    // path[-1] would be a write outside the buffer, so treat both as "no
    // directory known here" and fall through to looking the program up on PATH.
    if(length < 0 || length >= 200-1)
        length = 0;
    path[length] = '\0';
    epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #endif

    progpath = epath;
    epath = epath+fasttree_exec()+" </dev/null >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(helper_was_found(status))
        return true;

    progpath = "";
    status = system((fasttree_exec()+" </dev/null >/dev/null 2>/dev/null").c_str());

    if(Settings_handle::st.is("docker"))
        return true;

    return helper_was_found(status);

    #endif
}

string FastTree_tree::infer_phylogeny(std::vector<Fasta_entry> *sequences,bool is_protein,int n_threads)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    stringstream t_name;

    int r = rand();
    while(true)
    {

        f_name <<tmp_dir<<"d"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        /*
        t_name <<tmp_dir<<"d"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());
        */
        if(!f_file /*&& !t_file*/)
        {
            break;
        }
        r = rand();
    }

    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        f_output<<">"<<si->name<<"\n"<<si->sequence<<"\n";
    }
    f_output.close();

    stringstream command;
    command << progpath<<fasttree_exec()<<" -quiet -nopr -nosupport ";
    if(is_protein)
        command << f_name.str() << " 2>/dev/null";
    else
        command << "-nt "<<f_name.str() << " 2>/dev/null";

    Log_output::write_out("FastTree_tree: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with fasttree pipe.\nExiting.\n",0);
        exit(1);
    }

    char line[256];
     string tree = "";
    while ( fgets( line, sizeof line, fpipe))
    {
        Log_output::write_out("FastTree: "+string(line),2);
        tree += line;
    }
    pclose(fpipe);

    Log_output::write_out("FastTree_tree: "+tree+"\n",2);

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);

    return tree;
}


void FastTree_tree::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    f_name <<tmp_dir<<"d"<<r<<".fas";

    if ( remove( f_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);

    /*
    stringstream t_name;
    t_name <<tmp_dir<<"d"<<r<<".tre";

    if ( remove( t_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);
    */
}
