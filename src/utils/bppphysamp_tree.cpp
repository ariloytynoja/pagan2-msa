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

#include "utils/bppphysamp_tree.h"
#include "utils/log_output.h"
#include "utils/fasta_reader.h"
#include "utils/newick_reader.h"
#include "utils/tree_node.h"

using namespace std;
using namespace ppa;

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

BppPhySamp_tree::BppPhySamp_tree()
{
}

bool BppPhySamp_tree::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);
    path[length] = '\0';

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"bppphysamp >/dev/null 2>/dev/null";
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
    epath = epath+"bppphysamp >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 0)
        return true;

    bppdistpath = "";
    status = system("bppphysamp >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 0;

    #endif
}

void BppPhySamp_tree::reduce_sequences(set<string> *toremove,bool is_protein)
{
    string seqfile = "";
    string treefile = "";

    if(Settings_handle::st.is("ref-seqfile"))
        seqfile =  Settings_handle::st.get("ref-seqfile").as<string>();
    if(Settings_handle::st.is("ref-treefile"))
        treefile =  Settings_handle::st.get("ref-treefile").as<string>();

    if(seqfile != "" && treefile != "")
    {
        // Copy fasta file
        //
        string tmpseqfile = seqfile+".clean";
        try
        {
            Fasta_reader fr;
            vector<Fasta_entry> sequences;
            fr.read(seqfile, sequences, true, false);

            int datatype = fr.check_sequence_data_type(&sequences);
            fr.check_alphabet(&sequences,datatype);

            for(int i=0;i<(int)sequences.size();i++)
                sequences.at(i).comment="";
            fr.write(tmpseqfile,sequences,"fasta");
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error writing a temporary sequence file for BppPhySamp: '"+seqfile+"' or '"+tmpseqfile+"' fails.\nExiting.\n\n",0);
            exit(1);
        }
        tmpseqfile.append(".fas");

        // Copy newick file
        //
        Newick_reader nr;
        string tree;
        Node *tmp_root;
        try
        {
            tree = nr.read_tree(treefile);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the tree file '"+treefile+"'.\n",0);
            return;
        }
        try {
            tmp_root = nr.parenthesis_to_tree(tree);
        }
        catch(exception e)
        {
            Tree_node tn;
            tree = tn.get_rooted_tree(tree);
            tmp_root = nr.parenthesis_to_tree(tree);
        }
        string tmptreefile = treefile+".clean";
        tmp_root->write_nhx_tree(tmptreefile,"tre");
        tmptreefile.append(".tre");
        //

        stringstream command;
        command << bppdistpath<<"bppphysamp input.tree.file="<<tmptreefile<<" input.sequence.file="<<tmpseqfile<<" input.method=tree "
                   << " input.sequence.format=Fasta input.tree.format=Newick choice_criterion=length output.sequence.format=Fasta ";

        if(is_protein)
            command << "alphabet=Protein";
        else
            command << "alphabet=DNA";

        if( Settings_handle::st.is("prune-keep-threshold") )
        {
            float threshold = Settings_handle::st.get("prune-keep-threshold").as<float>();
            command << " threshold="<<threshold<<" deletion_method=threshold";
        }
        else
        {
            int number = Settings_handle::st.get("prune-keep-number").as<int>();
            command <<  " sample_size="<<number<<" deletion_method=sample";
        }

        Log_output::write_out("BppPhySamp_tree: command: "+command.str()+"\n",2);

//        cout<<endl<<command.str()<<endl;
        FILE *fpipe;
        if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
        {
            Log_output::write_out("Problems with bppphysamp pipe.\n",0);
            return;
        }

        char line[256];

        while ( fgets( line, sizeof line, fpipe))
        {
            Log_output::write_out("BppPhySamp: "+string(line),2);

            string linestr = string(line);
            if(linestr.find("Remove sequence") != string::npos)
            {
                linestr = linestr.substr(linestr.find(":")+2);
                linestr = linestr.substr(0,linestr.find("\n"));
                toremove->insert(linestr);
            }
        }
        pclose(fpipe);

        if( remove( tmpseqfile.c_str() ) != 0 )
           Log_output::write_out( "Error deleting file", 1);
        if( remove( tmptreefile.c_str() ) != 0 )
           Log_output::write_out( "Error deleting file", 1);
    }
}

// bppphysamp input.tree.file=reference_tree2.nhx input.sequence.file=reference_codon.fas input.method=tree sample_size=12 output.sequence.file=samples.fas alphabet=DNA deletion_method=sample input.sequence.format=Fasta input.tree.format=Newick choice_criterion=length output.sequence.format=Fasta

// bppphysamp input.tree.file=reference_tree2.nhx input.sequence.file=reference_codon.fas input.method=tree output.sequence.file=samples.fas alphabet=DNA deletion_method=threshold threshold=0.3 input.sequence.format=Fasta input.tree.format=Newick choice_criterion=length output.sequence.format=Fasta


void BppPhySamp_tree::delete_files(int r)
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
