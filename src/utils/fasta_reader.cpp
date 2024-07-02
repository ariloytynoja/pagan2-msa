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

//
// This file is partly based on the code from the Bio++ library.
// The following copyright information is given.
//

//
// File: Fasta.cpp
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Aug 21 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

Julien.Dutheil@univ-montp2.fr

This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/


#include <algorithm>
#include <iostream>
#include <set>
#include "utils/fasta_reader.h"
#include "utils/text_utils.h"
#include "utils/settings_handle.h"

using namespace ppa;

/****************************************************************************************/

void Fasta_reader::read(istream & input, vector<Fasta_entry> & seqs, bool short_names, bool degap) const 
{
    if (!input) { throw IOException ("Fasta_reader::read. Failed to open file"); }


    char c = input.get();
    while(c==' ' || c=='\n')
    {
        c = input.get();
    }

    if(c=='>')
    {
        input.unget();
        this->read_fasta(input,seqs,short_names,degap);
    }
    else if(c=='@')
    {
        input.unget();
        this->read_fastq(input,seqs);
    }
    else if(c=='#')
    {
        input.unget();
        this->read_graph(input,seqs,short_names);
    }
    else
    {
        Log_output::write_out("Input file format unrecognized. Only FASTA and FASTQ formats supported. Exiting.\n\n",0);
        exit(1);
    }

    map<string,int> copy_num;
    for(int i=0;i<(int)seqs.size();i++)
    {
        string name = seqs.at(i).name;

        map<string,int>::iterator it = copy_num.find(name);
        if(it != copy_num.end())
        {
            it->second++;

            stringstream ss;
            ss<<"."<<it->second;
            name.append(ss.str());

            Log_output::write_out("Warning: duplicate names found. '"+seqs.at(i).name+"' is renamed '"+name+"'.\n",2);

            seqs.at(i).name = name;
        }
        else
        {
            copy_num.insert( pair<string,int>(name,0) );
        }
    }

}

void Fasta_reader::read_fasta(istream & input, vector<Fasta_entry> & seqs, bool short_names, bool degap) const 
{

    string temp, name, comment, tmp_tid, sequence = "";  // Initialization
    int tmp_ndup;

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string


        // If first character is >
        if(temp[0] == '>')
        {
            temp = Text_utils::remove_last_whitespaces(temp);

            // If a name and a sequence were found
            if((name != "") && (sequence != ""))
            {
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
                sequence.erase( std::remove(sequence.begin(), sequence.end(), '\r'), sequence.end() );

                if(degap)
                    this->remove_gaps(&sequence);

                Fasta_entry fe;                               
                fe.name = name;
                fe.comment = comment;
                fe.sequence = sequence;
                fe.quality = "";
                fe.first_read_length = -1;
                fe.tid = tmp_tid;
                fe.cluster_attempts = 0;
                fe.num_duplicates = tmp_ndup;

                seqs.push_back(fe);
                name = "";
                sequence = "";
            }
            // Sequence name isolation
            if (! short_names)
            {
                name = temp;
                comment = "";
                tmp_tid = "";
                tmp_ndup = 1;
            }
            else
            {
                String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
                name = st->next_token();
                comment = "";
                tmp_tid = "";
                tmp_ndup = 1;
                while (st->has_more_token())
                {
                    string block = st->next_token();
                    block = Text_utils::remove_surrounding_whitespaces(block);
                    comment += " "+block;

                    if(block.substr(0,4)=="TID=")
                    {
                        block = block.substr(4);
                        tmp_tid = block;
                    }
                    if(block.substr(0,14)=="NumDuplicates=")
                    {
                        block = block.substr(14);
                        stringstream ss(block);
                        ss >> tmp_ndup;
                    }
                }
                delete st;
            }
            name.erase(name.begin());  // Character > deletion
        }
        else sequence += temp;  // Sequence isolation
    }

    // Addition of the last sequence in file
    if((name != "") && (sequence != ""))
    {
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
        sequence.erase( std::remove(sequence.begin(), sequence.end(), '\r'), sequence.end() );

        if(degap)
            this->remove_gaps(&sequence);

        Fasta_entry fe;
        fe.name = name;
        fe.comment = comment;
        fe.sequence = sequence;
        fe.quality = "";
        fe.first_read_length = -1;
        fe.tid = tmp_tid;
        fe.cluster_attempts = 0;
        fe.num_duplicates = tmp_ndup;

        seqs.push_back(fe);
    }

    if(Settings_handle::st.is("translate")
    || Settings_handle::st.is("mt-translate")
    || Settings_handle::st.is("find-orfs") )
    {
        if(this->check_sequence_data_type(&seqs) == Model_factory::dna)
        {
            vector<Fasta_entry>::iterator it = seqs.begin();
            for(;it != seqs.end(); it++)
            {
//                it->dna_sequence = it->sequence;

                string dna = it->sequence;
                this->rna_to_DNA(&dna);

                it->sequence = this->DNA_to_protein(&dna);

                it->dna_sequence = dna;
            }
        }
        else
        {
            Log_output::write_out("Option '--translate' cannot be used for proteins. Exiting.\n\n",0);
            exit(0);
        }
    }
}

void Fasta_reader::read_fastq(istream & input, vector<Fasta_entry> & seqs) const 
{
    string temp, name, comment = "";  // Initialization

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string
        // If first character is @
        if(temp[0] == '@')
        {
            temp = Text_utils::remove_last_whitespaces(temp);
            Fasta_entry fe;
            String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
            name = st->next_token();

            name.erase(name.begin());  // Character @ deletion
            fe.name = name;
            fe.num_duplicates = 1;

            comment = "";
            while (st->has_more_token())
            {
                string block = st->next_token();
                block = Text_utils::remove_surrounding_whitespaces(block);
                comment += block+" ";

                if(block.substr(0,4)=="TID=")
                {
                    block = block.substr(4);
                    fe.tid = block;
                }
                if(block.substr(0,14)=="NumDuplicates=")
                {
                    block = block.substr(14);
                    stringstream ss(block);
                    ss >> fe.num_duplicates;
                }
            }
            delete st;


            getline(input, temp, '\n');  // Copy current line in temporary string
            fe.sequence = Text_utils::to_upper( Text_utils::remove_last_whitespaces(temp) );

            getline(input, temp, '\n');  // Copy current line in temporary string
            temp = Text_utils::remove_last_whitespaces(temp);
            if(temp[0] != '+')
            {
                Log_output::write_out("Error in FASTQ comment:"+temp+"\nExiting.\n\n",0);
                exit(1);
            }

            temp.erase(temp.begin());  // Character + deletion

            if(temp.length()>0)
            {
                comment += " ; ";
                comment += temp;
            }
            fe.comment = comment;

            getline(input, temp, '\n');  // Copy current line in temporary string

            fe.quality = Text_utils::remove_last_whitespaces(temp);
            fe.first_read_length = -1;
            fe.cluster_attempts = 0;

            seqs.push_back(fe);
        }
        else if(temp != "")
        {
            Log_output::write_out("FASTQ file parse error. Expecting a line starting with '@':  \n"+temp+".\n\nExiting.\n\n",0);
            exit(1);
        }
    }

}

void Fasta_reader::read_graph(istream & input, vector<Fasta_entry> & seqs, bool short_names = false) const 
{

    string temp, name, comment, sequence, block = "";  // Initialization
    int prev_site = -1;
    vector<Seq_edge> edges;

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string

        // If first character is >
        if(temp[0] == '#')
        {
            temp.erase(temp.begin());  // Character > deletion
            temp = Text_utils::remove_first_whitespaces(temp);
            temp = Text_utils::remove_last_whitespaces(temp);

            // If a name and a sequence were found
            if((name != "") && (sequence != ""))
            {
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );

                Fasta_entry fe;
                fe.name = name;
                fe.comment = comment;
                fe.sequence = sequence;
                fe.edges = edges;
                seqs.push_back(fe);
                name = "";
                sequence = "";
            }       
            edges.clear();
            prev_site = -1;

            // Sequence name isolation
            if (! short_names)
            {
                name = temp;
                comment = "";
            }
            else
            {
                String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
                name = st->next_token();
                comment = "";
                while (st->has_more_token())
                {
                  comment += st->next_token()+" ";
                }
                delete st;
            }
        }
        else
        {
            temp = Text_utils::remove_surrounding_whitespaces(temp);
            if(temp.size()==0)
                continue;

            String_tokenizer * st = new String_tokenizer(temp, ";", true, false);
            block = st->next_token();

            String_tokenizer * bt = new String_tokenizer(block, " ", true, false);

            temp = bt->next_token();
            int site = Text_utils::to_int( temp );
            if(site != prev_site+1)
            {
                if(prev_site == -2)
                {
                    Log_output::write_out("Error reading the graph input: 'end' site is not the last site.\nExiting.\n\n",0);
                }
                else
                {
                    Log_output::write_out("Error reading the graph input: previous site "+Log_output::itos(prev_site)+
                                          " and this site "+Log_output::itos(site)+".\nExiting.\n\n",0);
                }
                exit(1);
            }
            prev_site++;

            temp = bt->next_token();
            temp = Text_utils::remove_surrounding_whitespaces(temp);
            if(temp == "start")
            {
                if(site != 0)
                {
                    Log_output::write_out("Error reading the graph input: 'start' is not the site 0.\nExiting.\n\n",0);
                    exit(1);
                }
            }
            else if(temp == "end")
            {
                prev_site = -2;
            }
            else
            {
                sequence += temp[0];
            }

            delete bt;

            // Next block
            double sum_weight = 0;

            block = st->next_token();
            block = Text_utils::remove_surrounding_whitespaces(block);

            while(block != "")
            {

                String_tokenizer * bt = new String_tokenizer(block, " ", true, false);

                int start_site = Text_utils::to_int( bt->next_token() );
                int end_site = Text_utils::to_int( bt->next_token() );
                double weight = Text_utils::to_double( bt->next_token() );

                if(start_site < 0 || end_site < start_site || end_site > site)
                {
                    Log_output::write_out("Error reading the graph input: edge coordinates at site "+Log_output::itos(site)+" appear incorrect.\nExiting.\n\n",0);
                    exit(1);
                }

                if(weight < 0 || weight > 1 || weight + sum_weight > 1)
                {
                    Log_output::write_out("Warning reading the graph input: edge weight at site "+Log_output::itos(site)+" appear incorrect.\n\n",2);
                }

                if(weight + sum_weight > 1)
                {
                    Log_output::write_out("Warning reading the graph input: edge weight at site "+Log_output::itos(site)+" appear incorrect.\n\n",2);
                    weight = 1.0-sum_weight;
                    sum_weight = 1.0;
                }

                Seq_edge edge;
                edge.start_site = start_site;
                edge.end_site = end_site;
                edge.weight = weight;

                edges.push_back(edge);

                delete bt;

                block = st->next_token();
                block = Text_utils::remove_surrounding_whitespaces(block);
            }
        }
    }


    // Addition of the last sequence in file
    if((name != "") && (sequence != ""))
    {
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );

        Fasta_entry fe;
        fe.name = name;
        fe.comment = comment;
        fe.sequence = sequence;
        fe.edges = edges;
        seqs.push_back(fe);
    }

}

/****************************************************************************************/

void Fasta_reader::write(ostream & output, const vector<Fasta_entry> & seqs, string format) const 
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Fasta_reader::write. Failed to open file"); }

    string seq, temp = "";  // Initialization

    multimap<string,int> copy_num;
    vector<Fasta_entry>::const_iterator vi1 = seqs.begin();
    for (; vi1 != seqs.end(); vi1++)
    {
        vector<Fasta_entry>::const_iterator vi2 = vi1;
        vi2++;

        int copy = 1;
        for (; vi2 != seqs.end(); vi2++)
        {
            if(vi1->name == vi2->name)
            {
                if(copy==1)
                    copy_num.insert( pair<string,int>(vi1->name,copy) );

                copy++;
                copy_num.insert( pair<string,int>(vi2->name,copy) );
            }
        }
    }

    vector<Fasta_entry> seqs2;

    vector<Fasta_entry>::const_iterator vi = seqs.begin();
    for (; vi != seqs.end(); vi++)
    {
        Fasta_entry e = *vi;
        multimap<string,int>::iterator cit = copy_num.find(vi->name);
        if(cit!=copy_num.end())
        {
            stringstream ss;
            ss << e.name << "/" << cit->second;
            e.name = ss.str();
            copy_num.erase(cit);
        }
        seqs2.push_back(e);
    }

    if(format == "fasta")
        write_fasta(output,seqs2);
    else if(format == "raxml")
        this->write_long_sequential(output,seqs2);
    else if (format == "phylipi")
        this->write_interleaved(output,seqs2);
    else if (format == "phylip")
        this->write_interleaved(output,seqs2);
    else if (format == "phylips")
        this->write_sequential(output,seqs2,true);
    else if (format == "nexus")
        this->write_simple_nexus(output,seqs2);
    else if (format == "paml")
        this->write_sequential(output,seqs2,false);
    else
    {
        Log_output::write_out("Outformat '"+format+"' not recognised.\nAccepted formats: fasta, raxml, paml, phylips, phylipi, nexus.\n Using 'fasta'.",0);
        write_fasta(output,seqs2);
    }
}

string Fasta_reader::get_format_suffix(string format) const 
{
    if(format == "raxml")
        return ".phy";
    else if (format == "phylip")
        return ".phy";
    else if (format == "phylipi")
        return ".phy";
    else if (format == "phylips")
        return ".phy";
    else if (format == "nexus")
        return ".nex";
    else if (format == "paml")
        return ".phy";
    else
        return ".fas";

}
void Fasta_reader::write_fasta(ostream & output, const vector<Fasta_entry> & seqs) const 
{
    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        output << ">" << vi->name;

        if(vi->comment != "")
            output << vi->comment;

        output << endl;

        // Sequence cutting to specified characters number per line
        string seq = vi->sequence;
        string temp;
        while (seq != "")
        {
            if (seq.size() > chars_by_line)
            {
                temp = string(seq.begin(), seq.begin() + chars_by_line);
                output << temp  << endl;
                seq.erase(seq.begin(), seq.begin() + chars_by_line);
            }
            else
            {
                output << seq << endl;
                seq = "";
            }
        }
    }
}


void Fasta_reader::write_interleaved(ostream & output, const vector<Fasta_entry> & seqs) const 
{
    vector<Fasta_entry>::const_iterator vi = seqs.begin();
    int length = vi->sequence.length();
    output << seqs.size() << " " << length << endl;

    for (int offset = 0; offset<length; offset+=chars_by_line)
    {
        vi = seqs.begin();

        for (; vi!=seqs.end(); vi++)
        {
            string tmp = vi->name.substr(0,10)+"          ";
            if (offset > 0)
            {
                tmp = "           ";
            }
            output << tmp.substr(0,10)<<" ";

            output<<vi->sequence.substr(offset,chars_by_line)<<endl;
        }
    }
}

void Fasta_reader::write_sequential(ostream & output, const vector<Fasta_entry> & seqs, bool truncate) const 
{

    vector<Fasta_entry>::const_iterator vi = seqs.begin();
    int length = vi->sequence.length();
    output << seqs.size() << " " << length << endl;

    for (; vi!=seqs.end(); vi++)
    {
        if (truncate)
            output << (vi->name+"          ").substr(0,10)<<" "<< endl;
        else
            output << vi->name << endl;

        // Sequence cutting to specified characters number per line
        string seq = vi->sequence;
        string temp;
        while (seq != "")
        {
            if ((int)seq.size() > (int)chars_by_line)
            {
                temp = string(seq.begin(), seq.begin() + chars_by_line);
                output << temp  << endl;
                seq.erase(seq.begin(), seq.begin() + chars_by_line);
            }
            else
            {
                output << seq << endl;
                seq = "";
            }
        }
    }
}

void Fasta_reader::write_long_sequential(ostream & output, const vector<Fasta_entry> & seqs) const 
{

    vector<Fasta_entry>::const_iterator vi = seqs.begin();
    int length = vi->sequence.length();
    output << seqs.size() << " " << length << endl;

    for (; vi!=seqs.end(); vi++)
        output << vi->name<< endl<<vi->sequence<<endl;
}

void Fasta_reader::write_simple_nexus(ostream & output, const vector<Fasta_entry> & seqs) const 
{

    vector<Fasta_entry>::const_iterator vi = seqs.begin();
    int length = vi->sequence.length();

    string datatype = "protein";

    if(this->check_sequence_data_type(&seqs) == Model_factory::dna)
        datatype = "dna";

    output<<"#NEXUS\nbegin data;\ndimensions ntax="<<seqs.size()<<" nchar="<<length<<";\nformat datatype="<<datatype<<" interleave=yes gap=-;\nmatrix\n"<<endl;


    for (int offset = 0; offset<length; offset+=chars_by_line)
    {
        output<<endl;

        vi = seqs.begin();

        for (; vi!=seqs.end(); vi++)
        {
            string tmp = vi->name.substr(0,20)+"'                    ";
            output << "'"<<tmp.substr(0,21)<<"     ";

            output<<vi->sequence.substr(offset,chars_by_line)<<endl;
        }
    }
    output<<";\nend;";

}

/****************************************************************************************/

void Fasta_reader::get_DNA_seqs(Node *root, const vector<Fasta_entry> *org_seqs, map<string,string> *dna_seqs)
{
    Text_utils tu;
    vector<Fasta_entry>::const_iterator it = org_seqs->begin();

    for (; it != org_seqs->end(); it++)
    {
        string seq = it->dna_sequence;
        tu.replace_all(seq,"-","");

        dna_seqs->insert(pair<string,string>(it->name,seq));
    }

    vector<Node*> read_nodes;
    root->get_read_nodes_below(&read_nodes);

    vector<Node*>::iterator nit = read_nodes.begin();
    for (; nit != read_nodes.end(); nit++)
    {
        if((*nit)->has_Orf())
        {
            dna_seqs->insert(pair<string,string>((*nit)->get_name(),(*nit)->get_Orf().dna_sequence));
        }
    }

}

void Fasta_reader::backtranslate_dna(const vector<Fasta_entry> & seqs, const map<string,string> *dna_seqs, vector<Fasta_entry> &outseqs, bool include_mock_ancestors) const 
{
    bool dna_seq_missing = false;

    map<string,string>::const_iterator it = dna_seqs->begin();

    for (; it != dna_seqs->end(); it++)
    {
        if( it->second.length() == 0 )
            dna_seq_missing = true;
    }

    if(dna_seq_missing)
    {
        Log_output::write_out("No DNA for all sequences. Back-translation failed.\n",0);
        return;
    }

    vector<Fasta_entry>::const_iterator vi;

    // Main loop : for all sequences in vector container
    for (vi = seqs.begin(); vi != seqs.end(); vi++)
    {
        Fasta_entry os;
        os.name = vi->name;
        os.comment = vi->comment;

        map<string,string>::const_iterator it2 = dna_seqs->find(vi->name);
        if(it2 == dna_seqs->end())
        {
            if(include_mock_ancestors)
            {
                string prot = vi->sequence;
                os.sequence = protein_to_mockDNA(&prot);
                os.num_duplicates = 1;
            }
            else
            {
                Log_output::write_out("No matching DNA sequence for "+vi->name+". Back-translation failed.\n",1);
                return;
            }
        }
        else
        {
            string dna = it2->second;
            string prot = vi->sequence;

            os.sequence = protein_to_DNA(&dna,&prot);
            os.num_duplicates = 1;
        }

        outseqs.push_back(os);
    }
}

void Fasta_reader::write_dna(ostream & output, const vector<Fasta_entry> & seqs, const vector<Fasta_entry> & org_seqs, Node *root, int output_type) const 
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Fasta_reader::write_dna. Failed to open file"); }

    bool dna_seq_missing = false;
    vector<Fasta_entry>::const_iterator it = org_seqs.begin();
    map<string,string> dna_seqs;

    Text_utils tu;

    for (; it != org_seqs.end(); it++)
    {
        if( it->dna_sequence.length() == 0 )
            dna_seq_missing = true;

        string seq = it->dna_sequence;
        tu.replace_all(seq,"-","");

        dna_seqs.insert(pair<string,string>(it->name,seq));
    }

    if(dna_seq_missing)
    {
        Log_output::write_out("No DNA for all sequences. Back-translation failed.\n",1);
        return;
    }

    if(Settings_handle::st.is("find-orfs"))
    {
        root->get_read_dna_sequences(&dna_seqs);
    }

    multimap<string,int> copy_num;
    vector<Fasta_entry>::const_iterator vi1 = seqs.begin();
    for (; vi1 != seqs.end(); vi1++)
    {
        vector<Fasta_entry>::const_iterator vi2 = vi1;
        vi2++;

        int copy = 1;
        for (; vi2 != seqs.end(); vi2++)
        {
            if(vi1->name == vi2->name)
            {
                if(copy==1)
                    copy_num.insert( pair<string,int>(vi1->name,copy) );

                copy++;
                copy_num.insert( pair<string,int>(vi2->name,copy) );
            }
        }
    }

    vector<Fasta_entry> outseqs;

    vector<Fasta_entry>::const_iterator vi;

    // Main loop : for all sequences in vector container
    for (vi = seqs.begin(); vi != seqs.end(); vi++)
    {
        Fasta_entry os;
        os.name = vi->name;

        multimap<string,int>::iterator cit = copy_num.find(vi->name);
        if(cit!=copy_num.end())
        {
            os.name += "/"+Text_utils::to_string(cit->second);
            copy_num.erase(cit);
        }
        if(vi->comment != "")
            os.comment = vi->comment;


        map<string,string>::iterator it2 = dna_seqs.find(vi->name);
        if(it2 == dna_seqs.end())
        {
            Log_output::write_out("No matching DNA sequence for "+vi->name+". Back-translation failed.\n",1);
            return;
        }
        else
        {
            string dna = it2->second;
            string prot = vi->sequence;

            os.sequence = protein_to_DNA(&dna,&prot);
            os.num_duplicates = 1;
        }

        outseqs.push_back(os);
    }

    set<string> read_names;
    root->get_read_node_names(&read_names);

    Fasta_entry entry;
    entry.name = "consensus";
    entry.sequence = "";

    int min_num_seqs = int( root->get_weighted_number_of_leaves() * (float)Settings_handle::st.get("consensus-minimum-proportion").as<float>() );

    if(min_num_seqs < Settings_handle::st.get("consensus-minimum").as<int>())
        min_num_seqs = Settings_handle::st.get("consensus-minimum").as<int>();

    if(output_type != Fasta_reader::plain_alignment)
    {

        for(int i=0;i<(int)outseqs.at(0).sequence.length();i++)
        {
            int sA=0; int sC=0; int sG=0; int sT=0;
            bool included_in_reference = false;
            char parent_char = '-';

            for (vi = outseqs.begin(); vi != outseqs.end(); vi++)
            {
                char x = vi->sequence.at(i);
                int nd = vi->num_duplicates;

                if(read_names.find(vi->name) != read_names.end())
                {
                    if(x == 'A')
                        sA += nd;
                    else if(x == 'C')
                        sC += nd;
                    else if(x == 'G')
                        sG += nd;
                    else if(x == 'T')
                        sT += nd;
                }
                else
                {
                    if(x != '-')
                    {
                        included_in_reference = true;
                        if(Settings_handle::st.is("show-contig-ancestor"))
                            parent_char = tolower(x);
                    }
                }
            }

            if(!included_in_reference && sA+sC+sG+sT < min_num_seqs)
            {
                entry.sequence.append("-");
            }
            else
            {
                if(sA==0 && sC==0 && sG==0 && sT==0)
                    entry.sequence.append(string(1,parent_char));
                else if(sA>sC && sA>sG && sA>sT)
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

    if(output_type != Fasta_reader::consensus_only || Settings_handle::st.is("inlude-parent-in-contig") )
    {
        for (vi = outseqs.begin(); vi != outseqs.end(); vi++)
        {
            if(read_names.find(vi->name) == read_names.end())
            {
                this->print_fasta_entry(output,&(*vi));
            }
        }
    }

     if(output_type != Fasta_reader::plain_alignment)
            this->print_fasta_entry(output,&entry);

     if(output_type != Fasta_reader::consensus_only)
     {
        for (vi = outseqs.begin(); vi != outseqs.end(); vi++)
        {
            if(read_names.find(vi->name) != read_names.end())
            {
                this->print_fasta_entry(output,&(*vi));
            }
        }
     }

}

void Fasta_reader::print_fasta_entry(ostream & output, const Fasta_entry *entry) const
{
    output << ">" << entry->name << " " <<entry->comment<< endl;
    string seq = entry->sequence;

    string temp = "";  // Initialization
    while (seq != "")
    {
        if (seq.size() > chars_by_line)
        {
            temp = string(seq.begin(), seq.begin() + chars_by_line);
            output << temp << endl;
            seq.erase(seq.begin(), seq.begin() + chars_by_line);
        }
        else
        {
            output << seq << endl;
            seq = "";
        }
    }
}

/****************************************************************************************/

void Fasta_reader::write_fastq(ostream & output, const vector<Fasta_entry> & seqs) const 
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Fasta_reader::write_fastq. Failed to open file"); }

    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        output << "@" << vi->name << endl;
        output << vi->sequence << endl;
        output << "+" << vi->comment << endl;
        output << vi->quality << endl;
    }

}

void Fasta_reader::write_graph(ostream & output, Node * root) const 
{

    Sequence *sequence = root->get_sequence();
    output<<"# root node\n";

    string alpha = sequence->get_full_alphabet();

    output<<"0 start;\n";
    for(int i=1; i<sequence->sites_length()-1; i++)
    {
        Site *site = sequence->get_site_at(i);
        output<<i<<" "<<alpha.at(site->get_state())<<"; ";

        Edge *edge = site->get_first_bwd_edge();
        output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();
            output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";
        }
        output<<endl;

    }
    output<<sequence->sites_length()-1<<" end; ";

    Site *site = sequence->get_site_at(sequence->sites_length()-1);

    Edge *edge = site->get_first_bwd_edge();
    output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";

    while(site->has_next_bwd_edge())
    {
        edge = site->get_next_bwd_edge();
        output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";
    }
    output<<endl;

}

/****************************************************************************************/

void Fasta_reader::remove_gap_only_columns(vector<Fasta_entry> *sequences)  
{
    int length = sequences->at(0).sequence.length();
    vector<Fasta_entry>::iterator vi = sequences->begin();
    for (; vi != sequences->end(); vi++)
    {
        if((int)vi->sequence.length() != length)
        {
            Log_output::write_out("Error: aligned sequences of different length. Removal of gap only columns failed.\n",0);
            return;
        }
    }


    for(int i=0;i<length;)
    {
        bool gap_only = true;
        for (vi = sequences->begin(); vi != sequences->end(); vi++)
        {
            if(vi->sequence.at(i) != '-')
                gap_only = false;
        }
        if(gap_only)
        {
            for (vi = sequences->begin(); vi != sequences->end(); vi++)
            {
                vi->sequence.erase(i,1);
            }
            length = sequences->at(0).sequence.length();
        }
        else
        {
            i++;
        }
    }
}

/****************************************************************************************/

void Fasta_reader::remove_gaps(string *seq) const 
{
    /*
    for(int i=0;i<(int)seq->length();)
    {
        if(seq->at(i) == '-')
            seq->erase(i,1);
        else
            i++;
    }*/
    string& sequence = *seq;
    stringstream seqss;
    for(int i = 0; i < (int)sequence.length(); i++)
    {
        if(sequence[i] != '-')
        {
            seqss << sequence[i];
        }
    }
    sequence = seqss.str();
}

void Fasta_reader::remove_gaps(vector<Fasta_entry> *seqs) const 
{
    for(int i=0;i<(int)seqs->size();i++)
        this->remove_gaps(&seqs->at(i).sequence);

}

/****************************************************************************************/

bool Fasta_reader::check_alphabet(vector<Fasta_entry> * sequences,int data_type) 
{

    bool allow_gaps = Settings_handle::st.is("ref-seqfile");

    if(data_type<0)
        data_type = this->check_sequence_data_type(sequences);

    // Check that alphabet is correct but use (faster?) build-in alphabet.
    //
    if(data_type == Model_factory::dna)
    {

        string full_alphabet = Model_factory::get_dna_full_char_alphabet();

        dna_pi[0] = dna_pi[1] = dna_pi[2] = dna_pi[3] = 0.0;

        bool chars_ok = true;
        vector<Fasta_entry>::iterator vi = sequences->begin();

        // Main loop : for all sequences in vector container
        for (; vi != sequences->end(); vi++)
        {
            vi->data_type = Model_factory::dna;

            // Convert U -> T and all uppercase
            this->rna_to_DNA(&vi->sequence);

            string::iterator si = vi->sequence.begin();
            for (;si != vi->sequence.end();si++)
            {
                char c = *si;
                switch (c)
                {
                    case 'A':
                        dna_pi[0]++;
                        break;
                    case 'C':
                        dna_pi[1]++;
                        break;
                    case 'G':
                        dna_pi[2]++;
                        break;
                    case 'T':
                        dna_pi[3]++;
                        break;
                    case '-':
                        if(!allow_gaps)
                        {
                            vi->sequence.erase(si);
                            si--;
                            chars_ok = false;
                        }
                        break;
                    default:
                        // Remove characters not in full alphabet
                        if(full_alphabet.find(c) == string::npos) {
                            vi->sequence.erase(si);
                            si--;
                            if(c!=' ')
                            {
                                chars_ok = false;
                            }
                       }
                }
            }
        }

        float tot = dna_pi[0]+dna_pi[1]+dna_pi[2]+dna_pi[3];

        dna_pi[0] /= tot;
        dna_pi[1] /= tot;
        dna_pi[2] /= tot;
        dna_pi[3] /= tot;

        return chars_ok;
    }

    else if(data_type == Model_factory::protein)
    {
        string full_alphabet = Model_factory::get_protein_full_char_alphabet();

        bool chars_ok = true;
        vector<Fasta_entry>::iterator vi = sequences->begin();

        // Main loop : for all sequences in vector container
        for (; vi != sequences->end(); vi++)
        {
            vi->data_type = Model_factory::protein;

            string::iterator si = vi->sequence.begin();
            for (;si != vi->sequence.end();si++)
            {
                char c = *si;

                // Remove characters not in full alphabet
                //
                if(full_alphabet.find(c) == string::npos) {

                    if(c=='U')
                        *si='X';

                    else if(c!='-' || !allow_gaps)
                    {
                        vi->sequence.erase(si);
                        si--;
                        if(c!=' ')
                        {
                            chars_ok = false;
                        }
                    }
                }
            }
        }

        return chars_ok;
    }

    return false;
}

/****************************************************************************************/

int Fasta_reader::check_sequence_data_type(const vector<Fasta_entry> *sequences) const
{

    vector<Fasta_entry>::const_iterator vi = sequences->begin();

    int dna = 0;
    int protein = 0;
    string dna_alphabet = "ACGTUN";
    string protein_alphabet = Model_factory::get_protein_char_alphabet();

    // Main loop : for all sequences in vector container
    for (; vi != sequences->end(); vi++)
    {
        string::const_iterator si = vi->sequence.begin();
        for (;si != vi->sequence.end();si++)
        {
            char c = *si;
            if(dna_alphabet.find(c) != string::npos)
            {
                dna++;
            }
            if(protein_alphabet.find(c) != string::npos)
            {
                protein++;
            }
        }
    }

    if( ((float)dna )/ (float)protein > 0.9)
        return Model_factory::dna;
    else
        return Model_factory::protein;
}


/****************************************************************************************/

void Fasta_reader::rna_to_DNA(string *sequence) const
{
//    transform( sequence->begin(), sequence->end(), sequence->begin(), (int(*)(int))toupper );
    size_t si = 0;
    si=sequence->find('U',si);
    while(si != string::npos) {
        sequence->replace(si,1,1,'T');
        si=sequence->find('U',si);
    }
}

void Fasta_reader::define_translation_tables()
{
    string codon[66] = {"TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                        "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                        "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                        "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                        "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
                        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
                        "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
                        "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
                        "NNN", "---"
                       };
    string unaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "X", "W", "R", "R", "R", "R", "S", "S", "R", "R", "G", "G", "G", "G",
                       "X", "-"
                      };
    string mtaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "M", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "W", "W", "R", "R", "R", "R", "S", "S", "X", "X", "G", "G", "G", "G",
                       "X", "-"
                      };

    if(Settings_handle::st.is("mt-translate"))
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],mtaa[i]));
            aa_to_codon.insert(make_pair(mtaa[i],codon[i]));
        }
    }
    else
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],unaa[i]));
            aa_to_codon.insert(make_pair(unaa[i],codon[i]));
        }
    }
}

string Fasta_reader::DNA_to_protein(string *sequence) const
{
    string prot;

    for (unsigned int j=0; j<sequence->length(); j+=3)
    {
        string codon = sequence->substr(j,3);
        if (codon_to_aa.find(codon) == codon_to_aa.end())
        {
            sequence->replace(j,3,"NNN");
            prot += "X";
        }
        else
        {
            prot += codon_to_aa.find(codon)->second;
        }
    }

    return prot;
}


string Fasta_reader::protein_to_DNA(string *dna,string *prot) const
{
    string out;
    int pos = 0;
    for (unsigned int j=0; j<prot->length(); j++)
    {
        string aa = prot->substr(j,1);
        if(aa=="-")
            out += "---";
        else
        {
            if(pos+3<=(int)dna->length())
                out += dna->substr(pos,3);
            else
                out += "NNN";
            pos += 3;
        }
    }

    return out;
}

string Fasta_reader::protein_to_mockDNA(string *prot) const
{
    string out;
    for (unsigned int j=0; j<prot->length(); j++)
    {
        string aa = prot->substr(j,1);
        if(aa=="-")
            out += "---";
        else
        {
            out += "NNN";
        }
    }

    return out;
}

/****************************************************************************************/

bool Fasta_reader::check_sequence_names(const vector<Fasta_entry> *sequences,const vector<Node*> *leaf_nodes) const
{
    unsigned int names_match = 0;
    vector<Fasta_entry>::const_iterator si = sequences->begin();
    for (; si != sequences->end(); si++)
    {
        string s_name = si->name;
        vector<Node*>::const_iterator ni = leaf_nodes->begin();
        for (; ni != leaf_nodes->end(); ni++)
        {
            if((*ni)->get_name() == s_name)
                names_match++;
        }
    }

    // All the leafs in the guidetree need a sequence. Not all sequences need to be in the tree.

    int overlap = 0;

    if(names_match != leaf_nodes->size())
    {
        if(Settings_handle::st.is("ref-seqfile") && Settings_handle::st.is("ref-treefile"))
        {
            Log_output::write_out("\nAll leaf node names in the guidetree file '"+Settings_handle::st.get("ref-treefile").as<string>()+"' do not match\n"+
                " with a sequence name in the sequence file '"+Settings_handle::st.get("ref-seqfile").as<string>()+"'.\n",2);
        }
        else
        {
            Log_output::write_out("\nAll leaf node names in the guidetree file '"+Settings_handle::st.get("treefile").as<string>()+"' do not match\n"+
                " with a sequence name in the sequence file '"+Settings_handle::st.get("seqfile").as<string>()+"'.\n",2);
        }

        set<string> snames;
        for (si = sequences->begin(); si != sequences->end(); si++)
        {
            snames.insert(si->name);
        }

        set<string> tnames;
        for (vector<Node*>::const_iterator ni = leaf_nodes->begin(); ni != leaf_nodes->end(); ni++)
        {
            tnames.insert((*ni)->get_name());
        }


        stringstream ss;
        ss<<endl<<"Leaf names not in the sequences:"<<endl;
        for (vector<Node*>::const_iterator ni = leaf_nodes->begin(); ni != leaf_nodes->end(); ni++)
        {
            if(snames.find((*ni)->get_name()) == snames.end() )
            {
                ss<<" "<<(*ni)->get_name()<<endl;
                (*ni)->has_sequence(false);
            }
            else
            {
                (*ni)->has_sequence(true);
                overlap++;
            }
        }
        Log_output::write_out(ss.str(),2);
    }
    else
    {
        overlap = names_match;
    }

    // Not all sequences need to be in the tree but give a warning that names don't match.
    if(names_match < 1)
    {
        Log_output::write_out("\nSequence names do not match the names in the tree! Exiting.\n\n",0);
        exit(1);
    }

//    if((int)sequences->size() > overlap && overlap == (int)leaf_nodes->size())
    if((int)sequences->size() > overlap )
    {
        Log_output::write_out("\nWarning: "+Log_output::itos(leaf_nodes->size())+" leaf nodes and "+
                              Log_output::itos(sequences->size())+" sequences! Sequences not matching will be removed.\n\n",0);
//        return true;
    }
    if(overlap < (int)leaf_nodes->size())
    {
        Log_output::write_out("\nWarning: "+Log_output::itos(leaf_nodes->size())+" leaf nodes but "+
                              Log_output::itos(overlap)+" matching sequences! Excess branches will be removed.\n\n",0);
        return false;
    }
    return true;
}

/****************************************************************************************/

void Fasta_reader::place_sequences_to_nodes(const vector<Fasta_entry> *sequences,vector<Node*> *leaf_nodes, bool gapped, int data_type)
{

    if(data_type<0)
        data_type = this->check_sequence_data_type(sequences);

    vector<Fasta_entry>::const_iterator si = sequences->begin();
    for (; si != sequences->end(); si++)
    {
        string s_name = si->name;

        vector<Node*>::iterator ni = leaf_nodes->begin();
        for (; ni != leaf_nodes->end(); ni++)
        {
            if((*ni)->get_name() == s_name) {
                (*ni)->add_name_comment( si->comment );
                (*ni)->add_sequence( *si, data_type, gapped);
            }
        }
    }
}

/****************************************************************************************/

void Fasta_reader::read_bpp_phylip(const char* filename,map<string,string> *sequences)
{
    ifstream input(filename, ios::in);

    if (! input) { throw IOException ("Fasta_reader::read_bpp_phylip. Failed to open file"); }

    this->read_bpp_phylip(input,sequences);
}

void Fasta_reader::read_bpp_phylip(istream & input,map<string,string> *sequences)
{
    int nseq = -1;
    int length = -1;

    string temp, name, sequence = "";
    getline(input, temp, '\n');

    stringstream nums(temp);
    nums>>nseq>>length;

    if (nseq<1 || length<1) { throw IOException ("Fasta_reader::read_bpp_phylip. Error reading the file"); }

    for (int i=0; i<nseq; i++)
    {
        getline(input, temp, '\n');

        stringstream rows(temp);
        rows>>name>>sequence;

        name = Text_utils::remove_last_whitespaces(name);
        sequence = Text_utils::remove_whitespaces(sequence);

        do
        {
            getline(input, temp, '\n');

            temp = Text_utils::remove_whitespaces(temp);
            sequence += temp;
        }
        while (temp.length()>0);

        sequences->insert(sequences->begin(),pair<string,string>("#"+name+"#",sequence) );
    }
}
