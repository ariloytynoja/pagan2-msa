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

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <time.h>
#include "utils/exceptions.h"
#include "main/node.h"
#include "utils/fasta_entry.h"

using namespace std;

namespace ppa
{

class Fasta_reader
{
    unsigned int chars_by_line;
    float dna_pi[4];

    std::map<std::string,std::string> codon_to_aa;
    std::map<std::string,std::string> aa_to_codon;

    void rna_to_DNA(string *sequence) const;
    void define_translation_tables();

    string DNA_to_protein(string *sequence) const;
    string protein_to_mockDNA(string *prot) const;
    string protein_to_DNA(string *dna,string *prot) const;
    void ensure_writeable_file(ofstream& output, string& suffix) const
    {
        if(! output)
        {

            time_t timestamp;
            time(&timestamp);
            stringstream fnamess;
            fnamess << "output" << timestamp << suffix;
            string filename = fnamess.str();

            Log_output::write_out("Invalid directory or file name. Automatic outputname used instead: "+filename+"\n",0);

            output.open((filename).c_str(), ios::out);

        }
    }
public:

    enum output_mode {plain_alignment,contig_alignment,consensus_only};

    Fasta_reader() : chars_by_line(60) {
        this->define_translation_tables();
    }

    void set_chars_by_line(int n) { chars_by_line = n; }

    void read(istream & input, vector<Fasta_entry> & seqs, bool short_names, bool degap) const ;
    void read(const string & path, vector<Fasta_entry> & seqs, bool short_names=false, bool degap=false) const 
    {
        ifstream input(path.c_str(), ios::in);
        read(input, seqs, short_names,degap);
        input.close();
    }
    void read_fasta(istream & input, vector<Fasta_entry> & seqs, bool short_names=false, bool degap=false) const ;
    void read_fastq(istream & input, vector<Fasta_entry> & seqs) const ;
    void read_graph(istream & input, vector<Fasta_entry> & seqs, bool short_names) const ;

    void trim_fastq_reads(vector<Fasta_entry> * seqs) const ;

    void write(ostream & output, const vector<Fasta_entry> & seqs, string format) const ;
    void write(const string & path, const vector<Fasta_entry> & seqs, string format, bool overwrite=true) const 
    {
        string suffix = this->get_format_suffix(format);
        ofstream output( (path+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        ensure_writeable_file(output, suffix);
        write(output, seqs, format);
        output.close();
    }

    string get_format_suffix(string format) const ;
    void write_fasta(ostream & output, const vector<Fasta_entry> & seqs) const ;
    void write_interleaved(ostream & output, const vector<Fasta_entry> & seqs) const ;
    void write_sequential(ostream & output, const vector<Fasta_entry> & seqs, bool truncate) const ;
    void write_long_sequential(ostream & output, const vector<Fasta_entry> & seqs) const ;
    void write_simple_nexus(ostream & output, const vector<Fasta_entry> & seqs) const ;

    void write_dna(ostream & output, const vector<Fasta_entry> & seqs, const vector<Fasta_entry> & org_seqs, Node *root, int output_type=Fasta_reader::plain_alignment) const ;
    void write_dna(const string & path, const vector<Fasta_entry> & seqs, const vector<Fasta_entry> & org_seqs, Node *root, bool overwrite=true, int output_type=Fasta_reader::plain_alignment) const 
    {
        string suffix = ".fas";
        ofstream output( (path+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app) );
        ensure_writeable_file(output, suffix);
        write_dna(output, seqs, org_seqs, root, output_type);
        output.close();
    }

    void get_DNA_seqs(Node *root, const vector<Fasta_entry> *org_seqs, map<string,string> *dna_seqs);
    void backtranslate_dna(const vector<Fasta_entry> & seqs, const map<string,string> *dna_seqs, vector<Fasta_entry> &outseqs, bool include_mock_ancestors=false) const ;

    void print_fasta_entry(ostream & output, const Fasta_entry *entry) const;

    void write_fastq(ostream & output, const vector<Fasta_entry> & seqs) const ;
    void write_fastq(const string & path, const vector<Fasta_entry> & seqs, bool overwrite=true) const 
    {
        string suffix = ".fastq";
        ofstream output( (path+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        ensure_writeable_file(output, suffix);
        write_fastq(output, seqs);
        output.close();
    }

    void write_anctree(const string outfile, Node *root, bool overwrite=true)
    {
        string suffix = ".anctree";
        ofstream output( (outfile+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        ensure_writeable_file(output, suffix);
        if (! output) { throw IOException ("Fasta_reader::write_anctree. Failed to open file"); }

        output<<root->print_tree(true);
        output.close();
    }

    void write_graph(ostream & output, Node * root) const ;
    void write_graph(const string & path, Node * root, bool overwrite=true) const 
    {
        string suffix = ".grp";
        ofstream output( (path+suffix).c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        ensure_writeable_file(output, suffix);
        write_graph(output, root);
        output.close();
    }

    void remove_gap_only_columns(vector<Fasta_entry> *sequences)  ;
    void remove_gaps(string *seq) const ;
    void remove_gaps(vector<Fasta_entry> *seqs) const ;

    bool check_alphabet(vector<Fasta_entry> *sequences, int data_type = -1)  ;
    bool check_sequence_names(const vector<Fasta_entry> *sequences,const vector<Node*> *leaf_nodes) const;

    float* base_frequencies() { return dna_pi; }
    int check_sequence_data_type(const vector<Fasta_entry> * sequences) const;

    void place_sequences_to_nodes(const vector<Fasta_entry> *sequences,vector<Node*> *leaf_nodes, bool gapped = false, int data_type = -1);

    void read_bpp_phylip(const char* filename,map<string,string> *sequences);
    void read_bpp_phylip(istream & input,map<string,string> *sequences);

};
}

#endif // FASTA_READER_H

