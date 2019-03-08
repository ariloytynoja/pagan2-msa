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

#include "input_output_parser.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"
#include "utils/newick_reader.h"
#include "utils/fasta_reader.h"
#include "utils/mafft_alignment.h"
#include "utils/bppancestors.h"
#include "utils/bppdist_tree.h"
#include "utils/bppphysamp_tree.h"
#include "utils/fasttree_tree.h"
#include "utils/raxml_tree.h"
#include "utils/exonerate_queries.h"
#include "utils/xml_writer.h"
#include "utils/model_factory.h"
#include "utils/evol_model.h"
#include "main/node.h"
#include "utils/tree_node.h"
#include "main/reads_aligner.h"
#include "utils/substring_hit.h"
#include "utils/codon_translation.h"
#include <iostream>
#include <vector>

using namespace ppa;
using namespace std;


Input_output_parser::Input_output_parser()
{
}

/************************************************************************************/

void Input_output_parser::parse_input_sequences(Fasta_reader *fr,vector<Fasta_entry> *sequences, bool *reference_alignment)
{

    /***********************************************************************/
    /*  Read the sequences                                                 */
    /***********************************************************************/

    if(Settings_handle::st.is("seqfile") && Settings_handle::st.is("queryfile"))
    {
        Log_output::write_out("Incorrect command: extension ('--queryfile') requires a reference alignment ('--ref-seqfile').\nExiting.\n\n",0);
        exit(1);
    }

    if(Settings_handle::st.is("seqfile"))
    {
        string seqfile =  Settings_handle::st.get("seqfile").as<string>();
        Log_output::write_out("Sequence file: "+seqfile+"\n",1);

        try
        {
            fr->read(seqfile, *sequences, true, true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the sequence file '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }
    }
    else if(Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("ref-seqfile").as<string>();
        Log_output::write_out("Reference alignment file: "+seqfile+"\n",1);

        try
        {
            Log_output::write_header("Reading reference alignment",0);
            fr->read(seqfile, *sequences, true, false);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the reference alignment  file '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        *reference_alignment = true;
    }
    else if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("queryfile").as<string>();
        Log_output::write_out("Reference sequence from: "+seqfile+"\n",1);

        try
        {
            Log_output::write_header("Reading reference alignment",0);
            fr->read(seqfile, *sequences, true,true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the queryfile '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        if(Settings_handle::st.is("use-duplicate-weigths"))
        {
            Fasta_entry heaviest = *sequences->begin();

            vector<Fasta_entry>::iterator it = sequences->begin();
            for(;it!=sequences->end();)
            {
                if(it->num_duplicates > heaviest.num_duplicates)
                    heaviest = *it;

                sequences->erase(it);
            }

            sequences->push_back(heaviest);
        }
        else
        {
            vector<Fasta_entry>::iterator it = sequences->begin();
            it++;
            for(;it!=sequences->end();)
                sequences->erase(it);
        }

       *reference_alignment = true;
    }
    else
    {
        Log_output::write_out("\nError: No sequence file defined.\n",0);
        Settings_handle::st.info();

        exit(1);
    }
}

/************************************************************************************/

Node *Input_output_parser::parse_input_tree(Fasta_reader *fr,vector<Fasta_entry> *sequences,bool reference_alignment,int n_threads)
{
    /***********************************************************************/
    /*  Read the guidetree                                                 */
    /***********************************************************************/

    Node *root;
    if(Settings_handle::st.is("treefile"))
    {
        string treefile =  Settings_handle::st.get("treefile").as<string>();
        Log_output::write_out("Tree file: "+treefile+"\n",1);

        Newick_reader nr;
        string tree;
        try
        {
            tree = nr.read_tree(treefile);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the guide tree file '"+treefile+"'.\nExiting.\n\n",0);
            exit(1);
        }
        try {
            root = nr.parenthesis_to_tree(tree);
        }
        catch(exception e)
        {
            Log_output::write_out("The guide tree should be a rooted binary tree. Trying mid-point rooting.\n\n",0);

            Tree_node tn;
            tree = tn.get_rooted_tree(tree);
            root = nr.parenthesis_to_tree(tree);
        }

        root = nr.parenthesis_to_tree(tree);
    }
    else if(Settings_handle::st.is("ref-treefile"))
    {
        string treefile =  Settings_handle::st.get("ref-treefile").as<string>();
        Log_output::write_out("Reference tree file: "+treefile+"\n",1);

        Newick_reader nr;
        string tree;
        try
        {
            tree = nr.read_tree(treefile);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the reference tree file '"+treefile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        try {
            root = nr.parenthesis_to_tree(tree);
        }
        catch(exception e)
        {
            Log_output::write_out("The reference guide tree should be a rooted binary tree. \nExiting.\n\n",0);
            exit(1);
        }
    }
    else if(reference_alignment && sequences->size()==1)
    {
        int data_type = fr->check_sequence_data_type(sequences);

        root = new Node();
        root->set_name(sequences->at(0).name);
        root->add_name_comment( sequences->at(0).comment );
        root->set_distance_to_parent(0);

        root->add_sequence( sequences->at(0), data_type);
    }
    else
    {
        if(Settings_handle::st.is("raxml-tree"))
        {
            RAxML_tree rt;
            Mafft_alignment ma;

            int data_type = fr->check_sequence_data_type(sequences);

            if(reference_alignment && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing RAxML guidetree from the given input alignment.",0);
                string tree = rt.infer_phylogeny(sequences,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else if(ma.test_executable() && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing an initial alignment with MAFFT.",0);
                ma.align_sequences(input);

                Log_output::write_msg("Computing RAxML guidetree from the initial alignment.",0);
                string tree = rt.infer_phylogeny(input,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else
            {
                Log_output::write_out("\nError: No tree file defined. MAFFT and/or RAxML not found.\n",0);
                Settings_handle::st.info();

                exit(1);
            }
        }
        else if(Settings_handle::st.is("bppdist-tree"))
        {
            BppDist_tree rt;
            Mafft_alignment ma;

            int data_type = fr->check_sequence_data_type(sequences);

            if(reference_alignment && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing BppDist guidetree from the given input alignment.",0);
                string tree = rt.infer_phylogeny(sequences,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else if(ma.test_executable() && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing an initial alignment with MAFFT.",0);
                ma.align_sequences(input);

                Log_output::write_msg("Computing BppDist guidetree from the initial alignment.",0);
                string tree = rt.infer_phylogeny(input,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else
            {
                Log_output::write_out("\nError: No tree file defined. MAFFT and/or BppDist not found.\n",0);
                Settings_handle::st.info();

                exit(1);
            }
        }
        else
        {
            FastTree_tree rt;
            Mafft_alignment ma;

            int data_type = fr->check_sequence_data_type(sequences);

            if(reference_alignment && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing FastTree guidetree from the given input alignment.",0);
                string tree = rt.infer_phylogeny(sequences,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else if(ma.test_executable() && rt.test_executable())
            {
                bool is_protein = (data_type==Model_factory::protein);

                vector<Fasta_entry> *input = sequences;
                vector<Fasta_entry> translated;

                if(Settings_handle::st.is("codons") && data_type==Model_factory::dna)
                {
                    this->translate_codons(sequences,&translated);
                    input = &translated;
                    is_protein = true;
                }

                Log_output::write_msg("Computing an initial alignment with MAFFT.",0);
                ma.align_sequences(input);

                Log_output::write_msg("Computing FastTree guidetree from the initial alignment.",0);
                string tree = rt.infer_phylogeny(input,is_protein,n_threads);

                Tree_node tn;
                tree = tn.get_rooted_tree(tree);

                Newick_reader nr;
                root = nr.parenthesis_to_tree(tree);
            }
            else
            {
                Log_output::write_out("\nError: No tree file defined. MAFFT and/or FastTree not found.\n",0);
                Settings_handle::st.info();

                exit(1);
            }
        }

        string outfile =  "outfile";
        if(Settings_handle::st.is("outfile"))
            outfile =  Settings_handle::st.get("outfile").as<string>();
        outfile += ".tre";

        ofstream output( outfile.c_str(), (ios::out));
        output<<root->print_tree()<<endl;
        output.close();

    }

    return root;
}

/************************************************************************************/

void Input_output_parser::match_sequences_and_tree(Fasta_reader *fr, std::vector<Fasta_entry> *sequences,Node *root,bool reference_alignment,int *data_type)
{
    /***********************************************************************/
    /*  Check that input is fine and place the sequences to nodes          */
    /***********************************************************************/

    // Check that the guidetree and sequences match

    vector<Node*> leaf_nodes;
    root->get_leaf_nodes(&leaf_nodes);


    bool tree_branches_ok = fr->check_sequence_names(sequences,&leaf_nodes);
    if(!tree_branches_ok)
    {
        Log_output::write_out("Attempting to prune the tree: ",2);
        Log_output::flush();

        root->prune_tree();
        Log_output::write_out("pruning done.\n",2);
        Log_output::write_out("New tree:\n"+root->print_tree()+"\n\n",3);
    }

    leaf_nodes.clear();
    root->get_leaf_nodes(&leaf_nodes);


    // Check the data type

    *data_type = fr->check_sequence_data_type(sequences);

    if(*data_type == Model_factory::protein && Settings_handle::st.is("use-consensus"))
    {
        Log_output::write_out("Option '--use-consensus' not supported for proteins. Exiting.\n\n",0);
        exit(0);
    }

    if(!fr->check_alphabet(sequences,*data_type))
        Log_output::write_out(" Warning: Illegal characters in input sequences removed!\n",2);


    //  Place the sequences to nodes
    fr->place_sequences_to_nodes(sequences,&leaf_nodes,reference_alignment,*data_type);

    if( ! ( Settings_handle::st.is("no-anchors") || Settings_handle::st.is("use-prefix-anchors") ) )
    {
        Exonerate_queries er;
        if(!er.test_executable())
            Log_output::write_out("The executable for Exonerate not found! Alignment anchoring not used!\n",0);
    }
}

/************************************************************************************/

void Input_output_parser::define_alignment_model(Fasta_reader *fr,Model_factory *mf,int data_type)
{
    /***********************************************************************/
    /*  Define the alignment model                                         */
    /***********************************************************************/

    if( ( data_type==Model_factory::dna && Settings_handle::st.is("codons") ) || data_type==Model_factory::codon)
    {
        // Create a codon alignment model using KHG.
        Log_output::write_out("Model_factory: creating a codon model\n",3);
        mf->codon_model(&Settings_handle::st); // does it need the handle????
    }
    else if(data_type==Model_factory::dna)
    {
        // Create a DNA alignment model using empirical base frequencies.
        Log_output::write_out("Model_factory: creating a DNA model\n",3);
        float *dna_pi = fr->base_frequencies();
        mf->dna_model(dna_pi,&Settings_handle::st);
    }
    else if(data_type==Model_factory::protein)
    {
        // Create a protein alignment model using WAG.
        Log_output::write_out("Model_factory: creating a protein model\n",3);
        mf->protein_model(&Settings_handle::st); // does it need the handle????
    }

}

/************************************************************************************/

void Input_output_parser::output_aligned_sequences(Fasta_reader *fr,std::vector<Fasta_entry> *sequences,Node *root)
{

    /***********************************************************************/
    /*  Collect the results and output them                                */
    /***********************************************************************/

    fr->set_chars_by_line(70);

//    cout<<"root: "<<root->get_name()<<endl;
//    root->show_seqs();



    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,true);

    Log_output::clean_output();

    // See if any sequences were placed
    //
    if(Settings_handle::st.is("queryfile") && root->get_number_of_query_leaves()==0)
    {
        Log_output::write_out("Failed to extend the alignment. No output created.\n",0);
    }
    else
    {

        // Save results in output file
        //
        string outfile =  "outfile";

        if(Settings_handle::st.is("outfile"))
            outfile =  Settings_handle::st.get("outfile").as<string>();

        string format = "fasta";
        if(Settings_handle::st.is("outformat"))
            format = Settings_handle::st.get("outformat").as<string>();

        if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
            Log_output::write_out("Alignment files: "+outfile+fr->get_format_suffix(format)+", "+outfile+".xml\n",0);
        else
            Log_output::write_out("Alignment file: "+outfile+fr->get_format_suffix(format)+"\n",0);

        if(!Settings_handle::st.is("treefile") && !Settings_handle::st.is("ref-treefile"))
            Log_output::write_out("Guidetree file: "+outfile+".tre\n",0);


        bool do_ancestors = Settings_handle::st.is("events") || Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors")
                || Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx");

        BppAncestors bppa;
        bool infer_bppa_ancestors = ( bppa.test_executable() && not Settings_handle::st.is("no-bppancestors") &&
                                      ( do_ancestors ) &&
                                      (int)sequences->size()<500 );

        if(infer_bppa_ancestors)
        {
            bool worked = bppa.infer_ancestors(root,&aligned_sequences);

            bppa.count_events(root,&aligned_sequences,outfile);

            if(Settings_handle::st.is("events"))
            {
                if( Settings_handle::st.is("translate") || Settings_handle::st.is("mt-translate")
                    || Settings_handle::st.is("find-orfs") )
                    Log_output::write_out("Inferred evolutionary events: "+outfile+".events, "+outfile+".codon.events\n",0);
                else
                    Log_output::write_out("Inferred evolutionary events: "+outfile+".events\n",0);
            }
        }
        else if( do_ancestors )
        {
            Log_output::write_out("\nWarning: BppAncestors not used. Performing approximate ancestor reconstruction.\n\n",0);
        }


        // Write alignment as flatfile (fasta by default)
        //
        if(Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors"))
        {
            fr->write(outfile, aligned_sequences, format, true);
        }
        else
        {
            vector<Fasta_entry> leaf_sequences;
            set<string> names;
            root->get_leaf_node_names(&names);
            vector<Fasta_entry>::iterator si = aligned_sequences.begin();
            for(;si!=aligned_sequences.end();si++)
            {
                if(names.find(si->name) != names.end())
                    leaf_sequences.push_back(*si);
            }
            fr->write(outfile, leaf_sequences, format, true);
        }

        // Write alignment as HSAML
        //
        int count = 1;
        root->set_name_ids(&count);

        if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
        {
            Xml_writer xw;
            xw.write(outfile, root, aligned_sequences, true);
        }

        if( Settings_handle::st.is("translate") || Settings_handle::st.is("mt-translate")
            || Settings_handle::st.is("find-orfs") )
        {
            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append(".codon");
            if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
                Log_output::write_out("Back-translated alignment files: "+outfile+fr->get_format_suffix(format)+", "+outfile+".xml\n",0);
            else
                Log_output::write_out("Back-translated alignment file: "+outfile+fr->get_format_suffix(format)+"\n",0);


            map<string,string> dna_seqs;
            fr->get_DNA_seqs(root, sequences, &dna_seqs);


            vector<Fasta_entry> dna_sequences;
            bool ancestors_done = false;

            if( not do_ancestors  )
            {
                vector<Fasta_entry> leaf_sequences;
                set<string> names;
                root->get_leaf_node_names(&names);
                vector<Fasta_entry>::iterator si = aligned_sequences.begin();
                for(;si!=aligned_sequences.end();si++)
                {
                    if(names.find(si->name) != names.end())
                    {
                        leaf_sequences.push_back(*si);
                    }
                }
                fr->backtranslate_dna(leaf_sequences,&dna_seqs,dna_sequences,infer_bppa_ancestors);
            }
            else if(infer_bppa_ancestors )
            {
                fr->backtranslate_dna(aligned_sequences,&dna_seqs,dna_sequences,infer_bppa_ancestors);

                ancestors_done = bppa.infer_ancestors(root,&dna_sequences,true);

                if( ancestors_done )
                {
                    bppa.count_events(root,&dna_sequences,outfile,true);

                    if( aligned_sequences.size() == dna_sequences.size() )
                        fr->write(outfile, dna_sequences, format, true);
                }
            }

            if( do_ancestors && not ancestors_done)
            {
                Newick_reader nr;
                Node *croot;
                croot = nr.parenthesis_to_tree(root->print_nhx_tree_with_intIDs());

//                cout<<root->print_nhx_tree_with_intIDs()<<endl;
//                cout<<croot->print_nhx_tree_with_intIDs()<<endl;

                vector<Fasta_entry> leaf_sequences;
                set<string> names;
                croot->get_leaf_node_names(&names);

                vector<Fasta_entry>::iterator si = aligned_sequences.begin();
                for(;si!=aligned_sequences.end();si++)
                {
                    if(names.find(si->name) != names.end())
                    {
                        leaf_sequences.push_back(*si);
                    }
                }

                dna_sequences.clear();
                fr->backtranslate_dna(leaf_sequences,&dna_seqs,dna_sequences,false);


                vector<Node*> leaf_nodes;
                croot->get_leaf_nodes(&leaf_nodes);

                fr->place_sequences_to_nodes(&dna_sequences,&leaf_nodes,true,int(Model_factory::codon));

                Model_factory mf(int(Model_factory::codon));
                mf.codon_model(&Settings_handle::st);

                croot->read_reference_alignment(&mf,true);

                dna_sequences.clear();
                croot->get_alignment(&dna_sequences,true);


                if( aligned_sequences.size() == dna_sequences.size() )
                    fr->write(outfile, dna_sequences, format, true);

                Log_output::write_out("\nReconstructing ML ancestral codon sequences failed. Outputting parsimony ancestors.\n\n",0);

                if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
                {
                    int count = 1;
                    croot->set_name_ids(&count);

                    Xml_writer xw;
                    xw.write(outfile, croot, dna_sequences, true);
                }

                delete croot;
            }
            else
            {
                if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
                {
                    Xml_writer xw;
                    xw.write(outfile, root, dna_sequences, true);
                }
            }

            if( not ancestors_done )
            {
                vector<Fasta_entry> leaf_sequences;
                set<string> names;
                root->get_leaf_node_names(&names);
                vector<Fasta_entry>::iterator si = dna_sequences.begin();
                for(;si!=dna_sequences.end();si++)
                {
                    if(names.find(si->name) != names.end())
                        leaf_sequences.push_back(*si);
                }
                fr->write(outfile, leaf_sequences, format, true);
            }


            /***************************************************/

            if(Settings_handle::st.is("build-contigs"))
            {
                string outfile =  "outfile";
                if(Settings_handle::st.is("outfile"))
                    outfile =  Settings_handle::st.get("outfile").as<string>();

                outfile.append("_contigs.dna");
                Log_output::write_out("Back-translated contig file: "+outfile+".fas\n",1);

                fr->set_chars_by_line(70);
                fr->write_dna(outfile, aligned_sequences, *sequences,root,true,Fasta_reader::contig_alignment);
            }
            if(Settings_handle::st.is("output-consensus"))
            {
                string outfile =  "outfile";
                if(Settings_handle::st.is("outfile"))
                    outfile =  Settings_handle::st.get("outfile").as<string>();

                outfile.append("_consensus.dna");
                Log_output::write_out("Back-translated consensus file: "+outfile+".fas\n",1);

                fr->set_chars_by_line(70);
                fr->write_dna(outfile, aligned_sequences, *sequences,root,true,Fasta_reader::consensus_only);
            }

            /***************************************************/

        }

        if(Settings_handle::st.is("prune-extended-alignment"))
            this->prune_extended_alignment(fr,root,&aligned_sequences);

        //////////////////

        if(Settings_handle::st.is("trim-extended-alignment"))
        {
            set<string> readnames;
            root->get_read_node_names(&readnames);

            int last_site = 0;
            int first_site = aligned_sequences.at(0).sequence.length();

            vector<Fasta_entry> trimmed_sequences;

            vector<Fasta_entry>::iterator fit = aligned_sequences.begin();
            for(;fit!=aligned_sequences.end();fit++)
            {
                if(readnames.find(fit->name)!=readnames.end())
                {
                    for(int i=0;i<(int)fit->sequence.length();i++)
                    {
                        if(fit->sequence.at(i)!='-' && i < first_site)
                            first_site = i;
                        if(fit->sequence.at(i)!='-' && i > last_site)
                            last_site = i;
                    }

                }
            }

            first_site = max(first_site-Settings_handle::st.get("trim-keep-sites").as<int>(),0);
            last_site  = min(last_site+Settings_handle::st.get("trim-keep-sites").as<int>(),(int)aligned_sequences.at(0).sequence.length());

            fit = aligned_sequences.begin();
            for(;fit!=aligned_sequences.end();fit++)
            {
                Fasta_entry seq;
                seq.name = fit->name;
                seq.comment = fit->comment;
                seq.sequence = fit->sequence.substr(first_site,last_site-first_site);
                trimmed_sequences.push_back(seq);
            }


            this->output_pruned_alignment(fr,root,root,&trimmed_sequences,"Trimmed",".trimmed");

        }

        //////////////////

        if(Settings_handle::st.is("output-ancestors") || Settings_handle::st.is("ancestors"))
        {
            fr->write_anctree(outfile, root);
        }

        if(Settings_handle::st.is("output-nhx-tree") || Settings_handle::st.is("guidetree"))
        {
            root->write_nhx_tree(outfile,"nhx_tree");
        }

        if( Settings_handle::st.is("scale-branches") ||
             Settings_handle::st.is("truncate-branches") ||
              Settings_handle::st.is("fixed-branches") )
        {
            Log_output::write_out("Modified guide tree: " +root->print_tree()+"\n",2);
        }


        /***************************************************/

        if(Settings_handle::st.is("build-contigs"))
        {
            vector<Fasta_entry> contigs;
            root->reconstruct_contigs(&contigs,false);

            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append("_contigs");
            Log_output::write_out("Contig file: "+outfile+".fas\n",1);

            fr->set_chars_by_line(70);
            fr->write(outfile, contigs, "fasta", true);
        }

        if(Settings_handle::st.is("output-consensus"))
        {
            vector<Fasta_entry> contigs;
            root->reconstruct_contigs(&contigs,false,true);

            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append("_consensus");
            Log_output::write_out("Consensus file: "+outfile+".fas\n",1);

            fr->remove_gap_only_columns(&contigs);

            fr->set_chars_by_line(70);
            fr->write(outfile, contigs, "fasta", true);
        }

        /***************************************************/

        if(Settings_handle::st.is("output-graph"))
        {
            fr->write_graph(outfile, root, true);
        }

        if(Settings_handle::st.is("mpost-graph-file")){
            root->write_sequence_graphs();
        }

        /***************************************************/
    }

}


void Input_output_parser::prune_extended_alignment(Fasta_reader *fr,Node *root,vector<Fasta_entry> *aligned_sequences)
{
    bool is_protein = root->get_sequence()->get_data_type()==Model_factory::protein;


    if(Settings_handle::st.is("prune-keep-number"))
    {

        Newick_reader nr;
        string tree = root->print_tree();
        Node *tmp_root = nr.parenthesis_to_tree(tree);

        set<string> removenames;
        root->get_all_terminal_node_names(&removenames);

        set<string> readnames;
        root->get_read_node_names(&readnames);

        if(Settings_handle::st.get("prune-keep-number").as<int>()>1)
        {
            BppPhySamp_tree bppphys;
            if(bppphys.test_executable())
            {
                removenames.clear();
                bppphys.reduce_sequences(&removenames,is_protein);

                tmp_root->set_has_sequence();
                tmp_root->unset_has_sequence(&removenames);
                tmp_root->set_has_sequence(&readnames);
                tmp_root->prune_tree();

                this->output_pruned_alignment(fr,root,tmp_root,aligned_sequences,"Pruned",".pruned");
            }
            else
            {
                Log_output::write_out("The executable for BppPhySamp not found! The alignment cannot be pruned!",0);
            }
        }
        else
        {
            if((int)readnames.size()>1)
            {
                tmp_root->unset_has_sequence();
                tmp_root->set_has_sequence(&readnames);
                tmp_root->prune_tree();

                this->output_pruned_alignment(fr,root,tmp_root,aligned_sequences,"Pruned",".pruned");
            }
            else
            {
                Log_output::write_out("Only one query sequence: pruned alignment without reference not meaningful.\n",0);
            }
        }
        delete tmp_root;
    }

    if(Settings_handle::st.is("prune-keep-closest"))
    {
        Newick_reader nr;
        string tree = root->print_tree();
        Node *tmp_root = nr.parenthesis_to_tree(tree);

        set<string> removenames;
        root->get_all_terminal_node_names(&removenames);

        set<string> readnames;
        root->get_read_node_names(&readnames);

        set<string> keepnames;
        root->get_closest_reference_leaves(&keepnames);

        tmp_root->unset_has_sequence();

        tmp_root->set_has_sequence(&keepnames);
        tmp_root->set_has_sequence(&readnames);
        tmp_root->prune_tree();

        set<string> leaves_kept;
        tmp_root->get_leaf_node_names(&leaves_kept);

        this->output_pruned_alignment(fr,root,tmp_root,aligned_sequences,"Pruned-with-closest",".pruned_closest");

        delete tmp_root;
    }

}

void Input_output_parser::output_pruned_alignment(Fasta_reader *fr,Node *root,Node *tmp_root,vector<Fasta_entry> *aligned_sequences,string desc,string prefix)
{

    set<string> readnames;
    root->get_read_node_names(&readnames);


    set<string> leaves_kept;
    tmp_root->get_leaf_node_names(&leaves_kept);

    string outfile =  "outfile";

    if(Settings_handle::st.is("outfile"))
        outfile =  Settings_handle::st.get("outfile").as<string>();

    outfile.append(prefix);

    string format = "fasta";
    if(Settings_handle::st.is("outformat"))
        format = Settings_handle::st.get("outformat").as<string>();

    if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
        Log_output::write_out(desc+" alignment files: "+outfile+fr->get_format_suffix(format)+", "+outfile+".xml\n",0);
    else
        Log_output::write_out(desc+" alignment file: "+outfile+fr->get_format_suffix(format)+"\n",0);
    Log_output::write_out(desc+" tree file: "+outfile+".tre\n",0);


    vector<Fasta_entry> pruned_sequences;
    vector<Fasta_entry>::iterator fit = aligned_sequences->begin();
    for(;fit!=aligned_sequences->end();fit++)
    {
        if(leaves_kept.find(fit->name)!=leaves_kept.end())
        {
            pruned_sequences.push_back(*fit);
        }
    }
    fr->remove_gap_only_columns(&pruned_sequences);

    if(Settings_handle::st.is("trim-extended-alignment"))
    {
        int last_site = 0;
        int first_site = pruned_sequences.at(0).sequence.length();

        vector<Fasta_entry>::iterator fit = pruned_sequences.begin();
        for(;fit!=pruned_sequences.end();fit++)
        {
            if(readnames.find(fit->name)!=readnames.end())
            {
                for(int i=0;i<(int)fit->sequence.length();i++)
                {
                    if(fit->sequence.at(i)!='-' && i < first_site)
                        first_site = i;
                    if(fit->sequence.at(i)!='-' && i > last_site)
                        last_site = i;
                }

            }
        }

        first_site = max(first_site-Settings_handle::st.get("trim-keep-sites").as<int>(),0);
        last_site  = min(last_site+Settings_handle::st.get("trim-keep-sites").as<int>(),(int)pruned_sequences.at(0).sequence.length());

        fit = pruned_sequences.begin();
        for(;fit!=pruned_sequences.end();fit++)
        {
            fit->sequence = fit->sequence.substr(first_site,last_site-first_site);
        }
    }

    fr->set_chars_by_line(70);
    fr->write(outfile, pruned_sequences, format, true);

    tmp_root->write_nhx_tree(outfile,"tre");

    if(Settings_handle::st.is("xml") || Settings_handle::st.is("xml-nhx"))
    {
        vector<Node*> leaf_nodes;
        tmp_root->get_leaf_nodes(&leaf_nodes);

        bool tree_branches_ok = fr->check_sequence_names(&pruned_sequences,&leaf_nodes);
        if(!tree_branches_ok)
        {
            Log_output::write_out("Attempting to prune the tree: ",2);
            Log_output::flush();

            root->prune_tree();
            Log_output::write_out("pruning done.\n",2);
            Log_output::write_out("New tree:\n"+root->print_tree()+"\n\n",3);
        }

        leaf_nodes.clear();
        tmp_root->get_leaf_nodes(&leaf_nodes);

        //  Place the sequences to nodes
        fr->place_sequences_to_nodes(&pruned_sequences,&leaf_nodes,true);
        int count = 1;
        tmp_root->name_internal_nodes(&count);

        count = 1;
        tmp_root->set_name_ids(&count);

        Xml_writer xw;
        xw.write(outfile, tmp_root, pruned_sequences, true, true);
    }
}
