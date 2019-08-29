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
// File: TreeTemplateTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct  13 13:00 2006
// From file TreeTools.cpp
// Created on: Wed Aug  6 13:45:28 2003
//

/*
Copyright or <A9> or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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


#include <iostream>
#include <fstream>
#include <vector>

#include "utils/newick_reader.h"
#include "utils/text_utils.h"

using namespace std;

using namespace ppa;

Newick_reader::Newick_reader()
{
    node_index = 1;
}

/*******************************************/

string Newick_reader::read_tree(const string & filename) throw (IOException)
{
    ifstream ist(filename.c_str());
    if (!ist) throw IOException("Newick_reader::read_tree. Can't open input file: "+filename);

    string tree;
    string row;
    while (ist>>row)
    {
        tree += Text_utils::remove_whitespaces(row);
    }

    return tree;
}

/*******************************************/

Node * Newick_reader::parenthesis_to_node(const string & description) throw (Exception)
{
//    cout<<description<<endl;
    Newick_reader::Element elt = Newick_reader::get_element(description);

    //New node:
    Node * node = new Node();
    if(!Text_utils::is_empty(elt.length))
        node->set_distance_to_parent(Text_utils::to_double(elt.length));
    else
        node->set_distance_to_parent(0);

    stringstream nhx_tag("");
    if(!Text_utils::is_empty(elt.nhx))
    {
        String_tokenizer * st = new String_tokenizer(elt.nhx, ":", true, false);

        node->set_nhx_tid("");
        while (st->has_more_token())
        {
            string block = st->next_token();
            block = Text_utils::remove_surrounding_whitespaces(block);
            if(block.substr(0,4)=="TID=")
            {
                block = block.substr(4);
                node->set_nhx_tid(block);
            }
            else
            {
                if(nhx_tag.str().length()==0)
                    nhx_tag<<block;
                else
                    nhx_tag<<":"<<block;
            }
        }
        delete st;
//        node->set_distance_to_parent(Text_utils::to_double(elt.length));
    }
    node->set_nhx_tag(nhx_tag.str());

    Node_tokenizer nt(elt.content);
    vector<string> elements;
    while(nt.has_next())
    {
        elements.push_back(nt.next());
    }

    if(elements.size() == 1)
    {
        //This is a leaf:
        string name = Text_utils::remove_surrounding_whitespaces(elements[0]);
        node->set_name(name);

        node->is_leaf(true);
    }
    else if(elements.size()!=2)
    {
        string t("(");
        int i=0;
        for(;i<int(elements.size())-2;i++)
            t+=elements[i]+",";
        t+=elements[i]+")";

        //This is a node:
        try
        {
            Node * child_0 = parenthesis_to_node(t);
            node->add_left_child(child_0);

            Node * child_1 = parenthesis_to_node(elements[elements.size()-1]);
            node->add_right_child(child_1);

//            node->set_name(description);
            if(elt.nodeid!="")
            {
                node->set_name(elt.nodeid);
            }
            else
            {
                stringstream ss("");
                ss<<"node"<<node_index;
                node->set_name(ss.str());
            }
            node_index++;
            node->is_leaf(false);
        }
        catch(exception e)
        {
            throw e;
        }

        if(! this->has_warned)
        {
            Log_output::write_out("Warning: removing multifurcation\n",1);
            this->has_warned = true;
        }
    }

    else
    {

        //This is a node: 
        try
        {
            Node * child_0 = parenthesis_to_node(elements[0]);
            node->add_left_child(child_0);

            Node * child_1 = parenthesis_to_node(elements[1]);
            node->add_right_child(child_1);

//            node->set_name(description);
            if(elt.nodeid!="")
            {
                node->set_name(elt.nodeid);
            }
            else
            {
                stringstream ss("");
                ss<<"node"<<node_index;
                node->set_name(ss.str());
            }
            node_index++;
            node->is_leaf(false);
        }
        catch(exception e)
        {
            throw e;
        }
    }
    return node;
}

/*******************************************/

Node * Newick_reader::parenthesis_to_tree(const string & description) throw (Exception)
{
    string::size_type lastP  = description.rfind(')');
    if(lastP == string::npos)
        throw Exception("Newick_reader::parenthesis_to_tree(). Bad format: no closing parenthesis found.");
    string::size_type firstP = description.find('(');
    if(firstP == string::npos)
        throw Exception("Newick_reader::parenthesis_to_tree(). Bad format: no opening parenthesis found.");
    string::size_type semi = description.rfind(';');
    if(semi == string::npos)
        throw Exception("Newick_reader::parenthesis_to_tree(). Bad format: no semi-colon found.");
    if(lastP <= firstP)
        throw Exception("Newick_reader::parenthesis_to_tree(). Bad format: closing parenthesis before opening parenthesis.");

    string content = description.substr(firstP + 1, lastP - firstP - 1);
    string element = semi == string::npos ? description.substr(lastP + 1) : description.substr(lastP + 1, semi - lastP - 1);

    //New root node:
    Node * node = new Node();


    Newick_reader::Element elt = Newick_reader::get_element(description);

    if(!Text_utils::is_empty(elt.length))
        node->set_distance_to_parent(Text_utils::to_double(elt.length));
    else
        node->set_distance_to_parent(0);


    if(!Text_utils::is_empty(elt.nhx))
    {
        String_tokenizer * st = new String_tokenizer(elt.nhx, ":", true, false);

        node->set_nhx_tid("");
        while (st->has_more_token())
        {
            string block = st->next_token();
            block = Text_utils::remove_surrounding_whitespaces(block);
            if(block.substr(0,4)=="TID=")
            {
                block = block.substr(4);
                node->set_nhx_tid(block);
            }
        }

        delete st;
//        node->set_distance_to_parent(Text_utils::to_double(elt.length));
    }

    Node_tokenizer nt(content);
    vector<string> elements;
    while(nt.has_next())
    {
        elements.push_back(nt.next());
    }

    if(elements.size()!=2)
    {
        // Alignment requires a binary guidetree!
//        throw Exception("Newick_reader::parenthesis_to_tree(). Not a binary tree: " + content);
        throw Exception("The guidetree should be a rooted binary tree. Exiting.\n");
    }

    else
    {
        //This is a node:
        try
        {
            Node * child_0 = parenthesis_to_node(elements[0]);
            node->add_left_child(child_0);

            Node * child_1 = parenthesis_to_node(elements[1]);
            node->add_right_child(child_1);

            if(elt.nodeid!="")
            {
                node->set_name(elt.nodeid);
            }
            else
            {
                node->set_name("root");
            }
        }
        catch(exception e)
        {
            throw e;
        }

    }
    return node;
}

/*******************************************/

Newick_reader::Element Newick_reader::get_element(const string & elt) throw (Exception)
{
    Element element;
    element.length    = ""; //default
    element.nhx = "";
//cout<<"elt: "<<elt<<endl;
//    string::size_type colon = elt.rfind(':');
//    string::size_type colon = elt.find(':');
    try
    {
        string::size_type openNHX = elt.rfind("[&&NHX");
        string::size_type closeNHX = elt.rfind(']');
        string::size_type lastBracket = elt.rfind(')');

        string eltt = elt;
        if( (openNHX != string::npos && lastBracket != string::npos && lastBracket < openNHX)
            || (openNHX != string::npos && lastBracket == string::npos) )
        {
            if(closeNHX != string::npos)
            {
                element.nhx = elt.substr(openNHX+1,closeNHX-openNHX-1);
                eltt = elt.substr(0,openNHX);
            }
        }

        string::size_type colon = eltt.rfind(':');
        string::size_type endP = eltt.rfind(')');

        string elt2;
        if(colon != string::npos)
        {
            if(endP == string::npos || colon>endP){
                //this is an element with length:
                elt2 = eltt.substr(0, colon);
                element.length = eltt.substr(colon + 1);
            }
            else
            {
                //this is an element without length;
                elt2 = eltt;
            }
        }
        else
        {
            //this is an element without length;
            elt2 = eltt;
        }

//        cout<<"elt2: "<<elt2<<endl;

        string::size_type lastP = elt2.rfind(')');
        string::size_type firstP = elt2.find('(');
        if(firstP == string::npos)
        {
            //This is a leaf:
            element.content = elt2;
        }
        else
        {
            //This is a node:
            if(lastP < firstP)
                throw Exception("Newick_reader::get_element. Invalid format: bad closing parenthesis in " + elt2);

            element.content = elt2.substr(firstP + 1, lastP - firstP - 1);

            string::size_type firstH = elt2.find('#',lastP);
            string::size_type lastH = elt2.find('#',firstH+1);

            if(firstH != string::npos)
            {
                int num = Text_utils::to_int(elt2.substr(firstH+1,lastH-firstH));
                if(num>0)
                {
                    element.nodeid = elt2.substr(firstH,lastH-firstH+1);
                }
            }
        }
    }
    catch(exception e)
    {
        throw Exception("Newick_reader::get_element. Bad tree description: " + elt);
    }
    return element;
}

/*******************************************/

