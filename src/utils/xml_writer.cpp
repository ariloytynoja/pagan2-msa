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

#include <iostream>

#include "utils/xml_writer.h"

using namespace ppa;

Xml_writer::Xml_writer()
{

}

/****************************************************************************************/

void Xml_writer::write(ostream & output, const Node *root, const vector<Fasta_entry> & seqs, bool append_comment) const 
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Xml_writer::write. Failed to open file"); }


    output << "<ms_alignment>\n<newick>" << root->print_xml_tree() << "</newick>\n<nodes>\n";

    string seq, temp = "";  // Initialization

    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    char c1,c2; int i1;
    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        string id = root->get_id_for_name(vi->name);

        stringstream ss(vi->name);
        ss >> c1 >> i1 >> c2;
        if(c1 == '#' && c2 == '#' && i1 > 0)
        {
            output << "<node id=\"" << vi->name <<"\" name=\"" << vi->name << "\">\n";
            output << "  <sequence>\n    " << vi->sequence << "\n  </sequence>\n</node>\n";
        }
        else
        {
            string name = vi->name;
            if(append_comment)
                name.append(vi->comment);
            output << "<leaf id=\"" << id <<"\" name=\"" << name << "\">\n";
            output << "  <sequence>\n    " << vi->sequence << "\n  </sequence>\n</leaf>\n";
        }
    }

    output << "</nodes>\n</ms_alignment>\n";
}

/****************************************************************************************/
