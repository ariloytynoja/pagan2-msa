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

#ifndef XML_WRITER_H
#define XML_WRITER_H

#include <string>
#include <vector>
#include <fstream>
#include "utils/exceptions.h"
#include "main/node.h"
#include "utils/settings.h"
#include "utils/fasta_entry.h"

using namespace std;

namespace ppa
{

class Xml_writer
{
public:
    Xml_writer();

    void write(ostream & output, const Node *root, const vector<Fasta_entry> & seqs, bool append_comment = false) const ;
    void write(const string & path, const Node *root, const vector<Fasta_entry> & seqs, bool append_comment = false, bool overwrite=true) const 
    {
        ofstream output( (path+".xml").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        write(output, root, seqs, append_comment);
        output.close();
    }

};

}
#endif // XML_WRITER_H
