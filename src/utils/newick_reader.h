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
// File: TreeTemplateTools.h
// Created by:  Julien Dutheil
// Created on: Fri Oct  13 13:00 2006
// From file TreeTools.h
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

#ifndef NEWICK_READER_H
#define NEWICK_READER_H

#include <string>
#include "utils/exceptions.h"
#include "main/node.h"

using namespace std;

namespace ppa
{

class Newick_reader
{
    int node_index;
    bool has_warned = false;
public:
    Newick_reader();
    ~Newick_reader() {}

    string read_tree(const string & filename) ;
    Node * parenthesis_to_node(const string & description) ;
    Node * parenthesis_to_tree(const string & description) ;

    struct Element
    {
      string content;
      string length;
      string nhx;
      string nodeid;
    };

    bool removed_multifurcation() { return has_warned; }
    static Element get_element(const string & elt) ;
};

} //end of namespace ppa.

#endif // NEWICK_READER_H
