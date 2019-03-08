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

#ifndef CODON_TRANSLATION_H
#define CODON_TRANSLATION_H

#include <string>
#include <map>

namespace ppa {

class Codon_translation
{
    std::map<std::string,std::string> codon_to_aa;

public:
    Codon_translation();
    void define_translation_tables();
    std::string gapped_DNA_to_protein(std::string *sequence) const;
    std::string codon_to_amino(std::string codon) const
    {
        if (codon_to_aa.find(codon) == codon_to_aa.end())
        {
            return "X";
        }
        else
        {
            return codon_to_aa.find(codon)->second;
        }
    }
};
}

#endif // CODON_TRANSLATION_H
