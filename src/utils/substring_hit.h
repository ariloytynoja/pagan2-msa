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

#ifndef SUBSTRING_HIT_H
#define SUBSTRING_HIT_H

#include <iostream>
#include <sstream>

namespace ppa
{

using namespace std;

struct Substring_hit
{
    int start_site_1;
    int start_site_2;
    bool plus_strand_1;
    bool plus_strand_2;
    int length;
    int score;
    Substring_hit() : plus_strand_1(true), plus_strand_2(true) {

    }
    bool operator < (const Substring_hit& hit) const{
        return (score < hit.score);
    }
    void print(stringstream& s){
        s << "Strands: (" << (plus_strand_1 ? "plus" : "minus");
        s << ", " << (plus_strand_2 ? "plus" : "minus");
        s << "), starts: (" << start_site_1;
        s << ", " << start_site_2;
        s << "), length: " << length;
        s << ", score: " << score << endl;
    }
};

}
#endif // SUBSTRING_HIT_H
