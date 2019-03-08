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

#ifndef INT_MATRIX_H
#define INT_MATRIX_H

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#include <string>
#include <iostream>
#include <cassert>
#include "utils/log_output.h"

namespace ppa {

class Int_matrix{
private:
    int x;
    int y;
    int z;
    int w;
    bool xar;
    bool yar;
    bool zar;
    bool war;

    std::string name;
    int* data;
    int i,j,k,l;
public:
    Int_matrix(int x, std::string name="");
    Int_matrix(int x, int y, std::string name="");
    Int_matrix(int x, int y, int z, std::string name="");
    Int_matrix(int x, int y, int z, int w, std::string name="");

    ~Int_matrix();

    void allocate();
    void initialise(int v = 0);

    int g(int xa, int ya=0, int za = 0, int wa = 0) {  assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z); assert(wa>=0); assert(wa<w); return data[xa + ya*x + za*x*y + wa*x*y*z]; }
    void s(int v, int xa, int ya=0, int za = 0, int wa = 0);
    void a(int v, int xa, int ya=0, int za = 0, int wa = 0) { assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z);  assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] += v; }
    void printName() {
        std::stringstream ss;
        ss<<"Name "<<name<<": x = "<<x<<", y = "<<y<<", z = "<<z<<", w = "<<w<<std::endl;
        Log_output::write_out(ss.str(),0);
    }
    void print();
    void print(int m);

    int X() { return x; }

    void resize(int i);
    void copyData(int *tmp,int new_x,int new_y,int new_z,int new_w);

    void allowResize(bool xr, bool yr=false, bool zr=false, bool wr=false);
};

}
#endif
