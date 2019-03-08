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

#ifndef LOG_OUTPUT_H
#define LOG_OUTPUT_H

#include <ostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

namespace ppa
{

class Log_output
{
    static ofstream fs;
    static ostream *os;
    static bool newline;
    static int header_length;
    static int msg_length;
    static int msg2_length;
    static string prev_header;
    static string prev_msg;
    static string prev_msg2;
public:
    Log_output();
    static void open_stream();
    static void write_out(const string str,const int priority);
    static void write_msg(const string str,const int priority);
    static void append_msg(const string str,const int priority);
    static void write_header(const string str,const int priority);
    static void write_new_header(const string str,const int priority);
    static void write_warning(const string str,const int priority);
    static void flush() {os->flush();}
    static void write_out(const string str,const string option);
    static void clean_output();
    static string itos(int i)
    {
        stringstream s;
        s << i;
        return s.str();
    }

    static string ftos(float i)
    {
        stringstream s;
        s.precision(2);
        s << i;
        return s.str();
    }

};

}

#endif // LOG_OUTPUT_H
