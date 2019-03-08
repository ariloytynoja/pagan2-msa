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
#include <fstream>

#include "utils/log_output.h"
#include "utils/settings_handle.h"

using namespace std;
using namespace ppa;

Log_output::Log_output()
{
}

ofstream Log_output::fs;

ostream *Log_output::os = &cout;

bool Log_output::newline = false;

int Log_output::header_length = 0;
int Log_output::msg_length = 0;
int Log_output::msg2_length = 0;

string Log_output::prev_header;
string Log_output::prev_msg;
string Log_output::prev_msg2;

void Log_output::open_stream()
{
    if(Settings_handle::st.is("log-output-file"))
    {
        fs.open(Settings_handle::st.get("log-output-file").as<string>().c_str());
        os = &fs;
        newline = true;
    }
}

void Log_output::write_out(const string str,const string option)
{
    if(!Settings_handle::st.is(option))
        return;

    if(msg_length>0)
        *os<<endl;
    *os<<str;
    msg_length = 0;
    msg2_length = 0;
    header_length = 0;
}

void Log_output::write_out(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    for(int i=0;i<priority;i++)
        *os<<" ";
    *os<<str;
}

void Log_output::write_msg(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    if(newline || Settings::noise > 0)
    {
        Log_output::write_out(str+"\n",priority);
    }
    else
    {
         for(int i=0;i<msg_length+msg2_length;i++)
             *os<<'\b';
        for(int i=0;i<msg_length+msg2_length;i++)
                 *os<<' ';
        for(int i=0;i<msg_length+msg2_length;i++)
            *os<<'\b';

        *os<<str;
        os->flush();
        msg_length = str.length();
        msg2_length = 0;
        prev_msg = str;
        prev_msg2 = "";

    }
}

void Log_output::append_msg(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    if(newline || Settings::noise > 0)
    {
        Log_output::write_out(str+"\n",priority);
    }
    else
    {
         for(int i=0;i<msg2_length;i++)
             *os<<'\b';
        for(int i=0;i<msg2_length;i++)
                 *os<<' ';
        for(int i=0;i<msg2_length;i++)
            *os<<'\b';

        *os<<str;
        os->flush();
        msg2_length = str.length();
        prev_msg2 = str;

    }
}

void Log_output::write_header(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    if(newline || Settings::noise > 0)
    {
        Log_output::write_out(str+"\n",priority);
    }
    else
    {
        for(int i=0;i<msg_length+msg2_length;i++)
            *os<<"\b";

        if(msg_length+msg2_length>0)
            *os << "\e[A";

        for(int i=0;i<header_length;i++)
            *os<<"\b";
        for(int i=0;i<header_length;i++)
                 *os<<' ';
        for(int i=0;i<header_length;i++)
            *os<<'\b';

        if(msg_length+msg2_length==0)
        {
            prev_msg = " ";
            msg_length = 1;
        }

            *os<<str<<"\n"<<prev_msg<<prev_msg2;


        header_length = str.length();
        prev_header = str;
    }
}

void Log_output::write_new_header(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    if(newline || Settings::noise > 0)
    {
        Log_output::write_out(str+"\n",priority);
    }
    else
    {
        Log_output::clean_output();

        prev_msg = " ";
        msg_length = 1;
        prev_msg2 = "";
        msg2_length = 0;


        *os<<str<<"\n"<<prev_msg;


        header_length = str.length();
        prev_header = str;
    }
}

void Log_output::write_warning(const string str,const int priority)
{
    if(priority>Settings::noise)
        return;

    if(newline || Settings::noise > 0)
    {
        Log_output::write_out(str+"\n",priority);
    }
    else
    {
        Log_output::clean_output();
        *os<<str<<"\n"<<prev_header<<"\n"<<prev_msg<<prev_msg2;
    }

}

void Log_output::clean_output()
{
    for(int i=0;i<msg_length+msg2_length;i++)
        *os<<'\b';
   for(int i=0;i<msg_length+msg2_length;i++)
            *os<<' ';
   for(int i=0;i<msg_length+msg2_length;i++)
       *os<<'\b';

   if(msg_length+msg2_length>0)
       *os << "\e[A";

   for(int i=0;i<header_length;i++)
       *os<<"\b";
   for(int i=0;i<header_length;i++)
            *os<<' ';
   for(int i=0;i<header_length;i++)
       *os<<'\b';

}
