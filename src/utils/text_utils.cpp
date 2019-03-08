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
// File: TextTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Aug  8 12:57:50 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#include "utils/text_utils.h"

using namespace ppa;

#include <ctype.h>
#include <sstream>
#include <iomanip>

using namespace std;


/******************************************************************************/

bool Text_utils::is_empty(const string & s)
{
  for(unsigned int i = 0; i < s.size(); i++)
  {
    char c = s[i];
    if(c != ' ' && c != '\n' && c != '\t') return false;
  }
  return true;
}

/******************************************************************************/

string Text_utils::to_upper(const string & s)
{
  string result = "";
  for(unsigned int i = 0; i < s.size(); i++)
  {
    result += toupper(s[i]);
  }
  return result;
}

/******************************************************************************/

string Text_utils::to_lower(const string & s)
{
  string result = "";
  for(unsigned int i = 0; i < s.size(); i++)
  {
    result += tolower(s[i]);
  }
  return result;
}

/******************************************************************************/

bool Text_utils::is_whitespace_character(char c)
{
    return (c == ' ')
        || (c == '\t')
        || (c == '\n')
        || (c == '\r')
        || (c == '\f');
}

/******************************************************************************/

string Text_utils::remove_whitespaces(const string & s)
{
  // Copy sequence
  string st (s);

  // For all sequence's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if(is_whitespace_character(st[i]))
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

string Text_utils::remove_first_whitespaces(const string & s)
{
  // Copy sequence
  string st (s);

  while(st.size() > 0 && is_whitespace_character(st[0]))
  {
    st.erase(st.begin());
  }

  // Send result
  return st;
}

/******************************************************************************/

string Text_utils::remove_last_whitespaces(const string & s)
{
  // Copy sequence
  string st (s);

  while(st.size() > 0 && is_whitespace_character(st[st.size() - 1]))
  {
    st.erase(st.end() - 1);
  }

  // Send result
  return st;
}

/******************************************************************************/

string Text_utils::remove_surrounding_whitespaces(const string & s)
{
  return remove_first_whitespaces(remove_last_whitespaces(s));
}

/******************************************************************************/

bool Text_utils::is_newline_character(char c)
{
  return (c == '\n')
      || (c == '\r');
}

/******************************************************************************/

string  Text_utils::remove_newlines(const string & s)
{
  // Copy string
  string st (s);

  // For all string's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if (is_newline_character(st[i]))
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

string Text_utils::remove_last_newlines(const string & s)
{
  // Copy string
  string st (s);

  while (st.size() > 0 && is_newline_character(st[st.size() - 1]))
  {
    st.erase(st.end() - 1);
  }

  // Send result
  return st;
}

/******************************************************************************/

bool Text_utils::is_decimal_number(char c)
{
  if(c == '0' || c == '1' || c == '2' || c == '3' || c == '4'
  || c == '5' || c == '6' || c == '7' || c == '8' || c == '9') return true;
  else return false;
}

/******************************************************************************/

bool Text_utils::is_decimal_number(const string & s, char dec)
{
  unsigned int sep_count = 0;
  for(unsigned int i = 0; i < s.size(); i++)
  {
    char c = s[i];
    if(c == dec) sep_count++;
    else if(!is_decimal_number(c)) return false;
    if(sep_count > 1) return false;
  }
  return true;
}

/******************************************************************************/

string Text_utils::to_string(int i)
{
  ostringstream oss;
  oss << i;
  return oss.str();
}

/******************************************************************************/

string Text_utils::to_string(char c)
{
  ostringstream oss;
  oss << c;
  return oss.str();
}

/******************************************************************************/

string Text_utils::to_string(double d, int precision)
{
  ostringstream oss;
  oss << setprecision(precision) << d;
  return oss.str();
}

/******************************************************************************/

int Text_utils::to_int(const string & s)
{
  istringstream iss(s);
  int i;
  iss >> i;
  return i;
}

/******************************************************************************/

double Text_utils::to_double(const string & s)
{
  istringstream iss(s);
  double d;
  iss >> d;
  return d;
}

/******************************************************************************/

string Text_utils::resize_right(const string & s, unsigned int new_size, char fill)
{
  if(s.size() > new_size) return s.substr(0, new_size);
  else return s + string(new_size - s.size(), fill);
}

/******************************************************************************/

string Text_utils::resize_left(const string & s, unsigned int new_size, char fill)
{
  if(s.size() > new_size) return s.substr(s.size() - new_size);
  else return string(new_size - s.size(), fill) + s;
}

/******************************************************************************/

vector<string> Text_utils::split(const string & s, unsigned int n)
{
  vector<string> v;
  string tmp = s;
  while(tmp.size() > n)
  {
    v.push_back(tmp.substr(0, n));
    tmp = tmp.substr(n);
  }
  v.push_back(tmp);
  return v;
}

/******************************************************************************/

string Text_utils::remove_substrings(const string & s, char block_beginning, char block_ending)
throw (Exception)
{
  string t = "";
  int block_count = 0;
  int beg_pos = 0;
  for(unsigned int i = 0; i < s.size(); i++)
  {
    char current = s[i];
    if(current == block_beginning)
    {
      block_count++;
      t += s.substr(beg_pos, i);
    }
    else if(current == block_ending)
    {
      block_count--;
      if(block_count == 0) {
        beg_pos = i + 1;
      }
      else if(block_count < 0)
        throw Exception("Text_utils::remove_substrings(). " +
          string("Ending block character without corresponding beginning one at position ") + to_string((int)i) + ".");
    }
  }
  t += s.substr(beg_pos, s.npos);
  return t;
}

/******************************************************************************/

string Text_utils::remove_char(const string & s, char c)
{
  // Copy sequence
  string st(s);

  // For all sequence's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if (st[i] == c)
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

unsigned int Text_utils::count(const string & s, const string & pattern)
{
  unsigned int count = 0;
  string::size_type index = s.find(pattern);
  while(index != string::npos)
  {
    count++;
    index = s.find(pattern, index+1);
  }
  return count;
}

/******************************************************************************/

//
// File: String_tokenizer.cpp
// Author : Julien Dutheil
//          Sylvain Gaillard
// Last modification : Monday September 20 2004
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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


using namespace ppa;

String_tokenizer::String_tokenizer(const string & s, const string & delimiters, bool solid, bool allow_empty_tokens)
{
    tokens = deque<string>();
    if (!solid)
    {
        string::size_type index = s.find_first_not_of(delimiters, 0);
        while(index != s.npos)
        {
            string::size_type new_index = s.find_first_of(delimiters, index);
            if(new_index != s.npos)
            {
                tokens.push_back(string(s.begin() + index, s.begin() + new_index));
                if(!allow_empty_tokens)
                    index = s.find_first_not_of(delimiters, new_index);
                else
                    index = new_index + 1;
            }
            else
            {
                tokens.push_back(string(s.begin() + index, s.end()));
                index = new_index;
            }
        }
    }
    else
    {
        string::size_type index = 0;
        while(index != s.npos)
        {
            string::size_type new_index = s.find(delimiters, index);
            if(new_index != s.npos)
            {
                tokens.push_back(string(s.begin() + index, s.begin() + new_index));
                if(!allow_empty_tokens)
                {
                    index = new_index + delimiters.size();
                    while(index != string::npos && s.substr(index, delimiters.size()) == delimiters) index++;
                }
                else
                    index = new_index + delimiters.size();
            }
            else
            {
                tokens.push_back(string(s.begin() + index, s.end()));
                index = new_index;
            }
        }
    }
    current_position = 0;
}

string String_tokenizer::next_token() throw (Exception)
{
    if(!has_more_token()) throw Exception("No more token in tokenizer.");
    return tokens[current_position++];
}

void String_tokenizer::remove_empty_tokens()
{
    for(unsigned int i = tokens.size(); i > current_position; i--)
    {
        if(tokens[i-1] == "") tokens.erase(tokens.begin() + i - 1);
    }
}

