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
// File: TextTools.h
// Created by: Julien Dutheil
// Created on: Fri Aug  8 12:57:50 2003
//

/*
   Copyright or � or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide basal and
   utilitary classes. This file belongs to the Bio++ Project.

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


#ifndef TEXT_UTILS_H
#define TEXT_UTILS_H

#include "utils/exceptions.h"

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;

namespace ppa
{

/**
* @brief Some utilitary functions that work on strings.
*/
class Text_utils
{
    public:

    /**
    * @brief Tell if a string is empty.
    *
    * A string is considered to be 'empty' if it is only made of white spaces.
    *
    * @param s The string to check.
    * @return True if the string has only white characters.
    */
    static bool is_empty(const string & s);

    /**
    * @brief Make the string uppercase.
    *
    * @param s The string to analyse.
    * @return A copy of the string with all chars uppercase.
    */
    static string to_upper(const string & s);

    /**
    * @brief Make the string lowercase.
    *
    * @param s The string to analyse.
    * @return A copy of the string with all chars lowercase.
    */
    static string to_lower(const string & s);

    /**
    * @brief Tell if a character is a white space or not.
    *
    * @param c The character to check.
    * @return True if c is one of the following: ' ', '\\t', '\\n', '\\r' or '\\f'.
    */
    static bool is_whitespace_character(char c);

    /**
    * @brief Remove all white spaces characters in a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' without white spaces characters.
    */
    static string remove_whitespaces (const string & s);

    /**
    * @brief Remove all white spaces characters at the beginning of a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' beginning with the first non-white character.
    */
    static string remove_first_whitespaces (const string & s);

    /**
    * @brief Remove all white spaces characters at the end of a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' ending with the last non-white character.
    */
    static string remove_last_whitespaces (const string & s);

    /**
    * @brief Remove all white spaces characters at the beginning and the
    * end of a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' beginning with the first non-white character
    * and ending with the last one.
    */
    static string remove_surrounding_whitespaces(const string & s);

    /**
    * @brief Tell if a character is a new line character or not.
    *
    * @param c The character to check.
    * @return True if c is one of the following: '\\n' or '\\r'.
    */
    static bool is_newline_character(char c);

    /**
    * @brief Remove all new line characters in a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' without new line characters.
    */
    static string remove_newlines (const string & s);

    /**
    * @brief Remove all new line characters at the end of a string.
    *
    * @param s The string to parse.
    * @return A copy of 's' ending with the last non-new line character.
    */
    static string remove_last_newlines(const string & s);

    /**
    * @brief Tell is a given character describes a decimal number.
    *
    * @param c The character to check.
    * @return true if the given character is the reprensentation of a decimal number.
    */
    static bool is_decimal_number(char c);

    /**
    * @brief Tell is a given character string describes a decimal number.
    *
    * NB: for now, this parser will not recognize thosands delimiters, and not the scientific notation neither.
    * @param s The string to parse.
    * @param dec The decimal separator.
    * @return true if the given string is the representation of a decimal number.
    */
    static bool is_decimal_number(const string & s, char dec = '.');

    /**
    * @brief General template method to convert to a string.
    *
    * @param t The object to convert.
    * @return A string equal to t.
    */
    template<class T> static string to_string(T t)
    {
    ostringstream oss;
    oss << t;
    return oss.str();
    }

    /**
    * @brief Template string conversion.
    *
    * @param t The object to convert.
    * @param precision To use (for numbers).
    * @return A string equal to t.
    */
    template<class T>
    static string to_string(T t, int precision)
    {
      ostringstream oss;
      oss << setprecision(precision) << t;
      return oss.str();
    }

    /**
    * @brief General template method to convert from string.
    *
    * @param s The string to convert.
    * @return An object from string t.
    */
    template<class T> static T from_string(const string & s)
    {
    istringstream iss(s);
    T obj;
    iss >> obj;
    return obj;
    }

    /**
    * @brief Convert from int to string.
    *
    * @param i The integer to convert.
    * @return A string equal to i.
    */
    static string to_string(int i);

    /**
    * @brief Convert from char to string.
    *
    * @param c The character to convert.
    * @return A string equal to c.
    */
    static string to_string(char c);

    /**
    * @brief Convert from double to string.
    *
    * @param d The double to convert.
    * @param precision To use (for numbers).
    * @return A string equal to d.
    */
    static string to_string(double d, int precision = 6);

    /**
    * @brief Convert from string to int.
    *
    * @param s The string to parse.
    * @return The integer corresponding to s.
    */
    static int to_int(const string & s);

    /**
    * @brief Convert from string to double.
    *
    * @param s The string to parse.
    * @return The double corresponding to s.
    */
    static double to_double(const string & s);

    /**
    * @brief Template to string conversion.
    *
    * @param s The string to parse.
    * @return An object of class R corresponding to s.
    */
    template<class T>
    static T to(const string & s)
    {
      istringstream iss(s);
      T t;
      iss >> t;
      return t;
    }

    /**
    * @brief Send a string of size 'new_size', which is a copy of 's' truncated or
    * filled with character 'fill' at the end.
    *
    * @param s       The string to parse.
    * @param new_size The new string size.
    * @param fill    The character to use to fill the string id length < new_size.
    * @return A string of size newsize which is a copy from the left of s.
    */
    static string resize_right(const string & s, unsigned int new_size, char fill = ' ');

    /**
    * @brief Send a string of size 'new_size', which is a copy of 's' truncated or
    * filled with character 'fill' at the beginning.
    *
    * @param s       The string to parse.
    * @param new_size The new string size.
    * @param fill    The character to use to fill the string id length < new_size.
    * @return A string of size newsize which is a copy from the right of s.
    */
    static string resize_left(const string & s, unsigned int new_size, char fill = ' ');

    /**
    * @brief Split a string into parts of size 'n'.
    *
    * The last part may contain < n chars.
    *
    * @param s The string to parse.
    * @param n The number of tokens.
    * @return A vector of strings with all tokens.
    */
    static vector<string> split(const string & s, unsigned int n);

    /**
    * @brief Remove substrings from a string.
    *
    * All substrings beginning with block_beginning
    * and ending with block_ending will be removed.
    * Nesting blocks are allowed, the most extern block will be removed.
    *
    * @param s The string to parse.
    * @param block_beginning The character specifying the beginning of each block.
    * @param block_ending    The character specifying the end of each block.
    * @return The string with all blocks removed.
    * @throw Exception If some blocks are not well formed.
    */
    static string remove_substrings(const string & s, char block_beginning, char block_ending)
    ;

    /**
    * @brief Remove all occurences of a character in a string.
    *
    * @param s The string to parse.
    * @param c The character to remove.
    * @return The string with all specified chars removed.
    */
    static string remove_char(const string & s, char c);

    /**
    * @brief Count the occurences of a given pattern in a string.
    *
    * @param s The string to search.
    * @param pattern The pattern to use.
    * @return The number of occurences of 'pattern' in 's'.
    */
    static unsigned int count(const string & s, const string & pattern);

    string &replace_all(string& context, const string& from, const string& to)
    {
        size_t look = 0;
        size_t found;

        while( (found = context.find(from, look)) != string::npos) {
            context.replace(found, from.size(), to);
            look = found + to.size();
        }
        return context;
    }

    string &replace_all(string& context, char from, char to)
    {
        size_t look = 0;
        size_t found;

        while( (found = context.find(from, look)) != string::npos) {
            context.replace(found, 1, 1, to);
            look = found ++;
        }
        return context;
    }

};


/**
 * @brief Inner class for parsing strings in Newick format.
 */
class Node_tokenizer
{
    protected:
    vector<string> tokens;
    mutable unsigned int current_position;

    public:
    Node_tokenizer(const string & description) : tokens(), current_position(0)
    {
        unsigned int tok_count = 0;
        int par_count = 0;
        unsigned int i;
        for(i = 0; i < description.size(); i++)
        {
            if(description[i] == '(') par_count++; //Another open parenthesis
            if(description[i] == ')') par_count--; //Another close parenthesis
            if(par_count < 0) throw IOException("Invalid tree description: closing parenthesis with no opening one, in " + description);
            if(description[i] == ',' && par_count == 0)
            {
                //New token found:
                tokens.push_back(description.substr(tok_count, i - tok_count));
                tok_count = i + 1;
            }
        }
        //Add last token:
        tokens.push_back(description.substr(tok_count));

        current_position = 0;
    }

    public:
    string next() const
    {
        string s = tokens[current_position];
        current_position++;
        return s;
    }
    bool has_next() const
    {
        return current_position < tokens.size();
    }
};

} //end of namespace ppa.

#endif // TEXT_UTILS_H

//
// File: String_tokenizer.h
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

#ifndef STRING_TOKENIZER_H_
#define STRING_TOKENIZER_H_

#include <deque>
#include <string>

using namespace std;

#include "exceptions.h"

namespace ppa
{

/**
 * @brief A tokenizer for strings.
 *
 * Splits a string according to a given (set of) delimiter(s).
 */
class String_tokenizer
{
    protected:

        /** @brief Where the tokens are stored. */
        deque<string> tokens;

        /** @brief the current position in the token list. */
        unsigned int current_position;

    public:

        /**
         * @brief Build a new String_tokenizer from a string.
         *
         * @param s                The string to parse.
         * @param delimiters       Chars that must be considered as delimiters.
         * @param solid            If true, delimiters is considered as a single bloc delimiter.
     * @param allow_empty_tokens Tell if empty tokens are allowed or should be ignored.
         */
        String_tokenizer(const string & s, const string & delimiters = " \t\n\f\r", bool solid = false, bool allow_empty_tokens = false);

        virtual ~String_tokenizer() {}

    public:

        /**
         * @brief Get the next available token.
         * If no token is availbale, throw an Exception.
         *
         * @return The next token if there is one.
         */
        string next_token() ;

        /**
         * @brief Tell if some token are still available.
         * @return True if some token are still available.
         */
        bool has_more_token() const { return current_position < tokens.size(); }

        /**
         * @brief Tell how many tokens are available.
         *
         * @return the number of tokens available.
         */
        int number_of_remaining_tokens() const { return tokens.size() - current_position; }

        /**
         * @brief Get a particular token.
         *
         * Do not move the iterator.
         *
         * @param pos The index of the token.
         * @return the token at position 'pos'.
         */
        string get_token(unsigned int pos) const { return tokens[pos]; }

        /**
         * @brief Retrieve all tokens.
         *
         * @return A reference toward the vector of tokens.
         */
        const deque<string> & get_tokens() const { return tokens; }

    /**
     * @brief remove all empty token from the current position.
     */
    void remove_empty_tokens();

};

} //end of namespace bpp.

#endif	//_STRING_TOKENIZER_H_

