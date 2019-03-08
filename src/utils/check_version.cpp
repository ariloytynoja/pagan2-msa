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

#include "utils/check_version.h"
#include "utils/log_output.h"

#include <stdio.h>
#include <curl/curl.h>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace ppa;

#define PORT 80

Check_version::Check_version(float version)
{

    stringstream ss;
    ss<<"\nThis is PAGAN v."<<version<<".\nChecking if updates are available at https://github.com/ariloytynoja/pagan-msa.\n";
    Log_output::write_out(ss.str(),0);

    CURL *curl;
    CURLcode res;
    std::string readBuffer;

    curl = curl_easy_init();
    if(curl) {
      curl_easy_setopt(curl, CURLOPT_URL, "https://raw.githubusercontent.com/ariloytynoja/pagan-msa/master/VERSION_HISTORY");
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
      res = curl_easy_perform(curl);
      curl_easy_cleanup(curl);

      bool print_this = true;
      bool has_printed = false;
      string s;
      stringstream stream(readBuffer);
      while( getline(stream,s) )
      {
          istringstream ss(s);
          double d;
          char v,p;
          while( ss >> v >> p >> d )
          {
              if(v=='v' && p=='.' && int(d*10000) <= int(version*10000)+10)
              {
                 print_this = false;
              }
          }

          if(print_this)
          {
              if(!has_printed)
                  Log_output::write_out("\nFound updates. Changes in the more recent versions:\n\n",0);

              has_printed = true;
              Log_output::write_out(s+"\n",0);
          }
          else
          {
              break;
          }
      }

      if(!has_printed)
          Log_output::write_out("\nNo updates are available.\n\n",0);

    }
    return;

}
