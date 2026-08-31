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

#ifndef TEMP_DIR_H
#define TEMP_DIR_H

#include "utils/settings_handle.h"
#include <string>
#include <cstdlib>
#include <sys/stat.h>

namespace ppa {

/*
 * Directory for the scratch files PAGAN writes for its helper programs
 * (Exonerate, MAFFT, RAxML, FastTree, BppAncestor, BppDist, BppPhySamp).
 *
 * Precedence, highest first:
 *
 *   1. --temp-folder, when given: an explicit choice always wins.
 *   2. $TMPDIR, $TMP, $TEMP: the first one that is set and non-empty.
 *   3. /tmp: the historical hardcoded default.
 *
 * Consulting the environment matters wherever /tmp is not a good place to
 * write: it is a small RAM-backed tmpfs on many desktop installations, and on
 * a cluster execution node it is routinely a tiny shared partition while the
 * space the job actually reserved is exported as $TMPDIR by the batch system.
 * Every other program in a pipeline honours these variables, so a user who
 * sets them is entitled to expect PAGAN to stay inside them too.
 *
 * Returns the directory with a trailing '/', or an empty string if it does
 * not exist -- preserving the long-standing fallback of writing into the
 * current working directory rather than failing.
 */
inline std::string resolve_temp_dir()
{
    std::string tmp_dir = "/tmp/";

    static const char *env_names[] = { "TMPDIR", "TMP", "TEMP" };
    for(unsigned int i=0; i<sizeof(env_names)/sizeof(env_names[0]); i++)
    {
        const char *value = getenv(env_names[i]);
        if(value != 0 && *value != '\0')
        {
            tmp_dir = std::string(value)+"/";
            break;
        }
    }

    if(Settings_handle::st.is("temp-folder"))
        tmp_dir = Settings_handle::st.get("temp-folder").as<std::string>()+"/";

    struct stat st;
    if(stat(tmp_dir.c_str(),&st) != 0)
        tmp_dir = "";

    return tmp_dir;
}

}

#endif // TEMP_DIR_H
