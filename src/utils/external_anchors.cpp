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

#include "utils/external_anchors.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

#include "utils/log_output.h"
#include "utils/settings_handle.h"

using namespace std;
using namespace ppa;

External_anchors::External_anchors()
    : n_anchors(0)
{
    // Loaded here rather than lazily on first use.  PAGAN2 aligns node pairs
    // under OpenMP, and this object is a function-local static shared by all
    // of them; C++11 guarantees the constructor runs exactly once even when
    // several threads reach it together, whereas a lazy "if(!loaded) load()"
    // inside get_anchors() would be a data race.  After construction the
    // object is read-only, so get_anchors() needs no locking.
    if(Settings_handle::st.is("anchor-file"))
        load(Settings_handle::st.get("anchor-file").as<string>());
}

void External_anchors::load(const string& filename)
{
    ifstream input(filename.c_str());
    if(!input)
    {
        Log_output::write_warning("Anchor file '"+filename+"' could not be "
                                  "opened. Anchoring normally instead.",0);
        return;
    }

    string line;
    int line_no = 0;
    int n_bad = 0;

    while(getline(input,line))
    {
        line_no++;

        string::size_type hash = line.find('#');
        if(hash != string::npos)
            line = line.substr(0,hash);

        istringstream fields(line);
        string id1, id2;
        long start1 = 0, start2 = 0, length = 0, score = 0;

        if(!(fields >> id1 >> id2 >> start1 >> start2 >> length))
        {
            // A blank or comment-only line is not an error.
            if(!id1.empty())
            {
                stringstream ss;
                ss << "Anchor file '" << filename << "' line " << line_no
                   << ": expected 'seq1 seq2 start1 start2 length [score]'. "
                   << "Skipping.";
                Log_output::write_warning(ss.str(),0);
                n_bad++;
            }
            continue;
        }

        if(!(fields >> score))
            score = length;

        // Coordinates are 1-based and lengths positive; anything else would
        // silently place the tunnel somewhere it does not belong.
        if(start1 < 1 || start2 < 1 || length < 1)
        {
            stringstream ss;
            ss << "Anchor file '" << filename << "' line " << line_no
               << ": start coordinates must be >= 1 and length > 0 (got "
               << start1 << ", " << start2 << ", " << length << "). Skipping.";
            Log_output::write_warning(ss.str(),0);
            n_bad++;
            continue;
        }

        if(id1 == id2)
        {
            stringstream ss;
            ss << "Anchor file '" << filename << "' line " << line_no
               << ": both sequences are '" << id1 << "'. Skipping.";
            Log_output::write_warning(ss.str(),0);
            n_bad++;
            continue;
        }

        Substring_hit hit;
        // Stored 0-based, as the internal providers produce them; swapped if
        // needed so that site_1 always belongs to the lexicographically first
        // name, which is how the key is built.
        if(id1 < id2)
        {
            hit.start_site_1 = (int)(start1 - 1);
            hit.start_site_2 = (int)(start2 - 1);
        }
        else
        {
            hit.start_site_1 = (int)(start2 - 1);
            hit.start_site_2 = (int)(start1 - 1);
        }
        hit.length = (int)length;
        hit.score  = (int)score;

        anchors[make_key(id1,id2)].push_back(hit);
        n_anchors++;
    }

    stringstream ss;
    ss << "Read " << n_anchors << " external anchor(s) for "
       << anchors.size() << " sequence pair(s) from '" << filename << "'";
    if(n_bad > 0)
        ss << " (" << n_bad << " malformed line(s) skipped)";
    ss << ".";
    // Priority 0: one line for the whole run, and if the user passed
    // --anchor-file, confirming the file was found and parsed is exactly what
    // they need. A silent no-op would look identical to it working.
    Log_output::write_out(ss.str(),0);
}

bool External_anchors::get_anchors(const string& id1, const string& id2,
                                   vector<Substring_hit>* hits,
                                   bool to_protein,
                                   int len1, int len2) const
{
    map<Pair_key, vector<Substring_hit> >::const_iterator it =
        anchors.find(make_key(id1,id2));

    if(it == anchors.end())
        return false;

    // Anchors are stored under the sorted key, so site_1 belongs to the
    // lexicographically first name.  The caller's (id1,id2) is the alignment's
    // own left/right order, which need not be the same -- if it is not, the
    // coordinates have to be swapped back or each one lands on the wrong
    // sequence.  Getting this wrong does not misalign quietly: the tunnel
    // builder indexes per-site vectors directly and dies with
    // std::out_of_range.
    const bool swap_back = !(id1 < id2);

    const vector<Substring_hit>& supplied = it->second;
    int n_used = 0;
    int n_frame = 0;
    int n_range = 0;

    for(unsigned int i=0; i<supplied.size(); i++)
    {
        Substring_hit hit = supplied.at(i);

        if(to_protein)
        {
            // Nucleotide coordinates only map onto the translated sequence if
            // the anchor starts on a codon boundary and spans whole codons.
            // Rounding would shift the anchor into the wrong frame and put the
            // tunnel around an alignment that cannot be right, so refuse.
            if(hit.start_site_1 % 3 != 0 || hit.start_site_2 % 3 != 0
               || hit.length % 3 != 0)
            {
                n_frame++;
                continue;
            }
            hit.start_site_1 /= 3;
            hit.start_site_2 /= 3;
            hit.length      /= 3;
        }

        if(hit.length < 1)
            continue;

        if(swap_back)
        {
            int tmp = hit.start_site_1;
            hit.start_site_1 = hit.start_site_2;
            hit.start_site_2 = tmp;
        }

        // The anchor must fit inside both sequences.  A typo in the file, or
        // coordinates taken against a different version of a sequence, would
        // otherwise abort the run inside the tunnel builder.
        if(hit.start_site_1 < 0 || hit.start_site_2 < 0
           || hit.start_site_1 + hit.length > len1
           || hit.start_site_2 + hit.length > len2)
        {
            n_range++;
            continue;
        }

        hits->push_back(hit);
        n_used++;
    }

    if(n_range > 0)
    {
        stringstream ss;
        ss << "Skipped " << n_range << " external anchor(s) for '" << id1
           << "' / '" << id2 << "': they run past the end of the sequences ("
           << len1 << ", " << len2 << " sites). Check the coordinates are "
           << "1-based and taken against these sequences.";
        Log_output::write_warning(ss.str(),0);
    }

    if(n_frame > 0)
    {
        stringstream ss;
        ss << "Skipped " << n_frame << " external anchor(s) for '" << id1
           << "' / '" << id2 << "': not codon-aligned, and this is a codon "
           << "alignment. Coordinates must be nucleotide positions with "
           << "(start-1) and length both divisible by three.";
        Log_output::write_warning(ss.str(),0);
    }

    if(n_used == 0)
        return false;

    stringstream ss;
    ss << "Using " << n_used << " external anchor(s) for '" << id1
       << "' / '" << id2 << "'.";
    // Priority 2: per-pair detail, one line per aligned pair, so it stays
    // behind higher verbosity -- the load line above already confirms the
    // file took effect.
    Log_output::write_out(ss.str(),2);

    return true;
}
