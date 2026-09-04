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

#ifndef EXTERNAL_ANCHORS_H
#define EXTERNAL_ANCHORS_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "utils/substring_hit.h"

namespace ppa
{

using namespace std;

/**
 * Anchors supplied by the user instead of computed by Exonerate/BLAST.
 *
 * PAGAN2 already anchors well when it can find its own hits.  This exists for
 * the case where it cannot but the caller already knows the answer: the caller
 * has just run a search (BLAST, minimap2, an existing curated alignment) and
 * holds high-identity local hits that would take PAGAN2 a second search to
 * rediscover -- or that its own anchoring misses because the pair is too
 * divergent overall for a hit to survive the internal thresholds.
 *
 * The anchors are read once, cached, and served per aligned pair.  A pair with
 * no supplied anchors falls back to the normal provider, so a partial file is
 * useful: supply anchors for the hard pairs only and let PAGAN2 handle the
 * rest.
 *
 * File format -- whitespace/TAB separated, '#' comments and blank lines
 * ignored:
 *
 *     seq1_id  seq2_id  start1  start2  length  [score]
 *
 *  - *seq1_id*, *seq2_id*  names as they appear in the input alignment/tree.
 *    The pair may be given in either order; coordinates are swapped to match.
 *  - *start1*, *start2*    1-based start of the hit in each sequence.
 *  - *length*              hit length, in the same units as the starts.
 *  - *score*               optional; defaults to *length*.  Used only for
 *                          ranking hits against each other.
 *
 * Coordinates are always in the units of the INPUT sequences -- nucleotides
 * for nucleotide input, including when ``--codons`` is in use.  That is what
 * the caller has from a nucleotide search, so no conversion is asked of them.
 * When PAGAN2 anchors in protein space (``--codons``, or codon data) the
 * coordinates are converted internally; an anchor that is not codon-aligned
 * cannot be converted without inventing a frame, so it is reported and skipped
 * rather than rounded.
 */
class External_anchors
{
    typedef pair<string,string> Pair_key;

    //: Anchors per unordered sequence pair, stored under the sorted key.
    map<Pair_key, vector<Substring_hit> > anchors;
    int  n_anchors;

    static Pair_key make_key(const string& a, const string& b)
    {
        return (a < b) ? Pair_key(a,b) : Pair_key(b,a);
    }

    void load(const string& filename);

public:
    External_anchors();

    /**
     * Anchors for the pair (*id1*, *id2*), appended to *hits*.
     *
     * const, and the object is immutable after construction, so this is safe
     * to call concurrently from PAGAN2's OpenMP regions.
     *
     * Coordinates are returned oriented to the CALLER's argument order:
     * ``start_site_1`` belongs to *id1* and ``start_site_2`` to *id2*,
     * whichever way round the file listed them.
     *
     * *to_protein* converts nucleotide coordinates to protein ones, for the
     * codon models.  *len1* and *len2* are the sequence lengths in the space
     * the anchors will be used in; an anchor running past either end is
     * dropped with a warning rather than handed on, since the tunnel builder
     * indexes straight into per-site vectors and would abort.
     *
     * Returns true if at least one anchor was supplied for this pair, i.e. if
     * the caller should skip its own anchor search.
     */
    bool get_anchors(const string& id1, const string& id2,
                     vector<Substring_hit>* hits, bool to_protein,
                     int len1, int len2) const;

    //: Total anchors accepted across the file.
    int count() const { return n_anchors; }
};

}
#endif // EXTERNAL_ANCHORS_H
