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

#ifndef FIND_ANCHORS_H
#define FIND_ANCHORS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "substring_hit.h"
#include "main/sequence.h"

/*
Partly adapted from http://www.drdobbs.com/architecture-and-design/algorithm-alley/184404588
*/

namespace ppa{

struct Coord{
    int x;
    int y;
    Coord() : x(-1), y(-1){

    }
};
enum block_type{
    normal,
    begin,
    end,
    full_matrix
};
struct Tunnel_block{
    Coord start;
    Coord end;

    long long int size() const{
        return (long long int)(end.x - start.x) * (long long int)(end.y - start.y);
    }
    bool operator > (const Tunnel_block& tb) const{
        return (size() > tb.size());
    }
    bool operator < (const Tunnel_block& tb) const{
        return (size() < tb.size());
    }
    void print(){
        std::cout << "start (x,y): " << start.x << ", " << start.y << endl;
        std::cout << "upper end (x,y): " << end.x << ", " << end.y << endl;
        std::cout << "size: " << size() << endl;
        std::cout << "----------------------------" << endl;
    }
};

class Find_anchors
{
    static int pstrcmp(const void *p, const void *q)
    {
        const char **cp = (const char **)p;
        const char **cq = (const char **)q;
        return strcmp(*cp, *cq);
    }

    static bool sort_by_length(Substring_hit p, Substring_hit q)
    {
        return (p.length > q.length);
    }

    static bool sort_by_score(Substring_hit p, Substring_hit q)
    {
        return (p.score > q.score);
    }

    static bool sort_by_start_site_1(Substring_hit p, Substring_hit q)
    {
        if(p.start_site_1 == q.start_site_1)
            return (p.start_site_2 < q.start_site_2);

        return (p.start_site_1 < q.start_site_1);
    }

    int identical_prefix_length(char *p, char *q)
    {
        int i = 0;
        while (*p && (*p++ == *q++))
            i++;
        return i;
    }

    bool different_strings(char *p, char *q)
    {
        if(p-beg1>=0 && p-beg1<len1 && q-beg2>=0 && q-beg2<len2)
            return true;
        if(q-beg1>=0 && q-beg1<len1 && p-beg2>=0 && p-beg2<len2)
            return true;

        return false;
    }

    int position_string1(char *p, char *q)
    {
        if(p-beg1>=0 && p-beg1<len1)
            return p-beg1;
        if(q-beg1>=0 && q-beg1<len1)
            return q-beg1;

        return -1;
    }

    int position_string2(char *p, char *q)
    {
        if(p-beg2>=0 && p-beg2<len2)
            return p-beg2;
        if(q-beg2>=0 && q-beg2<len2)
            return q-beg2;

        return -1;
    }

    bool hits_overlap(Substring_hit h1,Substring_hit h2)
    {
        if(h1.start_site_1 <= h2.start_site_1 && h1.start_site_1+h1.length >= h2.start_site_1 )
            return true;
        if(h1.start_site_2 <= h2.start_site_2 && h1.start_site_2+h1.length >= h2.start_site_2 )
            return true;

        return false;
    }


    char *beg1;
    char *beg2;
    int len1;
    int len2;

    int overlapsAtBegin(Substring_hit& hit, Substring_hit& subject);
    bool probaplyBadHit(Substring_hit& hit, Substring_hit& subject);
    bool totallyOverlappingHit(Substring_hit& hit, Substring_hit& subject);
    bool partlyOverlappingHit(Substring_hit& hit, Substring_hit& subject);
    unsigned int distance(Substring_hit& hit, Substring_hit& subject);
    void plotR(std::vector<std::pair<int, int> >& hits, std::vector<int>& upper, std::vector<int>& lower, int l1, int l2);



public:
    Find_anchors();
    void find_long_substrings(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits,int min_length);
    void find_hmmer_anchors(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits);
    void check_hits_order_conflict(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits);
    void define_tunnel(std::vector<Substring_hit> *hits,std::vector<int> *upper_bound,std::vector<int> *lower_bound,std::string *str1, std::string *str2);
    void define_tunnel_with_overlapping_hits(std::vector<Substring_hit>& all_hits, std::vector<int>& upper, std::vector<int>& lower, std::string& sequence1, std::string& sequence2, int width, vector<Tunnel_block>& empty_blocks);
    void eliminate_bad_hits(std::vector<Substring_hit>& hits, unsigned int threshold_for_bads, unsigned int threshold_for_goods);
    void plotRFile(std::vector<std::pair<int, int> >& hits, std::vector<int>& upper, std::vector<int>& lower, int l1, int l2, string filename);
};

}
#endif // FIND_ANCHORS_H
