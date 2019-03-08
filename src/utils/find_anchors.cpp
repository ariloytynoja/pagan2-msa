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

#include "utils/find_anchors.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"
#include <iostream>
#include <algorithm>

using namespace ppa;
using namespace std;

Find_anchors::Find_anchors()
{
}


void Find_anchors::find_long_substrings(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits,int min_length)
{
    len1 = seq1->length();
    len2 = seq2->length();

//    cout<<*seq1<<endl<<*seq2<<endl<<endl;
//    cout<<">root\n"<<*seq1<<endl;

    char c1[len1];
    char c2[len2];
    char *a[len1+len2];

    int m=0; int n=0;
    for(; n<len1; n++,m++) {
        c1[n] = seq1->at(n);
        a[m] = &c1[n];
    }
    c1[n] = 0;

    int l = m;

    n=0;
    for(; n<len2; n++,m++) {
        c2[n] = seq2->at(n);
        a[m] = &c2[n];
    }
    c2[n] = 0;

    beg1 = a[0];
    beg2 = a[l];

    qsort(a, m, sizeof(char *), Find_anchors::pstrcmp);

    for (int i = 0; i < m-1; i++)
    {
        if(different_strings(a[i], a[i+1]))
        {
            int length = this->identical_prefix_length(a[i], a[i+1]);
            if (length >= min_length) {
                int pos1 = position_string1(a[i],a[i+1]);
                int pos2 = position_string2(a[i],a[i+1]);

                Substring_hit s;
                s.start_site_1 = pos1;
                s.start_site_2 = pos2;
                s.length = length;
                s.score = length;
                hits->push_back(s);
            }
        }
    }

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_length);

    vector<bool> hit_site1;
    hit_site1.reserve(len1);
    for(int i=0;i<len1;i++)
        hit_site1.push_back(false);

    vector<bool> hit_site2;
    hit_site2.reserve(len2);
    for(int i=0;i<len2;i++)
        hit_site2.push_back(false);

    vector<Substring_hit>::iterator it1 = hits->begin();
    for(;it1!=hits->end();)
    {
        bool overlap = false;
        for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
        {
            if(hit_site1.at(i) || hit_site2.at(j))
            {
                overlap = true;
                break;
            }
        }


        if(overlap)
        {
            hits->erase(it1);
        }
        else
        {
            for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
            {
                hit_site1.at(i)=true;
                hit_site2.at(j)=true;
            }
            it1++;
        }
    }
}

void Find_anchors::find_hmmer_anchors(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits)
{
    len1 = seq1->length();
    len2 = seq2->length();

    ofstream h_in("hmmer_in.fas");
    h_in<<">1\n"<<*seq1<<endl<<">2\n"<<*seq2<<endl;

    stringstream command;
    command <<"hmmsearch --max  pagan.hmm hmmer_in.fas 2>&1";

    Log_output::write_out("Hmmer_anchors: command: "+command.str()+"\n",2);

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with hmmer pipe.\nExiting.\n",0);
        exit(1);
    }

    char line[256];

    int a,g,h,j,k,m,n;
    string b,i,l,o;
    float c,d,e,f,p;
    Substring_hit hitS;
    Substring_hit hitE;
    hitS.length = 5;
    hitE.length = 5;

    while ( fgets( line, sizeof line, fpipe))
    {
        string str(line);

        cout<<line;

        if(str.find(">> 1") != string::npos)
        {
            float minE=10000;

            char *no_warn = fgets( line, sizeof line, fpipe);
            no_warn = fgets( line, sizeof line, fpipe);
            while ( fgets( line, sizeof line, fpipe))
            {
                str = line;
                if(str.length()<=1)
                    break;

                stringstream strstr(str);
                strstr >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m >> n >> o >> p;

                if(e<minE)
                {
                    hitS.start_site_1 = j-g+2;
                    hitE.start_site_1 = k-h+22;
                    minE = e;
                }
                cout<<line;
//                cout<<e<<" "<<g<<" "<<h<<" "<<j<<" "<<k<<endl;
            }
        }
        if(str.find(">> 2") != string::npos)
        {
            float minE=10000;

            char *no_warn =  fgets( line, sizeof line, fpipe);
            no_warn = fgets( line, sizeof line, fpipe);
            while ( fgets( line, sizeof line, fpipe))
            {
                str = line;
                if(str.length()<=1)
                    break;

                stringstream strstr(str);
                strstr >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m >> n >> o >> p;

                if(e<minE)
                {
                    hitS.start_site_2 = j-g+2;
                    hitE.start_site_2 = k-h+22;
                    minE = e;
                }

                cout<<line;
//                cout<<e<<" "<<g<<" "<<h<<" "<<j<<" "<<k<<endl;
            }
        }
    }
    hits->push_back(hitS);
    hits->push_back(hitE);
    cout<<"anchor: "<<hitS.start_site_1<<" "<<hitS.start_site_2<<endl;
    cout<<"anchor: "<<hitE.start_site_1<<" "<<hitE.start_site_2<<endl;
    pclose(fpipe);

}

void Find_anchors::check_hits_order_conflict(std::string *seq1,std::string *seq2,vector<Substring_hit> *hits)
{
    len1 = seq1->length();
    len2 = seq2->length();

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_score);

    vector<bool> hit_site1;
    hit_site1.reserve(len1);
    for(int i=0;i<len1;i++)
        hit_site1.push_back(false);

    vector<bool> hit_site2;
    hit_site2.reserve(len2);
    for(int i=0;i<len2;i++)
        hit_site2.push_back(false);

//    cout<<"\nHits "<<hits->size()<<endl;
    int trim = Settings_handle::st.get("exonerate-hit-trim").as<int>();

    vector<Substring_hit>::iterator it1 = hits->begin();
    for(;it1!=hits->end();)
    {
        it1->start_site_1+trim;
        if(it1->start_site_1>len1) it1->start_site_1 == len1-1;
        it1->start_site_2+trim;
        if(it1->start_site_2>len1) it1->start_site_2 == len2-1;
        it1->length -= trim*2;

        bool overlap = false;
        for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
        {
            if(hit_site1.at(i) || hit_site2.at(j))
            {
                overlap = true;
                break;
            }
        }


        if(overlap)
        {
            hits->erase(it1);
        }
        else
        {
            for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
            {
                hit_site1.at(i)=true;
                hit_site2.at(j)=true;
            }
            it1++;
        }
    }

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_start_site_1);

    it1 = hits->begin();
    vector<Substring_hit>::iterator it2 = it1;
    it2++;
    for(;it1!=hits->end() && it2!=hits->end();)
    {
        if(it1->start_site_2 > it2->start_site_2)
        {
            if(it1->score < it2->score)
            {
                hits->erase(it1);
                it2 = it1;
                it2++;
            }
            else
            {
                hits->erase(it2);
                it2 = it1;
                it2++;
            }
            continue;
        }
        it1++;it2++;
    }


//    cout<<"Hits "<<hits->size()<<endl;

//    if(hits->size()>0)
//    {
//        vector<Substring_hit>::iterator it1 = hits->begin();
//        for(;it1!=hits->end();it1++)
//        {
//            cout<<"REMAINS "<<it1->start_site_1<<".."<<it1->start_site_1+it1->length<<" "<<it1->start_site_2<<".."<<it1->start_site_2+it1->length<<" ("<<it1->score<<")\n";
//        }
//    }
}


void Find_anchors::define_tunnel(std::vector<Substring_hit> *hits,std::vector<int> *upper_bound,std::vector<int> *lower_bound,string *str1,string *str2)
{
    int length1 = str1->length();
    int length2 = str2->length();

//    cout<<"\n>S1\n"<<str1<<"\n>S2\n"<<str2<<"\n";

    // Has to consider skpped-over gaps as those are in the sequence objects
    //
    vector<int> index1;
    vector<int> index2;

    for(int i=0;i<length1;i++)
        if(str1->at(i)!='-')
            index1.push_back(i+1);

    for(int i=0;i<length2;i++)
        if(str2->at(i)!='-')
            index2.push_back(i+1);


    upper_bound->reserve(length1+1);
    lower_bound->reserve(length1+1);

    vector<int> diagonals;
    diagonals.reserve(length1+1);
    for(int i=0;i<length1+1;i++)
        diagonals.push_back(-1);

    for(vector<Substring_hit>::iterator it = hits->begin();it!=hits->end();it++)
    {
//        cout<<it->start_site_1<<" "<<it->start_site_2<<" "<<it->length<<" "<<length1<<" "<<length2<<"\n";
        int i=0;
        for(;i<it->length;i++)
        {
            diagonals.at(index1.at(it->start_site_1+i)) = index2.at(it->start_site_2+i);
        }

        if(it->start_site_1+i < (int)index1.size() && index1.at(it->start_site_1+i) < (int)diagonals.size())
            diagonals.at(index1.at(it->start_site_1+i)) = -2;

    }


    int width = Settings_handle::st.get("anchors-offset").as<int>();

    int y1 = 0;
    int y2 = 0;
    int y;

    int prev_y = 0;
    int m_count = 0;

    for(int i=0;i<=length1;i++)
    {
        if(i>=width && diagonals.at(i-width)>=0)
            y1 = diagonals.at(i-width)+0;

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)-width+0;
        }

        if(diagonals.at(i)>=0 && i>0 && diagonals.at(i-1)+1 == diagonals.at(i))
        {
            m_count++;
        }
        else if(diagonals.at(i)==-2)
        {
            m_count = 0;
        }

        y = min(y1,y2);
        y = max(y,0);

        if(diagonals.at(i)>=0 && i>0 && diagonals.at(i-1)+1 == diagonals.at(i) && m_count>=width)
        {
            prev_y = y;
        }

//        if(i>0)
//            cout<<i<<" "<<y1<<" "<<y2<<" "<<prev_y<<"; "<<diagonals.at(i)<<" "<<diagonals.at(i-1)<<" "<<m_count<<endl;

        y = min(y,prev_y);
        y = max(y,0);

        upper_bound->push_back(y);

    }

    y1 = length2;
    y2 = length2;

    prev_y = length2;
    m_count = 0;

    for(int i=length1;i>=0;i--)
    {
        if(i<=length1-width && diagonals.at(i+width)>=0)
            y1 = diagonals.at(i+width)+0;

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)+width+0;
        }

        if(diagonals.at(i)>=0 && i<length1 && diagonals.at(i+1)-1 == diagonals.at(i))
        {
            m_count++;
        }
        else if(diagonals.at(i)==-2)
        {
            m_count = 0;
        }

        y = max(y1,y2);
        y = min(y,length2);

        if(diagonals.at(i)>=0 && i<length1 && diagonals.at(i+1)-1 == diagonals.at(i)&& m_count>=width)
        {
            prev_y = y;
        }

        y = max(y,prev_y);
        y = min(y,length2);

        lower_bound->insert(lower_bound->begin(),y);
    }

    /*
    if(Settings::noise>0)
    {
        int sum = 0;
        for(int i=0;i<length1;i++)
            sum += lower_bound->at(i)-upper_bound->at(i);

        stringstream s;
        s.precision(4);
        s<<"Anchoring: Computing "<<((float)sum/(length1*length2))*100.0<<"% of DP matrix.\n";
        Log_output::write_out(s.str(),1);
    }
    //*/

    if(Settings_handle::st.is("plot-anchors-for-R")){
        cout<<"\n\nl1="<<str1->length()<<"; l2="<<str2->length()<<"; plot(1,1,type=\"n\",xlim=c(0,l1),ylim=c(0,l2))\n";
        cout<<"segments(0,0,l1,0,col=\"red\"); segments(0,l2,l1,l2,col=\"red\"); segments(0,0,0,l2,col=\"red\"); segments(l1,0,l1,l2,col=\"red\")\n";
        stringstream xh;
        stringstream yh;
        int lastx = 0,lasty = 0;
        for(int i=0;i<(int)diagonals.size();i++)
        {
            if(diagonals.at(i)>=0)
                {
                xh<<i<<",";
                yh<<diagonals.at(i)<<",";
                lastx=i;lasty=diagonals.at(i);
            }
        }
        xh<<lastx;yh<<lasty;
        cout<<"xhit=c("<<xh.str()<<");yhit=c("<<yh.str()<<");points(xhit,yhit,pch=\".\")\n";
        cout<<"upper=c("<<upper_bound->at(0);
        for(int i=1;i<(int)upper_bound->size();i++)
            cout<<","<<upper_bound->at(i);
        cout<<")\nlower=c("<<lower_bound->at(0);
        for(int i=1;i<(int)lower_bound->size();i++)
            cout<<","<<lower_bound->at(i);
        cout<<")\nlines(0:l1,upper,col=\"green\"); lines(0:l1,lower,col=\"blue\")\n\n";
    }

}

/**
 * @brief Find_anchors::eliminate_bad_hits eliminates hits based on given thresholds
 * @param hits
 * @param threshold_totally_overlapping - max distance of hits that are fully overlapping with any of already chosen hits
 * @param threshold_partly_overlapping - max distance of hits that are partly overlapping with any of already chosen hits (should be higher than totally overlapping, otherwise weird results may occur)
 */
void Find_anchors::eliminate_bad_hits(std::vector<Substring_hit>& hits, unsigned int threshold_totally_overlapping, unsigned int threshold_partly_overlapping){


    vector<Substring_hit*> good_hits;
    //vector<Substring_hit*> decent_hits;

    std::vector<Substring_hit>::iterator hit = hits.begin();

    bool bad_hit = false;
    bool decent_hit = false;

    for (; hit != hits.end();){
        //std::cout << "Hitti: " << hit->score << std::endl;

        bad_hit = false;
        decent_hit = false;

        for (std::vector<Substring_hit*>::iterator s = good_hits.begin(); s != good_hits.end(); ++s){
            if(probaplyBadHit(*hit, **s) || totallyOverlappingHit(*hit, **s)){
                //std::cout << "huono hitti" << std::endl;
                if(distance(*hit, **s) > threshold_totally_overlapping){
                    bad_hit = true;
                    break;
                }else{
                    decent_hit = true;
                }

            }else if(partlyOverlappingHit(*hit, **s)){
                if(distance(*hit, **s) > threshold_partly_overlapping){
                    bad_hit = true;
                    break;
                }
            }


        }

        if(bad_hit){
            hits.erase(hit);
        }else{
            if(!decent_hit){
                good_hits.push_back(&(*hit));
            }
            ++hit;
        }

    }

}
/**
 * @brief Find_anchors::overlapsAtBegin returns the length of the hit's head overlapping subject's tail
 * @param hit
 * @param subject
 * @return
 */
int Find_anchors::overlapsAtBegin(Substring_hit& hit, Substring_hit& subject){

    int overlap = 0;

    if(hit.start_site_1 >= subject.start_site_1 && hit.start_site_1 + hit.length > subject.start_site_1 + subject.length){
        overlap = std::max(overlap, subject.start_site_1 + subject.length - hit.start_site_1);
    }

    if(hit.start_site_2 >= subject.start_site_2 && hit.start_site_2 + hit.length > subject.start_site_2 + subject.length){
        overlap = std::max(overlap, subject.start_site_2 + subject.length - hit.start_site_2);
    }

    return std::max(0, overlap);
}

/**
 * @brief Find_anchors::distance calculates the distance (vertical/horizontal) between two hits
 * @param hit
 * @param subject
 * @return
 */
unsigned int Find_anchors::distance(Substring_hit& hit, Substring_hit& subject){
    return abs((subject.start_site_1 - subject.start_site_2) - (hit.start_site_1 - hit.start_site_2));
}

/**
 * @brief Find_anchors::probaplyBadHit returns true if the hit conflicts with the subject (probaply just worse anchor)
 * @param hit
 * @param subject
 * @return
 */
bool Find_anchors::probaplyBadHit(Substring_hit& hit, Substring_hit& subject){
    if(hit.start_site_1 < subject.start_site_1 && hit.start_site_2 > subject.start_site_2 && hit.start_site_1 + hit.length < subject.start_site_1 + subject.length){
        return true;
    }
    if(hit.start_site_1 > subject.start_site_1 && hit.start_site_2 < subject.start_site_2 && hit.start_site_2 + hit.length < subject.start_site_2 + subject.length){
        return true;
    }
    return false;
}
/**
 * @brief Find_anchors::totallyOverlappingHit returns true if the hit is completely inside the subject range on either axis
 * @param hit
 * @param subject
 * @return
 */
bool Find_anchors::totallyOverlappingHit(Substring_hit& hit, Substring_hit& subject){
    if(hit.start_site_1 >= subject.start_site_1 && hit.start_site_1 + hit.length <= subject.start_site_1 + subject.length){
        return true;
    }
    if(hit.start_site_2 >= subject.start_site_2 && hit.start_site_2 + hit.length <= subject.start_site_2 + subject.length){
        return true;
    }
    return false;
}
/**
 * @brief Find_anchors::partlyOverlappingHit
 * Checks if hits are party overlapping. Does not care if the hit is otherwise bad or not.
 *
 * @param hit
 * @param subject
 * @return
 */
bool Find_anchors::partlyOverlappingHit(Substring_hit& hit, Substring_hit& subject){
    if(overlapsAtBegin(hit, subject)){
        return true;
    }
    if(overlapsAtBegin(subject, hit)){
        return true;
    }
    return false;
}
/**
 * @brief Find_anchors::define_tunnel_with_overlapping_hits creates a search area around the given hits (anchors)
 * @param hits - anchors
 * @param upper - upper bounds of the tunnel (given that sequence 1 is on y-axis)
 * @param lower - lower bounds
 * @param l1 - length of sequence 1
 * @param l2 - length of sequence 2
 * @param width - width of the area to add on both sides
 */
void Find_anchors::define_tunnel_with_overlapping_hits(std::vector<Substring_hit>& all_hits, std::vector<int>& upper, std::vector<int>& lower, std::string& sequence1, std::string& sequence2, int width, vector<Tunnel_block>& empty_blocks){


    ///only accept plus strand hits
    std::vector<Substring_hit> hits;
    for(std::vector<Substring_hit>::iterator hit = all_hits.begin(); hit != all_hits.end(); ++hit){
        if(hit->plus_strand_1 && hit->plus_strand_2){
            hits.push_back(*hit);
        }
    }

    int l1 = sequence1.length();
    int l2 = sequence2.length();

    ///must take into account gaps in sequences
    std::vector<int> i1;
    std::vector<int> i2;

    for(int i = 0; i < l1; i++){
        if(sequence1.at(i) != '-'){
            i1.push_back(i+1);
        }
    }

    for(int i = 0; i < l2; i++){
        if(sequence2.at(i) != '-'){
            i2.push_back(i+1);
        }
    }



    ///highest and lowest positions (y-axis) for hits
    int lowest_points[l1+1];
    int highest_points[l1+1];


    int min_height = 0;
    int max_height = l2;

    for(int i = 0; i <= l1; i++){
        lowest_points[i] = max_height + 1;
        highest_points[i] = min_height - 1;
    }

    ///defining limits for the hits
    for(std::vector<Substring_hit>::iterator hit = hits.begin(); hit != hits.end(); ++hit){
        for(int a = 0; a < hit->length; a++){
            if(i2.at(hit->start_site_2 + a/* - width*/) < lowest_points[i1.at(hit->start_site_1 + a)]){
                lowest_points[i1.at(hit->start_site_1 + a)] = std::max(i2.at(hit->start_site_2 + a/* - width*/), min_height);
            }
            if(i2.at(hit->start_site_2 + a/* + width*/) > highest_points[i1.at(hit->start_site_1 + a)]){
                highest_points[i1.at(hit->start_site_1 + a)] = std::min(i2.at(hit->start_site_2 + a/* + width*/), max_height);
            }
        }
    }


    int previous_lowest = min_height;
    int previous_highest = max_height;


    ///must not go zigzag
    previous_highest = highest_points[0];
    for(int i = 0; i <= l1; i++){
        if(highest_points[i] > min_height){    //skip gaps
            if(highest_points[i] < previous_highest){
                highest_points[i] = previous_highest;
            }
            previous_highest = highest_points[i];
        }
    }
    previous_lowest = lowest_points[l1];
    for(int i = l1; i >= 0; i--){
        if(lowest_points[i] < max_height){     //skip gaps
            if(lowest_points[i] > previous_lowest){
                lowest_points[i] = previous_lowest;
            }
            previous_lowest = lowest_points[i];
        }
    }

    ///detect empty blocks in tunnel for (possibly) later use
    Tunnel_block current_block;
    current_block.start.x = 0;
    current_block.start.y = 0;
    for(int i = 1; i <= l1; i++){
        ///start of a new block
        if(highest_points[i-1] >= min_height && highest_points[i] < min_height){
            current_block.start.x = i;
            current_block.start.y = highest_points[i-1];
        }
        ///end of a block
        else if(highest_points[i] >= min_height && highest_points[i-1] < min_height){
            if(lowest_points[i] > current_block.start.y){
                current_block.end.x = i;
                current_block.end.y = lowest_points[i];
                if(current_block.size() > 10)
                    empty_blocks.push_back(current_block);
            }
        }
        ///special case: tunnel ends when in empty block
        else if(i == l1 && highest_points[i] < min_height){
            if(max_height > current_block.start.y){
                current_block.end.x = i;
                current_block.end.y = max_height;
                if(current_block.size() > 10)
                    empty_blocks.push_back(current_block);
            }
        }
    }

    std::sort(empty_blocks.begin(), empty_blocks.end());


    previous_lowest = min_height;
    previous_highest = max_height;

    ///lower bounds for gaps
    for(int i = 0; i <= l1; i++){
        if(lowest_points[i] >= max_height){
            lowest_points[i] = previous_lowest;
        }
        previous_lowest = lowest_points[i];
    }

    ///upper bounds for gaps
    for(int i = l1; i >= 0; i--){
        if(highest_points[i] <= min_height){
            highest_points[i] = previous_highest;
        }
        previous_highest = highest_points[i];
    }

    ///force tunnel corners to be in matrix corners
    lowest_points[0] = min_height;
    highest_points[l1] = max_height;


    ///widen tunnel to given width (x2)
    ///first make it thicker on y-axis
    for(int i = 0; i <= l1; i++){
        if(highest_points[i] >= min_height){
            highest_points[i] = std::min(max_height, highest_points[i] + width);
        }
    }
    for(int i = 0; i <= l1; i++){
        if(lowest_points[i] <= max_height){
            lowest_points[i] = std::max(min_height, lowest_points[i] - width);
        }
    }


    ///then ensure the thickness on x-axis
    vector<pair<int, bool> > overflow_highest;   //vector of areas(starting points) that need to be widened
    for(int i = 1; i <= l1; i++){
        if((i+1 > l1 || highest_points[i] == highest_points[i+1]) && highest_points[i-1] < highest_points[i] - 1){     //case gapped
            overflow_highest.push_back(pair<int, bool>(i, true));
        }else if(highest_points[i-1] < highest_points[i] - 1){  //case overlapping hits
            overflow_highest.push_back(pair<int, bool>(i, false));
        }
    }
    for(int a = 0; a < (int) overflow_highest.size(); a++){
        int i = overflow_highest.at(a).first;
        if(overflow_highest.at(a).second){
            for(int x = i-1; x >= i-width && x >= 0 && highest_points[x] >= min_height; x--){
                highest_points[x] = std::max(highest_points[x], highest_points[i]);
            }
        }else{
            for(int x = i-1; x >= i-width && x >= 0 && highest_points[x] >= min_height; x--){
                highest_points[x] = std::max(highest_points[x], highest_points[x+1] - 1);
            }
        }
    }

    vector<pair<int, bool> > overflow_lowest;    //vector of areas(starting points) that need to be widened
    for(int i = l1-1; i >= 0; i--){
        if((i-1 < 0 || lowest_points[i] == lowest_points[i-1]) && lowest_points[i+1] > lowest_points[i] + 1){
            overflow_lowest.push_back(pair<int, bool>(i, true));
        }else if(lowest_points[i+1] > lowest_points[i] + 1){
            overflow_lowest.push_back(pair<int, bool>(i, false));
        }
    }
    for(int a = 0; a < (int) overflow_lowest.size(); a++){
        int i = overflow_lowest.at(a).first;
        if(overflow_lowest.at(a).second){
            for(int x = i+1; x <=i+width && x <= l1 && lowest_points[x] <= max_height; x++){
                lowest_points[x] = std::min(lowest_points[x], lowest_points[i]);
            }
        }else{
            for(int x = i+1; x <=i+width && x <= l1 && lowest_points[x] <= max_height; x++){
                lowest_points[x] = std::min(lowest_points[x], lowest_points[x-1] + 1);
            }
        }
    }




    upper.reserve(l1 + 1);
    lower.reserve(l1 + 1);

    ///extra points at the beginning (apparently pagan assumes sequence's starting index is 1 instead of 0 or smth)
    //upper.push_back(lowest_points[0]);
    //lower.push_back(highest_points[0]);

    ///Arrays to vectors (sequence 1 at y-axis, not x-axis)
    for(int i = 0; i <= l1; i++){
        upper.push_back(lowest_points[i]);
        lower.push_back(highest_points[i]);
    }


    std::vector<std::pair<int, int> > hit_coords;

    ///defining hit coordinates for plotting
    for(std::vector<Substring_hit>::iterator hit = hits.begin(); hit != hits.end(); ++hit){
        for(int a = 0; a < hit->length; a++){
            hit_coords.push_back(pair<int, int>(i1.at(hit->start_site_1 + a), i2.at(hit->start_site_2 + a)));
        }
    }

    if(Settings_handle::st.is("plot-anchors-for-R"))
    {
        plotR(hit_coords, upper, lower, l1, l2);
    }
    //plotRFile(hit_coords, upper, lower, l1, l2, "find_anchors.R");

}
/**
 * @brief Find_anchors::plotR plot tunnel graph to R-file
 * @param hits
 * @param upper
 * @param lower
 * @param l1
 * @param l2
 */
void Find_anchors::plotR(std::vector<std::pair<int, int> >& hits, std::vector<int>& upper, std::vector<int>& lower, int l1, int l2){

    std::cout << "\n\nl1=" << l1 << "; l2=" << l2 << "; plot(1,1,type=\"n\",xlim=c(0,l1),ylim=c(0,l2),xlab=\"seq1\",ylab=\"seq2\")\n";
    std::cout << "segments(0,0,l1,0,col=\"red\"); segments(0,l2,l1,l2,col=\"red\"); segments(0,0,0,l2,col=\"red\"); segments(l1,0,l1,l2,col=\"red\")\n";
    std::stringstream xh;
    std::stringstream yh;

    for(std::vector<std::pair<int, int> >::iterator hit_coord = hits.begin(); hit_coord != hits.end(); ++hit_coord){
        xh << hit_coord->first;
        yh << hit_coord->second;

        if(!(hit_coord + 1 == hits.end())){
            xh << ",";
            yh << ",";
        }
    }

    std::cout << "xhit=c(" << xh.str() << ")\nyhit=c(" << yh.str() << ")\npoints(xhit,yhit,pch=\".\")\n";
    std::cout << "upper=c(" << upper.at(0);
    for(int i=1;i<(int)upper.size();i++)
        std::cout << "," << upper.at(i);
    std::cout << ")\nlower=c(" << lower.at(0);
    for(int i=1;i<(int)lower.size();i++)
        std::cout << "," << lower.at(i);
    std::cout << ")\nlines(0:l1,upper,col=\"green\"); lines(0:l1,lower,col=\"blue\")\n\n";
}
void Find_anchors::plotRFile(std::vector<std::pair<int, int> >& hits, std::vector<int>& upper, std::vector<int>& lower, int l1, int l2, string filename){

    std::stringstream f;
    f << "\n\nl1=" << l1 << "; l2=" << l2 << "; plot(1,1,type=\"n\",xlim=c(0,l1),ylim=c(0,l2),xlab=\"seq1\",ylab=\"seq2\")\n";
    f << "segments(0,0,l1,0,col=\"red\"); segments(0,l2,l1,l2,col=\"red\"); segments(0,0,0,l2,col=\"red\"); segments(l1,0,l1,l2,col=\"red\")\n";
    std::stringstream xh;
    std::stringstream yh;

    for(std::vector<std::pair<int, int> >::iterator hit_coord = hits.begin(); hit_coord != hits.end(); ++hit_coord){
        xh << hit_coord->first;
        yh << hit_coord->second;

        if(!(hit_coord + 1 == hits.end())){
            xh << ",";
            yh << ",";
        }
    }

    f << "xhit=c(" << xh.str() << ")\nyhit=c(" << yh.str() << ")\npoints(xhit,yhit,pch=\".\")\n";
    f << "upper=c(" << upper.at(0);
    for(int i=1;i<(int)upper.size();i++)
        f << "," << upper.at(i);
    f << ")\nlower=c(" << lower.at(0);
    for(int i=1;i<(int)lower.size();i++)
        f << "," << lower.at(i);
    f << ")\nlines(0:l1,upper,col=\"green\"); lines(0:l1,lower,col=\"blue\")\n\n";

    std::ofstream myfile;
    myfile.open (filename.c_str());
    myfile << f.str();
    myfile.close();

}
/*
void Find_anchors::plotRFile2(std::vector<std::pair<int, int> >& hits, std::vector<int>& upper, std::vector<int>& lower, int l1, int l2, Sequence* s1, Sequence* s2, string filename){

    vector<Site>& sites_1 = *(s1->get_sites());
    vector<Site>& sites_2 = *(s2->get_sites());

    std::stringstream f;
    f << "\n\nl1=" << l1 << "; l2=" << l2 << "; plot(1,1,type=\"n\",xlim=c(0,l1),ylim=c(0,l2),xlab=\"seq1\",ylab=\"seq2\")\n";
    f << "segments(0,0,l1,0,col=\"red\"); segments(0,l2,l1,l2,col=\"red\"); segments(0,0,0,l2,col=\"red\"); segments(l1,0,l1,l2,col=\"red\")\n";

    std::stringstream sh;
    for(int i = 0; i < (int)sites_1.size(); i++){

        if(sites_1[i].site_type == Site::non_real || sites_1[i].site_type == 5){
            sh << "segments(";
            sh << i << ",0," << i << "," << l2;
            sh << ",col=\"grey\")\n";
        }

    }

    for(int i = 0; i < (int)sites_2.size(); i++){


        if(sites_2[i].path_state == 5 || sites_2[i].path_state == 6){
            sh << "segments(";
            sh << "0," << i << "," << l1 << "," << i;
            sh << ",col=\"yellow\")\n";
        }
        if((sites_2[i].site_type == Site::non_real || sites_2[i].site_type == 5) && !(sites_2[i].path_state == 5 || sites_2[i].path_state == 6)){
            sh << "segments(";
            sh << "0," << i << "," << l1 << "," << i;
            sh << ",col=\"red\")\n";
        }else if(sites_2[i].site_type == Site::non_real || sites_2[i].site_type == 5){
            sh << "segments(";
            sh << "0," << i << "," << l1 << "," << i;
            sh << ",col=\"purple\")\n";
        }


    }

    f << sh.str();

    std::stringstream xh;
    std::stringstream yh;

    for(std::vector<std::pair<int, int> >::iterator hit_coord = hits.begin(); hit_coord != hits.end(); ++hit_coord){
        xh << hit_coord->first;
        yh << hit_coord->second;

        if(!(hit_coord + 1 == hits.end())){
            xh << ",";
            yh << ",";
        }
    }

    f << "xhit=c(" << xh.str() << ")\nyhit=c(" << yh.str() << ")\npoints(xhit,yhit,pch=\".\")\n";
    f << "upper=c(" << upper.at(0);
    for(int i=1;i<(int)upper.size();i++)
        f << "," << upper.at(i);
    f << ")\nlower=c(" << lower.at(0);
    for(int i=1;i<(int)lower.size();i++)
        f << "," << lower.at(i);
    f << ")\nlines(0:l1,upper,col=\"green\"); lines(0:l1,lower,col=\"blue\")\n\n";

    std::ofstream myfile;
    myfile.open (filename.c_str());
    myfile << f.str();
    myfile.close();

}*/
