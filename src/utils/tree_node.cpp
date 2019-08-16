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


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "utils/tree_node.h"
#include "utils/text_utils.h"

using namespace std;
using namespace ppa;

std::string mpTree;
float halfLength;
float maxSpan;

Tree_node::~Tree_node()
{

    if (node_has_left_child)
    {
        delete child0;
        node_has_left_child = false;
    }
    if (node_has_right_child)
    {
        delete child1;
        node_has_right_child = false;
    }
}

int Tree_node::count = 1;

string Tree_node::get_rooted_tree(std::string t)
{
    t = Text_utils::remove_last_whitespaces(t);

    mpTree = "";
    max_length = 0.0;
    maxSpan = 0.0;
    is_last = false;

    this->count = 1;
    stringstream ss;
    ss<<Tree_node::count;
    this->set_name(ss.str());
    Tree_node::count++;

    node_has_left_child = true;
    node_has_right_child = true;

    tree = t;

    divide_tree(tree,sub_trees,sub_distances);

    sub_distances[0] = abs(sub_distances[0]);
    sub_distances[1] = abs(sub_distances[1]);

    float tot = sub_distances[0]+sub_distances[1];

    char num[10];
    sprintf(num,"%.5f",tot);

    rev_trees[0] = sub_trees[1]+":"+num;
    rev_trees[1] = sub_trees[0]+":"+num;

    child0 = new Tree_node(sub_trees[0],this,0);
    child1 = new Tree_node(sub_trees[1],this,1);


    float currPair = sub_distances[0]+child0->max_length+sub_distances[1]+child1->max_length;
    if (currPair > maxSpan)
    {
        maxSpan = currPair;
    }

    find_middle_point();

    return mpTree;
}

Tree_node::Tree_node(string t,Tree_node* p,int branch)
{
    tree = t;
    parent = p;
    max_length = 0;
    is_last = true;

    node_has_left_child = false;
    node_has_right_child = false;
    this->set_name(tree);

    if (tree.find(",",0)>0 && tree.find(",",0)<tree.length())
    {
        is_last = false;

        node_has_left_child = true;
        node_has_right_child = true;

        stringstream ss;
        ss<<Tree_node::count;
        this->set_name(ss.str());
        Tree_node::count++;

        divide_tree(tree,sub_trees,sub_distances);

        sub_distances[0] = abs(sub_distances[0]);
        sub_distances[1] = abs(sub_distances[1]);

        char num0[10];
        sprintf(num0,"%.5f",sub_distances[0]);
        char num1[10];
        sprintf(num1,"%.5f",sub_distances[1]);

        rev_trees[0] = "("+parent->rev_trees[branch]+","+sub_trees[1]+":"+num1+"):"+num0;
        rev_trees[1] = "("+parent->rev_trees[branch]+","+sub_trees[0]+":"+num0+"):"+num1;

        child0 = new Tree_node(sub_trees[0],this,0);
        child1 = new Tree_node(sub_trees[1],this,1);


        float currPair = sub_distances[0]+child0->max_length+sub_distances[1]+child1->max_length;
        if (currPair > maxSpan)
        {
            maxSpan = currPair;
        }

        if (sub_distances[0]+child0->max_length > sub_distances[1]+child1->max_length)
        {
            max_length = sub_distances[0]+child0->max_length;
        }
        else
        {
            max_length = sub_distances[1]+child1->max_length;
        }
    }
}

void Tree_node::find_middle_point()
{
    halfLength = maxSpan/2;

    if (halfLength >= child0->max_length && halfLength <= child0->max_length+sub_distances[0]+sub_distances[1])
    {

        float b0 = halfLength-child0->max_length;
        float b1 = sub_distances[0]+sub_distances[1]-b0;

        char num0[10];
        sprintf(num0,"%.5f",b0);
        char num1[10];
        sprintf(num1,"%.5f",b1);

        mpTree = "("+child0->tree+":"+num0+","+child1->tree+":"+num1+");";

        return;
    }

    child0->find_middle(1);
    child1->find_middle(0);

}

void Tree_node::find_middle(int branch)
{
    if (!is_last)
    {
        if (branch==0)
        {

            if (halfLength >= child0->max_length && halfLength <= child0->max_length+sub_distances[0])
            {

                float b0 = halfLength-child0->max_length;
                float b1 = sub_distances[0]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",sub_distances[1]);

                mpTree = "("+child0->tree+":"+num0+",("+parent->rev_trees[1]+","+sub_trees[1]+":"+num+"):"+num1+");";

                return;
            }

            if (halfLength >= child1->max_length && halfLength <= child1->max_length+sub_distances[1])
            {

                float b0 = halfLength-child1->max_length;
                float b1 = sub_distances[1]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",sub_distances[0]);

                mpTree = "("+child1->tree+":"+num0+",("+parent->rev_trees[1]+","+sub_trees[0]+":"+num+"):"+num1+");";

                return;
            }
            child0->find_middle(1);
            child1->find_middle(0);

        }
        else
        {

            if (halfLength >= child0->max_length && halfLength <= child0->max_length+sub_distances[0])
            {

                float b0 = halfLength-child0->max_length;
                float b1 = sub_distances[0]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",sub_distances[1]);

                mpTree = "("+child0->tree+":"+num0+",("+parent->rev_trees[0]+","+sub_trees[1]+":"+num+"):"+num1+");";

                return;
            }


            if (halfLength >= child1->max_length && halfLength <= child1->max_length+sub_distances[1])
            {

                float b0 = halfLength-child1->max_length;
                float b1 = sub_distances[1]-b0;

                char num0[10];
                sprintf(num0,"%.5f",b0);
                char num1[10];
                sprintf(num1,"%.5f",b1);

                char num[10];
                sprintf(num,"%.5f",sub_distances[0]);

                mpTree = "("+child1->tree+":"+num0+",("+parent->rev_trees[0]+","+sub_trees[0]+":"+num+"):"+num1+");";

                return;
            }
            child0->find_middle(1);
            child1->find_middle(0);
        }
    }
}

void Tree_node::divide_tree(string tree,string* trees,float* distance)
{
    trees[0] = "";

    if ((tree.substr(tree.length()-1)).compare(";")==0)
    {
        tree = tree.substr(0,tree.find_last_of(")")+1); // remove last ';' and anything after the last bracket
    }
    tree = tree.substr(1,tree.length()-2);     // remove first & last '('

    if (tree.at(0)!='(')   // only one taxon before midpoint comma
    {

        string tmp = tree.substr(0,tree.find(",",0));
        trees[0] = tmp;
        distance[0] = 0;
        if(tmp.find(":")!=string::npos)
        {
            trees[0] = tmp.substr(0,tmp.find(":",0));
            distance[0] = atof((tmp.substr(tmp.find(":",0)+1).c_str()));
        }
        tree = tree.substr(tree.find(",",0)+1);

        bool trifurc = false;
        int open = 0;
        for (unsigned int j = 0; j<tree.length(); j++)
        {
            if (tree.at(j)=='(')
            {
                open++;
            }
            else if (tree.at(j)==')')
            {
                open--;
            }
            if (j>0 && open==0 && tree.substr(j).find(",",0)<=tree.length())
            {
                trifurc = true;
            }
        }

        // correction for trifurcating root
        if (trifurc)
        {
            is_unrooted = true;
            trees[1] = "("+tree+")";
            distance[0] = distance[0]/2;
            distance[1] = distance[0];
        }
        else
        {
            trees[1] = tree;
            distance[1] = 0;
            if(tree.find(":")!=string::npos)
            {
                trees[1] = tree.substr(0,tree.find_last_of(":"));
                tmp = tree.substr(tree.find_last_of(":")+1);
                distance[1] = atof(tmp.c_str());
            }
        }

    }
    else
    {

        int open = 0;

        for (unsigned int i=0; i<tree.length(); i++)
        {

            // count parentheses that are "open"
            if (tree.at(i)=='(')
            {
                open++;
            }
            else if (tree.at(i)==')')
            {
                open--;
            }
            trees[0].append(tree.substr(i,1));

            if (open<=0)
            {
                distance[0] = 0;
                if(tree.at(i+1) == ':') {
                    string tmp = tree.substr(i+2,tree.find(",",i+2));
                    distance[0] = atof(tmp.c_str());
                }
                tree = tree.substr(tree.find(",",i)+1);

                bool trifurc = false;
                open = 0;
                for (unsigned int j = 0; j<tree.length(); j++)
                {
                    if (tree.at(j)=='(')
                    {
                        open++;
                    }
                    else if (tree.at(j)==')')
                    {
                        open--;
                    }

                    if (open==0 && tree.find(",",j)<=tree.length())
                    {
                        trifurc = true;
                    }
                }

                // correction for trifurcating root
                if (trifurc)
                {

                    is_unrooted = true;
                    trees[1] = "("+tree+")";
                    distance[0] = distance[0]/2;
                    distance[1] = distance[0];

                }
                else
                {

                    trees[1] = tree;
                    distance[1] = 0;
                    if(tree.find(":")!=string::npos)
                    {
                        trees[1] = tree.substr(0,tree.find_last_of(":"));
                        string tmp = tree.substr(tree.find_last_of(":")+1);
                        distance[1] = atof(tmp.c_str());
                    }
                }
                break;
            }
        }
    }
}

