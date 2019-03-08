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

#ifndef TREE_NODE_H
#define TREE_NODE_H

/*
 * Clumsy way of re-rooting
 */

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

namespace ppa {

class Tree_node
{
    Tree_node(std::string t,Tree_node* p,int branch);

    std::string tree;
    std::string sub_trees[2];
    std::string rev_trees[2];

    Tree_node* parent;
    Tree_node* child0;
    Tree_node* child1;

    float sub_distances[2];
    float max_length;

    bool is_last;
    bool is_unrooted;

    void find_middle_point();
    void find_middle(int branch);
    void divide_tree(std::string tree,std::string* trees,float* distances);

    bool node_has_left_child;
    bool node_has_right_child;

    bool is_leaf() { return is_last; }
    void is_leaf(bool i) { is_last = i; }

    bool has_left_child() { return node_has_left_child; }
    void has_left_child(bool h) { node_has_left_child = h; }

    bool has_right_child() { return node_has_right_child; }
    void has_right_child(bool h) { node_has_right_child = h; }


    string name;
    void set_name(string n) { name = n; }
    string get_name() { return name; }

    double dist_to_parent;
    void set_distance_to_parent(double d) { dist_to_parent = d; }
    double get_distance_to_parent() { return dist_to_parent; }

    void add_left_child(Tree_node *child)
    {
        child0 = child;
        is_leaf(false);
        this->has_left_child(true);
    }

    void add_right_child(Tree_node *child)
    {
        child1 = child;
        is_leaf(false);
        this->has_right_child(true);
    }

    void delete_left_child() { delete child0; }
    void delete_right_child() { delete child1; }

    static int count;
public:
    Tree_node() {}
    ~Tree_node();

    string get_rooted_tree(std::string t);
};
}

#endif
