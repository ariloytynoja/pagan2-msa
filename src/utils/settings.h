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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <boost/program_options.hpp>
#include <string>

namespace ppa{

class Settings
{
    boost::program_options::variables_map vm;
    boost::program_options::options_description full_desc;
    boost::program_options::options_description desc;
    boost::program_options::options_description min_desc;
    boost::program_options::options_description max_desc;

    float version;
    std::string date;


public:
    Settings();
    int read_command_line_arguments(int argc, char *argv[]);

    bool is(std::string name) { return vm.count(name); }
    const boost::program_options::variable_value & get(const std::string & name) const { return vm[name]; }

    void help();
    void help_all();
    void info();
    void info_noexit();
    void check_version();
    void print_msg();
    std::string print_log_msg();

    enum Placement_target_nodes {tid_nodes,terminal_nodes,internal_nodes,all_nodes};

    static int noise;
    static float resize_factor;

    static int exonerate_local_keep_best;
    static int exonerate_gapped_keep_best;

    static float tunneling_coverage;

    static int placement_target_nodes;
    static bool placement_preselection;

};

}
#endif // SETTINGS_H
