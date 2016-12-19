/******************************************************************************
   Copyright (C) 2015-2016 Einar J.M. Baumann <einar.baumann@gmail.com>

   This file and the WellIndexCalculator as a whole is part of the
   FieldOpt project. However, unlike the rest of FieldOpt, the
   WellIndexCalculator is provided under the GNU Lesser General Public
   License.

   WellIndexCalculator is free software: you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation, either version 3 of
   the License, or (at your option) any later version.

   WellIndexCalculator is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with WellIndexCalculator.  If not, see
   <http://www.gnu.org/licenses/>.
******************************************************************************/

/*!
 * @brief This file contains helper methods for the main.cpp file in this folder.
 *
 * It contains methods to parse the program arguments, as well as well as methods to
 * print the output in various formats.
 */

#ifndef WIC_MAIN_H
#define WIC_MAIN_H

#include "intersected_cell.h"
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <iostream>
#include <stdlib.h>
#include <boost/filesystem/operations.hpp>

namespace po = boost::program_options;
using namespace Reservoir::WellIndexCalculation;
using namespace std;

void printCsv(vector<IntersectedCell> &well_blocks) {
    cout << "i,\tj,\tk,\twi" << endl;
    for (auto block : well_blocks) {
        auto line = boost::str(boost::format("%d,\t%d,\t%d,\t%s")
                %(block.ijk_index().i() + 1)         // %1
                %(block.ijk_index().j() + 1)         // %2
                %(block.ijk_index().k() + 1)         // %3
                %block.well_index());                // %4
        cout << line << endl;
    }
}

void printCompdat(vector<IntersectedCell> &well_blocks, string well_name, double wellbore_radius) {
    string head = "COMPDAT\n";
    string foot = "\n/";
    vector<string> body;
    for (auto block : well_blocks) {
        //                                      NAME  I    J  K1  K2 OP/SH ST WI  RAD BS
        auto entry = boost::str(boost::format("   %s  %d  %d  %d  %d OPEN  1  %s  %s  %s")
                % well_name             // %1
                %(block.ijk_index().i() + 1) // %2
                %(block.ijk_index().j() + 1) // %3
                %(block.ijk_index().k() + 1) // %4
	            %(block.ijk_index().k() + 1) // %5
                %block.well_index()          // %5
                %wellbore_radius
                %"/");           // %6
        body.push_back(entry);
    }
    string full = head + boost::algorithm::join(body, "\n") + foot;
    cout << full << endl;
}

po::variables_map createVariablesMap(int argc, const char **argv) {
    //
    // This function parses the runtime arguments and creates a boost::program_options::variable_map from them.
    // It also displays help if the --help flag is passed.
    //
    po::options_description desc("FieldOpt options");
    desc.add_options()
            ("help", "print help message")
            ("grid,g", po::value<string>(),
             "path to model grid file (e.g. *.GRID)")
            ("heel,h", po::value<vector<double>>()->multitoken(),
             "Heel coordinates (x y z)")
            ("toe,t", po::value<vector<double>>()->multitoken(),
             "Toe coordinates (x y z)")
            ("radius,r", po::value<double>(),
             "wellbore radius")
            ("compdat,c", po::value<int>()->implicit_value(0),
             "print in compdat format instead of CSV")
            ("well-name,w", po::value<string>(),
             "well name to be used when writing compdat")
            ;
	
    // Process arguments to variable map
    po::variables_map vm;
    
    // Parse the input arguments and store the values 
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    
    // ??
    po::notify(vm);

    // If called with --help or -h flag:
    if (vm.count("help")) { // Print help if --help present or input file/output dir not present
        cout << "Usage: ./WellIndexCalculator --grid gridpath --heel x1 y1 z1 --toe x2 y2 z2 --radius r [options]" << endl;
        cout << desc << endl;
        exit(EXIT_SUCCESS);
    }

    assert(vm.count("grid"));
    assert(vm.count("heel"));
    assert(vm.count("toe"));
    assert(vm.count("radius"));
    if (vm.count("compdat"))
        assert(vm.count("well-name"));
    assert(vm["heel"].as<vector<double>>().size() == 3);
    assert(vm["toe"].as<vector<double>>().size() == 3);
    assert(boost::filesystem::exists(vm["grid"].as<string>()));
    assert(vm["radius"].as<double>() > 0);

    return vm;
}

#endif // WIC_MAIN_H
