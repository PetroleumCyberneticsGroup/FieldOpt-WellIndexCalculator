/******************************************************************************
   Copyright (C) 2015-2016 Einar J.M. Baumann <einar.baumann@gmail.com>
   Modified by M.Bellout (2017) <mathias.bellout@ntnu.no>

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
 * @brief This file contains helper methods for the main.cpp file in
 * this folder.
 * It contains methods to parse the program arguments, as well as well
 * as methods to print the output in various formats.
 */

#ifndef WIC_MAIN_H
#define WIC_MAIN_H

#include "intersected_cell.h"
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <iostream>
#include <stdlib.h>

// OV
#if _WIN32
#include <sys/stat.h>
    inline bool exists(const std::string& name)
    {
        struct stat buffer;
        return (stat (name.c_str(), &buffer) == 0);
    }
#else
#include <boost/filesystem/operations.hpp>
#endif


namespace po = boost::program_options;
using namespace Reservoir::WellIndexCalculation;
using namespace std;

void printCsv(vector<IntersectedCell> &well_blocks) {
  cout << "i,\tj,\tk,\twi" << endl;
  for (auto block : well_blocks) {
    auto line = boost::str(boost::format("%d,\t%d,\t%d,\t%s")
                               % (block.ijk_index().i() + 1)        // %1
                               % (block.ijk_index().j() + 1)        // %2
                               % (block.ijk_index().k() + 1)        // %3
                               % block.well_index());               // %4
    cout << line << endl;
  }
}

void printCompdat(vector<IntersectedCell>& well_blocks,
                  string well_name,
                  double wellbore_radius) {
  string head = "COMPDAT\n";
  string foot = "\n/";
  vector<string> body;
  //                      NAME   I   J  K1  K2  OP/SH ST  WI   RADIUS
  string compdat_frmt = "   %s  %d  %d  %d  %d  OPEN  1*  %8s  %8s /";
  for (auto block : well_blocks) {
    auto entry = boost::str(boost::format(compdat_frmt)
                                % well_name                   // %1
                                % (block.ijk_index().i() + 1) // %2
                                % (block.ijk_index().j() + 1) // %3
                                % (block.ijk_index().k() + 1) // %4
                                % (block.ijk_index().k() + 1) // %5
                                % block.well_index()          // %6
                                % wellbore_radius);           // %7
    body.push_back(entry);
  }
  string full = head + boost::algorithm::join(body, "\n") + foot;
  cout << full << endl;
}

po::variables_map createVariablesMap(int argc, const char **argv) {
  //
  // This function parses the runtime arguments and creates
  // a boost::program_options::variable_map from them.
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
       "well name to be used when writing compdat");

  // Positional arguments <- OV
  po::positional_options_description p;
  p.add("grid", 1);

  // Process arguments to variable map
  po::variables_map vm;

  // Parse the input arguments and store the values
  po::store(po::parse_command_line(
      argc, argv, desc, po::command_line_style::unix_style ^
      po::command_line_style::allow_short), vm);

  // ??
  po::notify(vm);

  // Print help if --help/-h present or input file/output dir not present
  string usage_msg = "Usage: ./wicalc "
      "--grid gridpath "
      "--heel x1 y1 z1 "
      "--toe x2 y2 z2 "
      "--radius r [options]";

  if (vm.count("help")) {
    cout << usage_msg  << endl;
    cout << desc       << endl;
    exit(EXIT_SUCCESS);
  }

  // Check input parameters are well-defined
  if(!vm.count("grid")){
    cout << "grid parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };
  if(!vm.count("heel")){
    cout << "heel parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };
  if(!vm.count("toe")){
    cout << "toe parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };
  assert(vm.count("radius"));

  if(!vm.count("radius")){
    cout << "radius parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };

  if(vm.count("compdat")){
    if(!vm.count("radius")){
      cout << "well-name parameter missing..." << endl
           << usage_msg << endl;
      exit(EXIT_FAILURE);
    };
  };
  if(vm["heel"].as<vector<double>>().size() != 3){
    cout << "heel parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };
  if(vm["toe"].as<vector<double>>().size() != 3){
    cout << "heel parameter missing..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };

  if(vm["radius"].as<double>() <= 0.0){
    cout << "radius must be larger than zero..." << endl
         << usage_msg << endl;
    exit(EXIT_FAILURE);
  };

  // Checking that provided grid file actually exists
  // NoteMB: See definition of exists() at start of file
  // I assume boost::filesystem::exists is problematic in
  // Windows, and therefore OV has defined an equivalent
  // function based on <sys/stat.h>
  // Here we apply the function conditional on OS
  bool out;
#if _WIN32
  out = exists(vm["grid"].as<string>());
#else
  out = boost::filesystem::exists(vm["grid"].as<string>());
#endif
  if (!out){
    cout << "Grid file missing..." << endl;
    exit(EXIT_FAILURE);
  };

  return vm;
}

#endif // WIC_MAIN_H
