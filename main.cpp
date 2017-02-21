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
 * @brief This file contains the main function for the stand-alone
 * well index calculator executable.
 */

#include "main.hpp"
#include "wellindexcalculator.h"
#include <Reservoir/grid/eclgrid.h>

using namespace std;

int main(int argc, const char *argv[]) {
  // Initialize some variables from the runtime arguments
  Eigen::setNbThreads(1); // OV

  auto vm = createVariablesMap(argc, argv);
  auto heel = Vector3d(vm["heel"].as<vector<double>>().data());
  auto toe = Vector3d(vm["toe"].as<vector<double>>().data());
  string grid_path = vm["grid"].as<string>();
  double wellbore_radius = vm["radius"].as<double>();

  // Initialize the Grid and WellIndexCalculator objects
  auto grid = new Reservoir::Grid::ECLGrid(grid_path);
  auto wic = WellIndexCalculator(grid);

  // Compute well blocks
  auto well_blocks = wic.ComputeWellBlocks(heel, toe, wellbore_radius);

  // Print as a COMPDAT table if the --compdat/-c flag was given
  if (vm.count("compdat")) {
    string well_name = vm["well-name"].as<string>();
    printCompdat(well_blocks, well_name, wellbore_radius);
  }
    // Otherwise, print as a CSV table
  else {
    printCsv(well_blocks);
  }

  return 0;
}
