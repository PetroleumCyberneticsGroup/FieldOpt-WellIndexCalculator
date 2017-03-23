/******************************************************************************
   Copyright (C) 2015-2016 Einar J.M. Baumann <einar.baumann@gmail.com>
   Modified by Alin G. Chitu (2016-2017) <alin.chitu@tno.nl, chitu_alin@yahoo.com>

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

using namespace Reservoir::WellIndexCalculation;
using namespace std;

int main(int argc, const char *argv[])
{
	// Initialize some variables from the runtime arguments
	Eigen::setNbThreads(1); // OV //AGC not sure what this does

	auto vm = createVariablesMap(argc, argv);

	// Checking that provided grid file actually exists
	// NoteMB: See definition of exists() at start of main.hpp file
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

	if (!out)
	{
		cout << "Grid file missing..." << endl;
		exit(EXIT_FAILURE);
	};

	// Make sure that the user either specifies an input file or the data for a fingle well segment is specified correctly
	assert(vm.count("well-filedef")	||(
		vm.count("heel") && vm.count("toe") &&
			vm["heel"].as<vector<double>>().size() == 3 &&
			vm["toe"].as<vector<double>>().size() == 3 &&
			vm.count("radius") && vm["radius"].as<double>() > 0 &&
			(vm.count("compdat")? (vm.count("well-name")? true:false):true)
	)
	);

	// Get the path to the grid file
	string gridpth = vm["grid"].as<string>();

	// Initialize the Grid and WellIndexCalculator objects
	Reservoir::Grid::ECLGrid* grid;
	try
	{
		grid = new Reservoir::Grid::ECLGrid(gridpth);
	}
	catch (const std::runtime_error& e)
	{
		std::cout << "Error reading the Eclipse grid " << e.what();
		std::cout << std::endl << "The program will stop now";

		exit(EXIT_FAILURE);
	}

	auto wic = WellIndexCalculator(grid);
	vector<WellDefinition> wells;

	if (vm.count("well-filedef") == 1)
	{
		assert(boost::filesystem::exists(vm["well-filedef"].as<string>()));
		WellDefinition::ReadWellsFromFile(vm["well-filedef"].as<string>(), wells);
	}
	else
	{
		wells.push_back(WellDefinition());
		if (vm.count("well-name"))
		{
			wells.at(0).wellname = vm["well-name"].as<string>();
		}
		wells.at(0).heels.push_back(Eigen::Vector3d(vm["heel"].as<vector<double>>().data()));
		wells.at(0).toes.push_back(Eigen::Vector3d(vm["toe"].as<vector<double>>().data()));
		wells.at(0).radii.push_back(vm["radius"].as<double>());
	}

	// Compute the well blocks
	auto well_indices = wic.ComputeWellBlocks(wells);

	if (vm.count("compdat")) { // Print as a COMPDAT table if the --compdat/-c flag was given
		printCompdat(well_indices);
	}
	else { // Otherwise, print as a CSV table
		printCsv(well_indices);
	}

	if (vm.count("debug") > 0)
	{
		printDebug(well_indices);
	}

	exit(EXIT_SUCCESS);
}
