/***********************************************************
 Copyright (C) 2017
 Oleg Volkov <ovolkov@stanford.edu>
 Mathias C. Bellout <mathias.bellout@ntnu.no>

 Created by bellout on 22/3/18.

 This file is part of the FieldOpt project.

 FieldOpt is free software: you can redistribute it
 and/or modify it under the terms of the GNU General
 Public License as published by the Free Software
 Foundation, either version 3 of the License, or (at
 your option) any later version.

 FieldOpt is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the
 implied warranty of MERCHANTABILITY or FITNESS FOR
 A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

 You should have received a copy of the GNU
 General Public License along with FieldOpt.
 If not, see <http://www.gnu.org/licenses/>.
***********************************************************/

// ---------------------------------------------------------
// STD
#include <stdio.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

// ---------------------------------------------------------
// FieldOpt: WIC
#include "dllmain.hpp"
#include "main.hpp"
#include "wellindexcalculator.h"
#include "Reservoir/grid/eclgrid.h"
#include "intersected_cell.h"
#include "wicalc_rixx.h"

// ---------------------------------------------------------
using std::runtime_error;


// =========================================================
inline bool exists(const char* name)
{
  struct stat buffer;
  return (stat (name, &buffer) == 0);
}

using namespace Reservoir::WellIndexCalculation;

// ---------------------------------------------------------
#if _WIN32
BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
  switch (ul_reason_for_call)
  {
  case DLL_PROCESS_ATTACH:
    grid = NULL;
    if (lpReserved == NULL) // Dynamic load
    {
      // Initialize your stuff or whatever
      // Return FALSE if you don't want your module to be dynamically loaded
    }
    else // Static load
    {
      // Return FALSE if you don't want your module to be statically loaded
      return FALSE;
    }
    break;
  case DLL_THREAD_ATTACH:
  case DLL_THREAD_DETACH:
  case DLL_PROCESS_DETACH:
    break;
  }
  return TRUE;
}
#else

// =========================================================
//CP_BEGIN_EXTERN_C
__attribute__((constructor))
/**
 * initializer of the dylib.
 */
static void Initializer(int argc, char** argv, char** envp)
{
  grid = NULL;
  printf("DllInitializer: Loading WIClib\n");
}

// =========================================================
__attribute__((destructor))
/**
 * It is called when dylib is being unloaded.
 *
 */
static void Finalizer()
{
  printf("DllFinalizer: Loaded WIClib\n");
}

//CP_END_EXTERN_C
#endif

// ---------------------------------------------------------
// This is an example of an exported variable
//WELLINDEXCALCULATOR_API Grid::Grid grid = NULL;

// =========================================================
// This is an example of an exported function.
WELLINDEXCALCULATOR_API int
computeWellIndices(const char* basepth,
                   const double* heel,
                   const double* toe,
                   const double* wellbore_radius,
                   int* nblks,
                   int* i, int* j, int* k,
                   double* wi) {

  // -------------------------------------------------------
  try {

    // -----------------------------------------------------
    string gridpth = string(basepth) + ".EGRID";
    printf("Loading grid: %s\n", gridpth.c_str());

    // -----------------------------------------------------
    if ( !exists(gridpth.c_str()) ) {
      throw runtime_error("ComputeWellIndices: "
                              "file .EGRID does not exist");
    }

    // -----------------------------------------------------
    if ( *wellbore_radius <= 0 ) {
      throw runtime_error("ComputeWellIndices: "
                              "wellbore radius is negative");
    }

    // -----------------------------------------------------
    printf("Initializing Grid object.\n");
    Reservoir::Grid::ECLGrid gridnew(gridpth);

    // NOTE: Make this work so we don't have to
    // initialize a new grid for every well!
    /*if (grid == NULL)
    {
      grid = new Reservoir::Grid::ECLGrid(gridpth);
      cout << "Making new grid\n";
    }*/

    // -----------------------------------------------------
    printf("Initializing WellIndexCalculator object.\n");
    // Old WIC
    // auto wic = WellIndexCalculator((Reservoir::Grid::Grid*)grid);
    // auto wic = WellIndexCalculator(&gridnew);

    // -----------------------------------------------------
    // New WIC
    Settings::Model::Well well_settings_;
    well_settings_.name = "DEFWELL";
    well_settings_.wellbore_radius = *wellbore_radius;
    well_settings_.verb_vector_ = std::vector<int>(11,0);
    Reservoir::WellIndexCalculation::wicalc_rixx wicalc_rixx =
        Reservoir::WellIndexCalculation::wicalc_rixx(well_settings_,
                                                     &gridnew);

    // -----------------------------------------------------
    printf("Setting up well.\n");
    vector<WellDefinition> wells;
    wells.push_back(WellDefinition());
    wells.at(0).wellname = well_settings_.name.toStdString();
    wells.at(0).radii.push_back(well_settings_.wellbore_radius);
    wells.at(0).skins.push_back(0.0);

    auto heelV3d = Eigen::Vector3d(heel);
    auto toeV3d = Eigen::Vector3d(toe);
    wells.at(0).heels.push_back(heelV3d);
    wells.at(0).toes.push_back(toeV3d);

    wells.at(0).well_length.push_back(sqrt((toeV3d - heelV3d).norm()));
    wells.at(0).heel_md.push_back(heelV3d(2));
    wells.at(0).toe_md.push_back(wells.at(0).heel_md.back() + wells.at(0).well_length.back());

    // -----------------------------------------------------
    printf("Computing well blocks.\n");
    map<string, vector<IntersectedCell>> well_indices;
    // Old WIC
    // wic.ComputeWellBlocks(well_indices, wells);
    // New WIC
    wicalc_rixx.ComputeWellBlocks(well_indices, wells);

    // -----------------------------------------------------
    // printCompdat(well_indices);
    // Remember to remove blocks with very
    // low WIs to ev. reduce sim problems?
    vector<IntersectedCell>& well_blocks = well_indices.at("DEFWELL");
    *nblks = well_blocks.size();

    // -----------------------------------------------------
    printf("%s", "Well block check 1: ");
    if ( (i == NULL) || (j == NULL) || (k == NULL) || (wi == NULL) )
      throw runtime_error("ComputeWellIndices: I, J, K, WI not allocated");

    // -----------------------------------------------------
    printf("%s", "Well block check 2: ");
    try {

      for (size_t iblk = 0; iblk < *nblks; iblk++) {

        i[iblk] = well_blocks[iblk].ijk_index().i() + 1;
        j[iblk] = well_blocks[iblk].ijk_index().j() + 1;
        k[iblk] = well_blocks[iblk].ijk_index().k() + 1;
        wi[iblk] = well_blocks[iblk].cell_well_index_matrix();
      }

    } catch (...) {
      throw runtime_error("ComputeWellIndices: "
                              "problem with copying I, J, K, WI");
    }

    // -----------------------------------------------------
    std::stringstream str;
    str << "\x1b[33m" << "# of blocks: ";
    str << well_blocks.size() << ". \x1b[0m";
    printf("%s\n", str.str().c_str());
    //delete grid;

  } catch (std::runtime_error& e) {

    fprintf(stdout, "Failure: %s", e.what());
    fflush(stdout);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// =========================================================
WELLINDEXCALCULATOR_API int getBlockCenters(const char* basepth,
                                            const int* heel,
                                            const int* toe,
                                            double* heelxyz,
                                            double* toexyz) {

  try {

    // -------------------------------------------------------------
    string gridpth = string(basepth) + ".EGRID";
    if ( !exists(gridpth.c_str()) )
      throw runtime_error("getBlockCenters: file .GRID does not exist");

    // -------------------------------------------------------------
    // Initialize the Grid and WellIndexCalculator objects
    Reservoir::Grid::ECLGrid gridnew(gridpth);
    /*if (grid == NULL)
    {
      grid = new Reservoir::Grid::ECLGrid(gridpth);
      cout << "Making new grid\n";
    }*/

    // -------------------------------------------------------------
    if ( (heel == NULL) || (toe == NULL) || (heelxyz == NULL) || (toexyz == NULL) )
      throw runtime_error("getBlockCenters: input and output not allocated");

    // -------------------------------------------------------------
    try {

      /*Reservoir::Grid::Cell current_cell =
        ((Reservoir::Grid::Grid*)grid)->GetCell(heel[0] - 1, heel[1] - 1, heel[2] - 1);*/

      Reservoir::Grid::Cell current_cell;
      current_cell =
          gridnew.GetCell(heel[0] - 1, heel[1] - 1, heel[2] - 1);
      heelxyz[0] = current_cell.center()[0];
      heelxyz[1] = current_cell.center()[1];
      heelxyz[2] = current_cell.center()[2];

      //current_cell = ((Reservoir::Grid::Grid*)grid)->GetCell(toe[0] - 1, toe[1] - 1, toe[2] - 1);
      current_cell =
          gridnew.GetCell(toe[0] - 1, toe[1] - 1, toe[2] - 1);
      toexyz[0] = current_cell.center()[0];
      toexyz[1] = current_cell.center()[1];
      toexyz[2] = current_cell.center()[2];
    }
    catch (...) {
      throw runtime_error("getBlockCenters: error in computing X,Y,Z");
    }
  }
  catch (runtime_error& e)
  {
    printf("%s", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// -----------------------------------------------------------------
WELLINDEXCALCULATOR_API int getBoundaryVertices(const char* filepth,
                                                int* npnts,
                                                double* xes,
                                                double* yes,
                                                double* zes) {

  try {

    // -------------------------------------------------------------
    string bndrypth = string(filepth);
    if ( !exists(bndrypth.c_str()) )
      throw runtime_error("getBoundaryVertices: file does not exist\n");

    // -------------------------------------------------------------
    try {

      ifstream fbndry(bndrypth, ios::in);
      *npnts = 0;
      while ( fbndry >> xes[*npnts] )
      {
        fbndry >> yes[*npnts] >> zes[*npnts];
        *npnts = *npnts + 1;
      }
    }
    catch (...) {
      throw runtime_error("getBoundaryVertices: error in reading file\n");
    }

  }
  // ---------------------------------------------------------------
  catch (std::runtime_error& e) {
    printf("%s", e.what());
    return EXIT_FAILURE;
  }

  // ---------------------------------------------------------------
  return EXIT_SUCCESS;
}

// This is the constructor of a class that has been exported.
// see WellIndexCalculator.h for the class definition
//CWellIndexCalculator::CWellIndexCalculator()
//{
//  return;
//}

