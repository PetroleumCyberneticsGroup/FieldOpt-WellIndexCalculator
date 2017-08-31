
#include "dllmain.hpp"
#include "main.hpp"
#include "wellindexcalculator.h"
#include "Reservoir/grid/eclgrid.h"
#include "intersected_cell.h"

#include <stdio.h>
#include <sys/stat.h>

#include <stdbool.h>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

inline bool exists(const char* name)
{
  struct stat buffer;
  return (stat (name, &buffer) == 0);
}

using namespace Reservoir::WellIndexCalculation;

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

// This is an example of an exported variable
//WELLINDEXCALCULATOR_API Grid::Grid grid = NULL;

// This is an example of an exported function.
WELLINDEXCALCULATOR_API int computeWellIndices(const char* basepth,
    const double* heel, const double* toe, const double* wellbore_radius,
    int* nblks, int* i, int* j, int* k, double* wi)
{
  try
  {
    string gridpth = string(basepth) + ".EGRID";
    printf("Loading grid: %s\n", gridpth.c_str());

    if ( !exists(gridpth.c_str()) ) {
      throw std::runtime_error("ComputeWellIndices: file .EGRID does not exist");
    }

    if ( *wellbore_radius <= 0 ){
      throw std::runtime_error("ComputeWellIndices: wellbore radius is negative");
    }

    printf("Initializing Grid object.\n");
    Reservoir::Grid::ECLGrid gridnew(gridpth);
    /*if (grid == NULL)
    {
      grid = new Reservoir::Grid::ECLGrid(gridpth);
      cout << "Making new grid\n";
    }*/
    
    printf("Initializing WellIndexCalculator object.\n");
    // auto wic = WellIndexCalculator((Reservoir::Grid::Grid*)grid);
    auto wic = WellIndexCalculator(&gridnew);

    printf("Compute well blocks.\n");
    vector<WellDefinition> wells;
    wells.push_back(WellDefinition());
    wells.at(0).wellname = "unnamed_well";
    wells.at(0).heels.push_back(Eigen::Vector3d(heel));
    wells.at(0).toes.push_back(Eigen::Vector3d(toe));
    wells.at(0).radii.push_back(*wellbore_radius);
    wells.at(0).skins.push_back(0.0);

    map<string, vector<IntersectedCell>> well_indices;
    wic.ComputeWellBlocks(well_indices, wells);
    //printCompdat(well_indices);
    vector<IntersectedCell>& well_blocks = well_indices.at("unnamed_well");
    *nblks = well_blocks.size();

    if ( (i == NULL) || (j == NULL) || (k == NULL) || (wi == NULL) )
      throw std::runtime_error("ComputeWellIndices: I, J, K, WI not allocated");

    try
    {
      for (size_t iblk = 0; iblk < *nblks; iblk++)
      {
        i[iblk] = well_blocks[iblk].ijk_index().i() + 1;
        j[iblk] = well_blocks[iblk].ijk_index().j() + 1;
        k[iblk] = well_blocks[iblk].ijk_index().k() + 1;
        wi[iblk] = well_blocks[iblk].cell_well_index_matrix();
      }
    }
    catch (...)
    {
      throw std::runtime_error("ComputeWellIndices: problem with copying I, J, K, WI");
    }
    //delete grid;
  }
  catch (std::runtime_error& e)
  {
    fprintf(stdout, "%s", e.what());
    fflush(stdout);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

WELLINDEXCALCULATOR_API int getBlockCenters(const char* basepth,
    const int* heel,  const int* toe, double* heelxyz, double* toexyz)
{
  try
  {
    string gridpth = string(basepth) + ".EGRID";
    if ( !exists(gridpth.c_str()) )
      throw std::runtime_error("getBlockCenters: file .GRID does not exist");

    // Initialize the Grid and WellIndexCalculator objects
    Reservoir::Grid::ECLGrid gridnew(gridpth);
    /*if (grid == NULL)
    {
      grid = new Reservoir::Grid::ECLGrid(gridpth);
      cout << "Making new grid\n";
    }*/

    if ( (heel == NULL) || (toe == NULL) || (heelxyz == NULL) || (toexyz == NULL) )
      throw std::runtime_error("getBlockCenters: input and output not allocated");

    try
    {
      /*Reservoir::Grid::Cell current_cell =
        ((Reservoir::Grid::Grid*)grid)->GetCell(heel[0] - 1, heel[1] - 1, heel[2] - 1);*/
      Reservoir::Grid::Cell current_cell = gridnew.GetCell(heel[0] - 1, heel[1] - 1, heel[2] - 1);
      heelxyz[0] = current_cell.center()[0];
      heelxyz[1] = current_cell.center()[1];
      heelxyz[2] = current_cell.center()[2];

      //current_cell = ((Reservoir::Grid::Grid*)grid)->GetCell(toe[0] - 1, toe[1] - 1, toe[2] - 1);
      current_cell = gridnew.GetCell(toe[0] - 1, toe[1] - 1, toe[2] - 1);
      toexyz[0] = current_cell.center()[0];
      toexyz[1] = current_cell.center()[1];
      toexyz[2] = current_cell.center()[2];
    }
    catch (...)
    {
      throw std::runtime_error("getBlockCenters: error in computing X,Y,Z");
    }
  }
  catch (std::runtime_error& e)
  {
    printf("%s", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


WELLINDEXCALCULATOR_API int getBoundaryVertices(const char* filepth,
    int* npnts, double* xes, double* yes, double* zes)
{
  try
  {
    string bndrypth = string(filepth);
    if ( !exists(bndrypth.c_str()) )
      throw std::runtime_error("getBoundaryVertices: file does not exist");

    try
    {
      ifstream fbndry(bndrypth, ios::in);
      *npnts = 0;
      while ( fbndry >> xes[*npnts] )
      {
        fbndry >> yes[*npnts] >> zes[*npnts];
        *npnts = *npnts + 1;
      }
    }
    catch (...)
    {
      throw std::runtime_error("getBoundaryVertices: error in reading file");
    }

  }
  catch (std::runtime_error& e)
  {
    printf("%s", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// This is the constructor of a class that has been exported.
// see WellIndexCalculator.h for the class definition
//CWellIndexCalculator::CWellIndexCalculator()
//{
//  return;
//}

