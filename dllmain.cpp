
#include "dllmain.hpp"
#include "wellindexcalculator.h"
#include "Reservoir/grid/eclgrid.h"
#include "intersected_cell.h"

#include <stdio.h>
#include <sys/stat.h>

#include <stdbool.h>
#include <iostream>
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
  printf("DllInitializer\n");
}

__attribute__((destructor))
/**
 * It is called when dylib is being unloaded.
 *
 */
static void Finalizer()
{
  printf("DllFinalizer\n");
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
    string gridpth = string(basepth) + ".GRID";
    if ( !exists(gridpth.c_str()) )
      throw std::runtime_error("ComputeWellIndices: file .GRID does not exist");

    if ( *wellbore_radius <= 0 )
      throw std::runtime_error("ComputeWellIndices: wellbore radius is negative");

    // Initialize the Grid and WellIndexCalculator objects
    if (grid == NULL)
      grid = new Reservoir::Grid::ECLGrid(gridpth);
    auto wic = WellIndexCalculator((Reservoir::Grid::Grid*)grid);

    // Compute the well blocks
    auto well_blocks = wic.ComputeWellBlocks(Eigen::Vector3d(heel), Eigen::Vector3d(toe), *wellbore_radius);
    *nblks = well_blocks.size();

    if ( (i == NULL) || (j == NULL) || (k == NULL) || (wi == NULL) )
      throw std::runtime_error("ComputeWellIndices: allocate memory for I, J, K, WI");

    for (size_t iblk = 0; iblk < *nblks; iblk++)
    {
      i[iblk] = well_blocks[iblk].ijk_index().i() + 1;
      j[iblk] = well_blocks[iblk].ijk_index().j() + 1;
      k[iblk] = well_blocks[iblk].ijk_index().k() + 1;
      wi[iblk] = well_blocks[iblk].well_index();
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

WELLINDEXCALCULATOR_API int computeBlockCenter(const char* basepth,
                                               const int* heel,  const int* toe, double* heelxyz, double* toexyz)
{
  try
  {
    string gridpth = string(basepth) + ".GRID";
    if ( !exists(gridpth.c_str()) )
      throw std::runtime_error("ComputeBlockCenter: file .GRID does not exist");

    // Initialize the Grid and WellIndexCalculator objects
    if (grid == NULL)
      grid = new Reservoir::Grid::ECLGrid(gridpth);

    Reservoir::Grid::Cell current_cell =
        ((Reservoir::Grid::Grid*)grid)->GetCell(heel[0] - 1, heel[1] - 1, heel[2] - 1);
    heelxyz[0] = current_cell.center()[0];
    heelxyz[1] = current_cell.center()[1];
    heelxyz[2] = current_cell.center()[2];

    current_cell = ((Reservoir::Grid::Grid*)grid)->GetCell(toe[0] - 1, toe[1] - 1, toe[2] - 1);
    toexyz[0] = current_cell.center()[0];
    toexyz[1] = current_cell.center()[1];
    toexyz[2] = current_cell.center()[2];

    //delete grid;
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

