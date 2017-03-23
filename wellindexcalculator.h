/******************************************************************************
   Copyright (C) 2015-2016 Hilmar M. Magnusson <hilmarmag@gmail.com>
   Modified by Einar J.M. Baumann (2016) <einar.baumann@gmail.com>
   Modified by Alin G. Chitu (2016-2017) <alin.chitu@tno.nl, chitu_alin@yahoo.com>
   Modified by Einar J.M. Baumann (2017) <einar.baumann@gmail.com>

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

#ifndef WELLINDEXCALCULATOR_H
#define WELLINDEXCALCULATOR_H

#include <Eigen/Dense>
#include <vector>
#include <Eigen/Core>
#include "Reservoir/grid/grid.h"
#include "intersected_cell.h"

namespace Reservoir {
namespace WellIndexCalculation {
using namespace Eigen;
using namespace std;

class WellDefinition {
 public:
  std::string wellname;
  std::vector<Vector3d> heels;
  std::vector<Vector3d> toes;
  std::vector<double> radii;

 public:
  static void ReadWellsFromFile(std::string file_path, std::vector<WellDefinition>& wells);
};

/*!
 * \brief The WellIndexCalculation class deduces the well blocks
 * and their respective well indices/transmissibility factors for
 * one or more well splines defined by a heel and a toe.
 *
 * Note that some of the internal datastructures in this class seem
 * more complex than they need to be. This is because the internal
 * methods support well splines consisting of more than one point.
 * This is, however, not yet supported by the Model library and so
 * have been "hidden".
 *
 * Credit for computations in this class goes to @hilmarm.
 */
class WellIndexCalculator {
 public:
  WellIndexCalculator(){}
  WellIndexCalculator(Grid::Grid *grid);

  /*!
   * \brief Compute the well block indices for all wells
   * \param wells The list of wells
   * \return A map containing for each well given my its name the list of cells intersected by the well.
   * Each intersected cell has stored the well connectivity information.
   */
  std::map<std::string, std::vector<IntersectedCell>> ComputeWellBlocks(std::vector<WellDefinition> wells);


 private:
  /*!
   * \brief The Well struct holds the information needed to compute
   * the well blocks and their respective well indices for a well
   * spline consisting of a heel and a toe.
   */

  Grid::Grid *grid_; //!< The grid used in the calculations.

 public:
  /*!
   * \brief Given a reservoir with blocks and a line (start_point
   * to end_point), return global index of all blocks interesected
   * by the line, as well as the point where the line enters the
   * block
   *
   * ?? by the line and the points of intersection
   *
   * \param intersected_cells Vector in which to store the cells
   * \param start_point The start point of the well path.
   * \param end_point The end point of the well path.
   * \param grid The grid object containing blocks/cells.
   * \param bb_cells
   *
   * \return A pair containing global indices of intersected
   * cells and the points where it enters each cell (and thereby
   * leaves the previous cell) of the line segment inside each
   * cell.
   */
  void collect_intersected_cells(std::vector<IntersectedCell> &intersected_cells,
                                 Vector3d start_point, Vector3d end_point, double wellbore_radius,
                                 std::vector<int> bb_cells);

  /*!
   * \brief Find the point where the line between the start_point and end_point exits a cell.
   *
   * Takes as input an entry_point end_point which defines the well
   * path. Finds the two points on the path which intersects the
   * block faces and chooses the one that is not the entry point,
   * i.e. the exit point.
   *
   * \todo Find a better name for the exception_point and describe it better.
   *
   * \param cell The cell to find the well paths exit point in.
   * \param start_point The start point of the well path.
   * \param end_point The end point of the well path.
   * \param exception_point A specific point we don't
   * want the function to end up in.
   *
   * \return The point where the well path exits the cell.
   */
  Vector3d find_exit_point(std::vector<IntersectedCell> &cells, int cell_index,
                           Vector3d &start_point, Vector3d &end_point, Vector3d &exception_point);

  /*!
   * \brief Compute the well index (aka. transmissibility factor)
   * for a (one) single cell/block by using the Projection Well
   * Method (Shu 2005).
   *
   * Assumption: The block is fairly regular,
   * i.e. corners are straight angles.
   *
   * \note Corner points of Cell(s) are always listed in the same
   * order and orientation. (see Grid::Cell for illustration).
   *
   * \param icell Well block to compute the WI in.
   * \return Well index for block/cell
  */
  void compute_well_index(std::vector<IntersectedCell> &cells, int cell_index);

  /*!
   * \brief Auxilary function for compute_well_index function
   *
   * \param Lx lenght of projection in first direction
   * \param dy size block second direction
   * \param dz size block third direction
   * \param ky permeability second direction
   * \param kz permeability second direction
   *
   * \return directional well index
  */
  double dir_well_index(double Lx, double dy, double dz, double ky, double kz, double wellbore_radius);

  /*!
   * \brief Auxilary function(2) for compute_well_index function
   *
   * \param dx size block second direction
   * \param dy size block third direction
   * \param kx permeability second direction
   * \param ky permeability second direction
   *
   * \return directional wellblock radius
   */
  double dir_wellblock_radius(double dx, double dy,
                              double kx, double ky);
};
}
}

#endif // WELLINDEXCALCULATOR_H
