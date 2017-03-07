/******************************************************************************
   Copyright (C) 2015-2016 Hilmar M. Magnusson <hilmarmag@gmail.com>
   Modified by Einar J.M. Baumann (2016) <einar.baumann@gmail.com>

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

#include <iostream>
#include "wellindexcalculator.h"

namespace Reservoir {
namespace WellIndexCalculation {
using namespace std;

WellIndexCalculator::WellIndexCalculator(Grid::Grid *grid) {
    grid_ = grid;
}

vector<IntersectedCell>
WellIndexCalculator::ComputeWellBlocks(Vector3d heel,
                                       Vector3d toe,
                                       double wellbore_radius) {
    heel_ = heel;
    toe_ = toe;
    wellbore_radius_ = wellbore_radius;
    vector<IntersectedCell> intersected_cells = cells_intersected();

    for (int i = 0; i < intersected_cells.size(); ++i) {
        intersected_cells[i].set_well_index(compute_well_index(intersected_cells[i]));
    }
    return intersected_cells;
}

vector<IntersectedCell> WellIndexCalculator::cells_intersected() {
    vector<IntersectedCell> intersected_cells;

    // Find the heel cell and add it to the list
    intersected_cells.push_back(
        IntersectedCell( grid_->GetCellEnvelopingPoint(heel_) )
    );
    grid_->FillCellProperties(intersected_cells[0]); // Get poro and perm data
    intersected_cells[0].set_entry_point(heel_);

    // Find the toe cell -- removed by OV; why?
    Grid::Cell last_cell = grid_->GetCellEnvelopingPoint(toe_);

    // If the first and last blocks are the same, return the block and start+end points
    if (last_cell.global_index() == intersected_cells[0].global_index()) {
        intersected_cells[0].set_exit_point(toe_);
        return intersected_cells;
    }

    // Make sure we follow line in the correct direction. (i.e. dot product positive)
    Vector3d exit_point;
    exit_point = find_exit_point(intersected_cells[0], heel_, toe_, heel_);

    if ((toe_ - heel_).dot(exit_point - heel_) <= 0.0) {
        exit_point = find_exit_point(intersected_cells[0], heel_, toe_, exit_point);
    }
    intersected_cells[0].set_exit_point(exit_point);

    double epsilon = 0.02 / (toe_ - exit_point).norm();

    // Add previous exit point to list, find next exit point and all other up to the end_point
    while (true) {
        // Move into the next cell, add it to the list and set the entry point
        Vector3d move_exit_epsilon = exit_point * (1 - epsilon) + toe_ * epsilon;
        intersected_cells.push_back(
            IntersectedCell(grid_->GetCellEnvelopingPoint( move_exit_epsilon ))
        );
        grid_->FillCellProperties(intersected_cells.back()); // Get poro and perm data
        // The entry point of each cell is the exit point of the previous cell
        intersected_cells.back().set_entry_point(exit_point);

        // Terminate if we're in the last cell
        if (intersected_cells.back().global_index() == last_cell.global_index()) {
            intersected_cells.back().set_exit_point(toe_);
            break;
        }

        // Find the exit point of the cell and set it in the list
        exit_point = find_exit_point(intersected_cells.back(), exit_point, toe_, exit_point);
        intersected_cells.back().set_exit_point(exit_point);
        if (intersected_cells.size() > 500) {
            cout << "WARNING: More than 500 cells intersected. Returning empty array." << endl;
            return vector<IntersectedCell>();
        }
    }

    assert(intersected_cells.back().global_index() == last_cell.global_index());
    return intersected_cells;
}

Vector3d WellIndexCalculator::find_exit_point(Grid::Cell &cell,
                                              Vector3d &entry_point,
                                              Vector3d &end_point,
                                              Vector3d &exception_point) {
    Vector3d line = end_point - entry_point;

    // Loop through the cell faces until we find one that the line intersects
    for (Grid::Cell::Face face : cell.faces()) {
        if (face.normal_vector.dot(line) != 0) {
            // Check that the line and face are not parallel.
            auto intersect_point = face.intersection_with_line(entry_point, end_point);

            // Check that the intersect point is on the correct
            // side of all faces (i.e. inside the cell)
            bool feasible_point = true;
            for (auto p : cell.faces()) {
                if (!p.point_on_same_side(intersect_point, 10e-6)) {
                    feasible_point = false;
                    break;
                }
            }

            // Return intersection point if:
            // it is deemed feasible:
            auto cond_a = feasible_point;
            // + it is not identical to entry point
            auto cond_b = (exception_point - intersect_point).norm() > 10e-10;
            // + it is going in the correct direction.
            auto cond_c = (end_point - entry_point).dot(end_point - intersect_point) >= 0;
            if (cond_a && cond_b && cond_c) {
                return intersect_point;
            }
        }
    }
    // If all fails, the line intersects the cell in a
    // single point (corner or edge) -> return entry_point
    return entry_point;
}

double WellIndexCalculator::compute_well_index(IntersectedCell &icell) {
    double Lx = 0;
    double Ly = 0;
    double Lz = 0;

    for (int ii = 0; ii < icell.points().size() - 1; ++ii) { // Current segment ii
        // Compute vector from segment
        Vector3d current_vec = icell.points().at(ii+1) - icell.points().at(ii);

        /*
         * Projects segment vector to directional spanning vectors and determines the length.
         * of the projections. Note that we only only care about the length of the projection,
         * not the spatial position. Also adds the lengths of previous segments in case there
         * is more than one segment within the well.
         */
        Lx = Lx + (icell.xvec() * icell.xvec().dot(current_vec) / icell.xvec().dot(icell.xvec())).norm();
        Ly = Ly + (icell.yvec() * icell.yvec().dot(current_vec) / icell.yvec().dot(icell.yvec())).norm();
        Lz = Lz + (icell.zvec() * icell.zvec().dot(current_vec) / icell.zvec().dot(icell.zvec())).norm();
    }

    // Compute Well Index from formula provided by Shu (\todo Introduce ref/year)
    double well_index_x = (dir_well_index(Lx, icell.dy(), icell.dz(), icell.permy(), icell.permz()));
    double well_index_y = (dir_well_index(Ly, icell.dx(), icell.dz(), icell.permx(), icell.permz()));
    double well_index_z = (dir_well_index(Lz, icell.dx(), icell.dy(), icell.permx(), icell.permy()));

    double wi = sqrt(well_index_x * well_index_x +
        well_index_y * well_index_y +
        well_index_z * well_index_z);
    return wi;
}

double WellIndexCalculator::dir_well_index(double Lx, double dy,
                                           double dz, double ky, double kz) {
    double silly_eclipse_factor = 0.008527;
    double well_index_i = silly_eclipse_factor * (2 * M_PI * sqrt(ky * kz) * Lx) /
        (log(dir_wellblock_radius(dy, dz, ky, kz) / wellbore_radius_));
    return well_index_i;
}

double WellIndexCalculator::dir_wellblock_radius(double dx, double dy,
                                                 double kx, double ky) {
    double r = 0.28 * sqrt((dx * dx) * sqrt(ky / kx) + (dy * dy) * sqrt(kx / ky)) /
        (sqrt(sqrt(kx / ky)) + sqrt(sqrt(ky / kx)));
    return r;
}
}
}
