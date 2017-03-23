/******************************************************************************
   Copyright (C) 2015-2016 Hilmar M. Magnusson <hilmarmag@gmail.com>
   Modified by Einar J.M. Baumann (2016) <einar.baumann@gmail.com>
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

#include "wellindexcalculator.h"
#include <fstream>
#include <map>
#include <iostream>

namespace Reservoir {
namespace WellIndexCalculation {
WellIndexCalculator::WellIndexCalculator(Grid::Grid *grid) {
    grid_ = grid;
}

std::map<std::string, std::vector<IntersectedCell>>
WellIndexCalculator::ComputeWellBlocks(std::vector<WellDefinition> wells)
{
    std::map<std::string, std::vector<IntersectedCell>> well_indices;
    for (int iWell = 0; iWell < wells.size(); ++iWell)
    {
        // Compute an overall bounding box per well --> speed up cell searching
        double xi, yi, zi, xf, yf, zf;
        xi = std::numeric_limits<double>::max();
        xf = std::numeric_limits<double>::min();
        yi = std::numeric_limits<double>::max();
        yf = std::numeric_limits<double>::min();
        zi = std::numeric_limits<double>::max();
        zf = std::numeric_limits<double>::min();

        for (int iSegment = 0; iSegment < wells[iWell].radii.size(); ++iSegment )
        {
            xi = std::min(xi, std::min(wells[iWell].heels[iSegment].x(), wells[iWell].toes[iSegment].x()));
            xf = std::max(xf, std::max(wells[iWell].heels[iSegment].x(), wells[iWell].toes[iSegment].x()));
            yi = std::min(yi, std::min(wells[iWell].heels[iSegment].y(), wells[iWell].toes[iSegment].y()));
            yf = std::max(yf, std::max(wells[iWell].heels[iSegment].y(), wells[iWell].toes[iSegment].y()));
            zi = std::min(zi, std::min(wells[iWell].heels[iSegment].z(), wells[iWell].toes[iSegment].z()));
            zf = std::max(zf, std::max(wells[iWell].heels[iSegment].z(), wells[iWell].toes[iSegment].z()));
        }

        std::cout << xi  << " " << xf << " " << yi << " " << yf << " " << zi << " " << zf << std::endl;
        // Artificially and heuristically increase the size of the searching area
        xi = xi - 0.1*(xf-xi); xf = xf + 0.1*(xf-xi);
        yi = yi - 0.1*(yf-yi); yf = yf + 0.1*(yf-yi);
        zi = zi - 0.1*(zf-zi); zf = zf + 0.1*(zf-zi);

        // Get the list of all cell in the bounding box
        std::vector<int> bb_cells;
        bb_cells = grid_->GetBoundingBoxCellIndices( xi, yi, zi, xf, yf, zf );
        std::cout << bb_cells.size() << std::endl;

        // Compute cells intersected by all well segments
        std::vector<IntersectedCell> intersected_cells;
        for (int iSegment = 0; iSegment < wells[iWell].radii.size(); ++iSegment )
        {
            collect_intersected_cells(intersected_cells, wells[iWell].heels[iSegment],
                                      wells[iWell].toes[iSegment], wells[iWell].radii[iSegment],
                                      bb_cells);
        }

        // std::cout << "number of intersected cells: " << intersected_cells.size() << std::endl;

        // For all intersected cells compute well transmissibility index
        for (int iCell = 0; iCell < intersected_cells.size(); ++iCell)
        {
            compute_well_index(intersected_cells, iCell);
        }

        well_indices[wells[iWell].wellname] = intersected_cells;
    }

    // Return the objects
    return well_indices;
}


void WellIndexCalculator::collect_intersected_cells(std::vector<IntersectedCell> &intersected_cells,
                                                    Vector3d start_point, Vector3d end_point, double wellbore_radius,
                                                    std::vector<int> bb_cells)
{
    // If no cells are found in the bounding box it means this segment is completely out of the reservoir
    if (bb_cells.size() == 0)
        return;

    // std::cout << "number of cells in the bb " << bb_cells.size() << std::endl;

    // Find the heel cell
    Grid::Cell first_cell;
    Vector3d org_start_point = start_point;
    double epsilon = 0.01 / (end_point - start_point).norm();
    double step = 0;

    while (true)
    {
        try
        {
            first_cell = grid_->GetCellEnvelopingPoint(start_point, bb_cells);
            //std::cout << "first cell found " << first_cell.global_index() << std::endl;
            break;
        }
        catch (const std::runtime_error& e)
        {
            step += epsilon;
            //std::cout << "The starting point is outside the grid: " << e.what() << std::endl;
            start_point = org_start_point * (1 - step) + end_point * step;

            std::cout << " " << step << " " << (org_start_point - end_point).dot(start_point - end_point) << std::endl;

            // Check if we went too far
            //if ((end_point - org_start_point).dot(start_point - org_start_point) <= 0.0)
            if (step > 1.0)
            {
                // The entire segment is outside the grid;
                std::cout << "_";
                return;
            }
        }
    }

    // Find the toe cell
    Grid::Cell last_cell;
    Vector3d org_end_point = end_point;
    step = 0;
    while (true)
    {
        try
        {
            last_cell = grid_->GetCellEnvelopingPoint(end_point, bb_cells);
            //std::cout << "last cell found " << last_cell.global_index() << std::endl;
            break;
        }
        catch (const std::runtime_error& e)
        {
            step += epsilon;
            //std::cout << "The end point is outside the grid: " << e.what() << std::endl;
            end_point = org_end_point*(1 - step) + start_point*step;
            std::cout << " " << step << std::endl;
            // Check if we went too far
            //if ((end_point - start_point).dot(org_end_point - start_point) <= 0.0)
            if (step > 1.0)
            {
                // The entire segment is outside the grid;
                std::cout << ".";
                return;
            }
        }
    }

    // If the first and last blocks are the same, return the block and start+end points
    if (last_cell.global_index() == first_cell.global_index())
    {
        // Get the index of the interesected cell that corresponds to the cell where
        // the first point resides (if this is not yet in the list it will be added)
        int intersected_cell_index = IntersectedCell::GetIntersectedCellIndex(intersected_cells, first_cell);

        intersected_cells.at(intersected_cell_index).add_new_segment(start_point, end_point, wellbore_radius);
        return;
    }

    // Get the index of the interesected cell that corresponds to the cell where
    // the first point resides (if this is not yet in the list it will be added)
    int intersected_cell_index = IntersectedCell::GetIntersectedCellIndex(intersected_cells, first_cell);

    // Make sure we follow line in the correct direction. (i.e. dot product positive)
    Vector3d exit_point = find_exit_point(intersected_cells, intersected_cell_index, start_point, end_point, start_point);
    if ((end_point - start_point).dot(exit_point - start_point) <= 0.0)
    {
        exit_point = find_exit_point(intersected_cells, intersected_cell_index, start_point, end_point, exit_point);
    }
    intersected_cells.at(intersected_cell_index).add_new_segment(start_point, exit_point, wellbore_radius);

    // The step is different here
    epsilon = 0.01 / (end_point - exit_point).norm();
    step = 0;

    // Add previous exit point to list, find next exit point and all other up to the end_point
    while (true)
    {
        // Move into the next cell, add it to the list and set the entry point
        Vector3d move_exit_epsilon = exit_point * (1 - epsilon) + end_point * epsilon;

        // Find the next cell, but this might be inactive
        Grid::Cell new_cell;
        while (true)
        {
            step = epsilon;
            try
            {
                new_cell = grid_->GetCellEnvelopingPoint(move_exit_epsilon, bb_cells);
                break;
            }
            catch( const std::runtime_error& e)
            {
                step += epsilon;
                move_exit_epsilon = exit_point * (1 - step) + end_point * step;
                std::cout << " " << step << std::endl;
                // Check if we are not too far
                //if ((exit_point - end_point).dot(move_exit_epsilon - end_point) <= 0.0)
                if (step > 1.0)
                {
                    std::cout << "*";
                    return;
                }
            }
        }
        intersected_cell_index = IntersectedCell::GetIntersectedCellIndex(intersected_cells, new_cell);

        // Stop if we're in the last cell
        if (intersected_cells.at(intersected_cell_index).global_index() == last_cell.global_index())
        {
            // The entry point of each cell is the exit point of the previous cell
            intersected_cells.at(intersected_cell_index).add_new_segment(exit_point, end_point, wellbore_radius);
            break;
        }

        // Find the exit point of the cell and set it in the list
        Vector3d entry_point = exit_point;
        exit_point = find_exit_point(intersected_cells, intersected_cell_index, entry_point, end_point, exit_point);
        intersected_cells.at(intersected_cell_index).add_new_segment(entry_point, exit_point, wellbore_radius);

        // This should be removed
        // assert(intersected_cells.size() < 500);
    }

    assert(intersected_cells.at(intersected_cell_index).global_index() == last_cell.global_index());
}

Vector3d WellIndexCalculator::find_exit_point(std::vector<IntersectedCell> &cells, int cell_index,
                                              Vector3d &entry_point, Vector3d &end_point, Vector3d &exception_point)
{
    Vector3d line = end_point - entry_point;

    // Loop through the cell faces until we find one that the line intersects
    for (Grid::Cell::Face face : cells.at(cell_index).faces())
    {
        if (face.normal_vector.dot(line) != 0)
        { // Check that the line and face are not parallel.
            auto intersect_point = face.intersection_with_line(entry_point, end_point);

            // Check that the intersect point is on the correct side of all faces (i.e. inside the cell)
            bool feasible_point = true;
            for (auto p : cells.at(cell_index).faces())
            {
                if (!p.point_on_same_side(intersect_point, 10e-6))
                {
                    feasible_point = false;
                    break;
                }
            }

            // Return the point if it is deemed feasible, not identical to the entry point, and going in the correct direction.
            if (feasible_point && (exception_point - intersect_point).norm() > 10e-10
                && (end_point - entry_point).dot(end_point - intersect_point) >= 0)
            {
                return intersect_point;
            }
        }
    }

    // If all fails, the line intersects the cell in a single point (corner or edge) -> return entry_point
    return entry_point;
}

void WellIndexCalculator::compute_well_index(std::vector<IntersectedCell> &cells, int cell_index) {
    double well_index_x = 0;
    double well_index_y = 0;
    double well_index_z = 0;

    IntersectedCell &icell = cells.at(cell_index);

    for (int iSegment = 0; iSegment < icell.num_segments(); iSegment++){
        // Compute vector from segment
        Vector3d current_vec = icell.get_segment_exit_point(iSegment) - icell.get_segment_entry_point(iSegment);

        /* Projects segment vector to directional spanning vectors and determines the length.
         * of the projections. Note that we only care about the length of the projection,
         * not the spatial position. Also adds the lengths of previous segments in case there
         * is more than one segment within the well.
         */
        double current_Lx = (icell.xvec() * icell.xvec().dot(current_vec) / icell.xvec().dot(icell.xvec())).norm();
        double current_Ly = (icell.yvec() * icell.yvec().dot(current_vec) / icell.yvec().dot(icell.yvec())).norm();
        double current_Lz = (icell.zvec() * icell.zvec().dot(current_vec) / icell.zvec().dot(icell.zvec())).norm();

        // Compute Well Index from formula provided by Shu per Segment
        // (Note that this has a glich since segments from the same well could have different radius (e.g. radial well))
        double current_wx = dir_well_index(current_Lx, icell.dy(), icell.dz(), icell.permy(), icell.permz(), icell.get_segment_radius(iSegment));
        double current_wy = dir_well_index(current_Ly, icell.dx(), icell.dz(), icell.permx(), icell.permz(), icell.get_segment_radius(iSegment));
        double current_wz = dir_well_index(current_Lz, icell.dx(), icell.dy(), icell.permx(), icell.permy(), icell.get_segment_radius(iSegment));

        // Compute the sum of well index for each direction.
        // For segments with equal radius this will in the end calculate the well index based on the Shu formula in its original formulation
        well_index_x += current_wx;
        well_index_y += current_wy;
        well_index_z += current_wz;

        // Store data for later use
        icell.set_segment_calculation_data(iSegment, "dx", icell.dx());
        icell.set_segment_calculation_data(iSegment, "dy", icell.dy());
        icell.set_segment_calculation_data(iSegment, "dz", icell.dz());

        icell.set_segment_calculation_data(iSegment, "permx", icell.permx());
        icell.set_segment_calculation_data(iSegment, "permy", icell.permy());
        icell.set_segment_calculation_data(iSegment, "permz", icell.permz());

        icell.set_segment_calculation_data(iSegment, "Lx", current_Lx);
        icell.set_segment_calculation_data(iSegment, "Ly", current_Ly);
        icell.set_segment_calculation_data(iSegment, "Lz", current_Lz);

        icell.set_segment_calculation_data(iSegment, "wx", current_wx);
        icell.set_segment_calculation_data(iSegment, "wy", current_wy);
        icell.set_segment_calculation_data(iSegment, "wz", current_wz);
    }

    // Compute the combined well index as the Sum the segment
    // Compute Well Index from formula provided by Shu for the entire combined projections (this is the original formulation)
    icell.set_cell_well_index(sqrt(well_index_x * well_index_x + well_index_y * well_index_y + well_index_z * well_index_z));
}

double WellIndexCalculator::dir_well_index(double Lx, double dy, double dz, double ky, double kz, double wellbore_radius) {
    double silly_eclipse_factor = 0.008527;
    double well_index_i = silly_eclipse_factor * (2 * M_PI * sqrt(ky * kz) * Lx) /
        (log(dir_wellblock_radius(dy, dz, ky, kz) / wellbore_radius));
    return well_index_i;
}

double WellIndexCalculator::dir_wellblock_radius(double dx, double dy, double kx, double ky)
{
    double r = 0.28 * sqrt((dx * dx) * sqrt(ky / kx) + (dy * dy) * sqrt(kx / ky)) /
        (sqrt(sqrt(kx / ky)) + sqrt(sqrt(ky / kx)));
    return r;
}

void WellDefinition::ReadWellsFromFile(std::string file_path, std::vector<WellDefinition>& wells)
{
    std::ifstream infile(file_path);
    std::string previous_well_name = "";
    std::string well_name;
    double hx, hy, hz, tx, ty, tz;
    double radius;
    //try
    {
        while (infile >> well_name >> hx >> hy >> hz >> tx >> ty >> tz >> radius)
        {
            if (previous_well_name != well_name)
            {
                wells.push_back(WellDefinition());
                wells.back().wellname = well_name;

                previous_well_name = well_name;
            }

            wells.back().heels.push_back(Eigen::Vector3d(hx,hy,hz));
            wells.back().toes.push_back(Eigen::Vector3d(tx,ty,tz));
            wells.back().radii.push_back(radius);
        }
    }

}
}
}
