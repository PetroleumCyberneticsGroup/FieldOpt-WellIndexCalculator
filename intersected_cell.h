/******************************************************************************
   Copyright (C) 2015-2016 Hilmar M. Magnusson <hilmarmag@gmail.com>
   Modified by Einar J.M. Baumann (2016) <einar.baumann@gmail.com>
   Modified by Alin G. Chitu (2016) <alin.chitu@tno.nl, chitu_alin@yahoo.com>

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

#ifndef FIELDOPT_INTERSECTEDCELL_H
#define FIELDOPT_INTERSECTEDCELL_H

#include "Reservoir/grid/cell.h"
#include <map>
#include <vector>

namespace Reservoir {
namespace WellIndexCalculation {
    using namespace Eigen;
    /*!
     * \brief The IntersectedCell struct holds information about an intersected cell.
     *
     */
    class IntersectedCell : public Grid::Cell {
    public:
        IntersectedCell() {}
        IntersectedCell(const Grid::Cell &cell) : Grid::Cell(cell) {};

        /*!
         * \brief The cell x axis
         */
        Vector3d xvec() const;
        /*!
         * \brief The cell y axis
         */
        Vector3d yvec() const;
        /*!
         * \brief The cell z axis
         */
        Vector3d zvec() const;

        // Cell size
        double dx() const;
        double dy() const;
        double dz() const;

        void add_new_segment(Vector3d entry_point, Vector3d exit_point, double segment_radius);
        int num_segments() const;

        Vector3d get_segment_entry_point(int segment_index) const;
        Vector3d get_segment_exit_point(int segment_index) const;
        double get_segment_radius(int segment_index) const;

        double cell_well_index() const;
        void set_cell_well_index(double well_index);

        void set_segment_calculation_data(int segment_index, std::string name, double value);
        std::map<std::string, std::vector<double>>& get_calculation_data();

        // This is a class method
        static int GetIntersectedCellIndex(std::vector<IntersectedCell> &cells, Grid::Cell grdcell);

    private:
        // intersecting well segment definition
        std::vector<Vector3d> entry_points_;
        std::vector<Vector3d> exit_points_;
        std::vector<double> segment_radius_;

        // per segment well index calculation data
        std::map<std::string, std::vector<double>> calculation_data_;

        // well index
        double well_index_;
    };
}
}

#endif //FIELDOPT_INTERSECTEDCELL_H
