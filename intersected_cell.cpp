/******************************************************************************
   Copyright (C) 2015-2016 Einar J.M. Baumann <einar.baumann@gmail.com>
???? Copyright (C) 2015-2016 Hilmar M. Magnusson <hilmarmag@gmail.com>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#include "intersected_cell.h"

namespace Reservoir {
namespace WellIndexCalculation {

std::vector<Vector3d> IntersectedCell::points() const {
    return std::vector<Vector3d>({entry_point_, exit_point_});
}

Vector3d IntersectedCell::xvec() const {
    return corners()[5] - corners()[4];
}

Vector3d IntersectedCell::yvec() const {
    return corners()[6] - corners()[4];
}

Vector3d IntersectedCell::zvec() const {
    return corners()[0] - corners()[4];
}

double IntersectedCell::dx() const {
    return xvec().norm();
}

double IntersectedCell::dy() const {
    return yvec().norm();
}

double IntersectedCell::dz() const {
    return zvec().norm();
}

const Vector3d &IntersectedCell::entry_point() const {
    return entry_point_;
}

void IntersectedCell::set_entry_point(const Vector3d &entry_point) {
    entry_point_ = entry_point;
}

const Vector3d &IntersectedCell::exit_point() const {
    return exit_point_;
}

void IntersectedCell::set_exit_point(const Vector3d &exit_point) {
    exit_point_ = exit_point;
}

double IntersectedCell::well_index() const {
    return well_index_;
}

void IntersectedCell::set_well_index(double well_index) {
    well_index_ = well_index;
}
}
}
