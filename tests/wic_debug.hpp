/******************************************************************************
   Copyright (C) 2017

   Created by bellout on 6/27/17.

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

#ifndef FIELDOPT_WIC_DEBUG_H
#define FIELDOPT_WIC_DEBUG_H

#include "FieldOpt-WellIndexCalculator/wellindexcalculator.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

namespace WICDebug {

// WIC Debugging functions

// ---------------------------------------------------------------------
// Debug print functions

/*!
 * \brief Prints debug messages
 * @param debug_msg
 * @param dbg_loc
 */
inline void print_wic_dbg(bool dbg_mode, bool append,
                          string dbg_loc, string dbg_msg) {
    string dbg_file = "wic.dbg";
    fstream fs;
    if (dbg_mode) {
        if (append) {
            fs.open(dbg_file, std::fstream::out | std::fstream::app);
        } else {
            fs.open(dbg_file, std::fstream::out | std::fstream::trunc);
        }
        fs.write(dbg_loc.c_str(), dbg_loc.size());
        fs.write(dbg_msg.c_str(), dbg_msg.size());
        fs.close();
    }
};

/*!
 * \brief Format of stream
 */
inline ostringstream get_dbg_msg(ostringstream * dbg_msg) {
//    ostringstream dbg_msg;
    dbg_msg->precision(3);
    dbg_msg->setf(ios::fixed, ios::floatfield);
    dbg_msg->setf(ios::adjustfield, ios::right);
//    return dbg_msg;
};

/*!
 * \brief Get time stamp
 */
inline string get_time_stamp()
{
    time_t ltime;
    struct tm *Tm;
    char ts_char [50];
    string ts_str;

    ltime = time(NULL); /* get current cal time */
    sprintf(ts_char, "%s", asctime( localtime(&ltime) ) );
    ts_str = string(ts_char);

    // Remove newline character
    ts_str.erase(remove(ts_str.begin(), ts_str.end(), '\n'), ts_str.end());

    return ts_str;
}

// ---------------------------------------------------------------------
// wellindexcalculator.cpp

/*!
 * \brief Test numeric limit types used to find bounding box
 *
 * Use:

    // Debug --------------------------------
    WICDebug::dbg_ComputeWellBlocks_num_lims(
        dbg_mode, xi, yi, zi, xf, yf, zf);
 *
 */
inline void dbg_ComputeWellBlocks_num_lims(
    bool dbg_mode,
    double xi, double yi, double zi,
    double xf, double yf, double zf) {
    int wdth = 18;

    ostringstream dbg_msg;
    dbg_msg.precision(3);
    dbg_msg.setf(ios::fixed, ios::floatfield);
    dbg_msg.setf(ios::adjustfield, ios::right);
    dbg_msg << "(xi yi zi): ["
            << setw(wdth) << xi << setw(wdth) << yi << setw(wdth) << zi << "]\n"
            << "(xf yf zf): ["
            << setw(wdth) << xf << setw(wdth) << yf << setw(wdth) << zf << "]\n";
    print_wic_dbg(
        dbg_mode, true, "[TIME:" + get_time_stamp() +
                    "]\n[ComputeWellBlocks (WellIndexCalculator.cpp)] "
                            "Test numeric limit types: \n", dbg_msg.str());
};

/*!
 * \brief
 * This is the first function called by the WIC process
 *
 * Use:

    // Debug ------------------------------
    WICDebug::dbg_ComputeWellBlocks_bbox_i(
        dbg_mode, wells, iWell, iSegment,
        xi, yi, zi, xf, yf, zf);
 */
inline void dbg_ComputeWellBlocks_bbox_i(
    bool dbg_mode,
    vector<Reservoir::WellIndexCalculation::WellDefinition> wells,
    int iWell, int iSegment,
    double xi, double yi, double zi,
    double xf, double yf, double zf) {

    int wdth = 11;
    time_t rawtime;
    ostringstream dbg_msg;
    dbg_msg.precision(3);
    dbg_msg.setf(ios::fixed, ios::floatfield);
    dbg_msg.setf(ios::adjustfield, ios::right);

    dbg_msg << "Heel (xyz): ["
            << setw(wdth) << wells[iWell].heels[iSegment].x()
            << setw(wdth) << wells[iWell].heels[iSegment].y()
            << setw(wdth) << wells[iWell].heels[iSegment].z() << "]\n"
            << "Toe (xyz):  ["
            << setw(wdth) << wells[iWell].toes[iSegment].x()
            << setw(wdth) << wells[iWell].toes[iSegment].y()
            << setw(wdth) << wells[iWell].toes[iSegment].z() << "]\n";


    string start_str = "\n" + string(64,'=') + "\n" +
        "WIC DEBUG Version: " + get_time_stamp() + "\n" +
        "[ComputeWellBlocks (wellindexcalculator.cpp)] "
        "Starting xyz well segment:\n";
    print_wic_dbg(dbg_mode, true, start_str, dbg_msg.str());
    dbg_msg.str("");

    dbg_msg << "(xi yi zi): ["
            << setw(wdth) << xi << setw(wdth) << yi << setw(wdth) << zi << "]\n"
            << "(xf yf zf): ["
            << setw(wdth) << xf << setw(wdth) << yf << setw(wdth) << zf << "]\n";
    print_wic_dbg(
        dbg_mode, true, "[ComputeWellBlocks (wellindexcalculator.cpp)] "
            "Bounding box corresponding to xyz segment :\n", dbg_msg.str());
};

/*!
 * \brief Test bounding box after heuristic expansion
 *
 * Use:

    // Debug ------------------------------
    WICDebug::dbg_ComputeWellBlocks_bbox_f(
        dbg_mode, xi, yi, zi, xf, yf, zf);
 */
inline void dbg_ComputeWellBlocks_bbox_f(
    bool dbg_mode,
    double xi, double yi, double zi,
    double xf, double yf, double zf) {
    int wdth = 11;

    ostringstream dbg_msg;
    dbg_msg.precision(3);
    dbg_msg.setf(ios::fixed, ios::floatfield);
    dbg_msg.setf(ios::adjustfield, ios::right);
    dbg_msg << "(xi yi zi): ["
            << setw(wdth) << xi << setw(wdth) << yi << setw(wdth) << zi << "]\n"
            << "(xf yf zf): ["
            << setw(wdth) << xf << setw(wdth) << yf << setw(wdth) << zf << "]\n";
    print_wic_dbg(
        dbg_mode, true, "[ComputeWellBlocks (WellIndexCalculator.cpp)] "
            "B-box after heuristic increase:\n", dbg_msg.str());

};

// ---------------------------------------------------------------------
// eclgrid.cpp => GetBoundingBoxCellIndices

/*!
 * \brief Debug function vector<int> ECLGrid::GetBoundingBoxCellIndices
 * in eclgrid.cpp
 *
 * Use:

    // Debug -------------------------------
    WICDebug::dbg_GetBoundingBoxCellIndices(
        dbg_mode, bb_cells);
 */
inline void dbg_GetBoundingBoxCellIndices(bool dbg_mode,
                                          vector<int> indices_list) {
    ostringstream dbg_msg;
    dbg_msg.precision(3);
    dbg_msg.setf(ios::fixed, ios::floatfield);
    dbg_msg.setf(ios::adjustfield, ios::right);
    dbg_msg << "indices_list: [\n";
    for (int ii = 0; ii < indices_list.size(); ++ii) {
        dbg_msg << setw(7) << indices_list[ii] << " ";
        if (remainder(ii+1, 8) == 0) { dbg_msg << "\n"; }
    }
    dbg_msg << "]\n";
    dbg_msg << "Total indices in B-box: "
            << indices_list.size() << "\n";

    print_wic_dbg(
        dbg_mode, true, "[GetBoundingBoxCellIndices (eclgrid.cpp)] "
            "Indices of current b-box:\n", dbg_msg.str());
};

// ---------------------------------------------------------------------
// wellindexcalculator.cpp

/*!
 * \brief Test overall algorithm for finding intersected cells
 *
 * Use:

 *
 */
inline void dbg_collect_intersected_cells_well_outside_box(
        bool dbg_mode, string dbg_str) {
    print_wic_dbg(dbg_mode, true, "[collect_intersected_cells "
        "(WellIndexCalculator.cpp)] Error: \n", dbg_str);
};

// ---------------------------------------------------------------------
// wellindexcalculator.cpp
/*!
 * \brief
 *
 *
 * Use:

    // Debug -------------------------------
    WICDebug::dbg_FindHeelToeEndPoints(bool dbg_mode, string dbg_str);
 */
inline void dbg_FindHeelToeEndPoints(bool dbg_mode, string dbg_str) {
    print_wic_dbg(dbg_mode, true, "[if-statement: Find the heel and "
            "toe cells (WellIndexCalculator.cpp)] Error: \n", dbg_str);
};


// ---------------------------------------------------------------------
// wellindexcalculator.cpp
/*!
 * \brief
 *
 *
 * Use:

    // Debug -------------------------------
    WICDebug::dbg_FindHeelToeEndPoints(bool dbg_mode, string dbg_str);
 */

// bool WellIndexCalculator::findEndpoint(const vector<int> &bb_cells,
//                                        Vector3d &start_pt,
//                                        Vector3d end_point,
//                                        Grid::Cell &cell)


}
#endif //FIELDOPT_WIC_DEBUG_H
