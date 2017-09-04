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
#include "Reservoir/grid/cell.h"
#include "FieldOpt-WellIndexCalculator/intersected_cell.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

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
inline void print_wic_dbg(bool dbg_mode, bool append, int rank,
                          string dbg_loc, string dbg_msg) {

    stringstream dbg_file;
    dbg_file << "wic" << setw(3) << setfill('0') << rank << ".dbg";

    fstream fs;
    if (dbg_mode) {
        if (append) {
            fs.open(dbg_file.str(), std::fstream::out | std::fstream::app);
        } else {
            fs.open(dbg_file.str(), std::fstream::out | std::fstream::trunc);
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
    double xf, double yf, double zf,
    int rank) {

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
        dbg_mode, true, rank, "[TIME:" + get_time_stamp() +
                    "]\n[ComputeWellBlocks (WellIndexCalculator.cpp)] "
                            "Test numeric limit types: \n", dbg_msg.str());
};

/*!
 * \brief
 * This is the first function called by the WIC process
 */
inline void dbg_ComputeWellBlocks_bbox_i(
    bool dbg_mode,
    vector<Reservoir::WellIndexCalculation::WellDefinition> wells,
    int iWell, int iSegment,
    double xi, double yi, double zi,
    double xf, double yf, double zf,
    int rank) {

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
        "[[ WIC DEBUG ]] Timestamp: " + get_time_stamp() + "\n" +
        "@wellindexcalculator.cpp [ComputeWellBlocks ()]: "
        "Starting xyz well segment:\n";
    print_wic_dbg(dbg_mode, true, rank, start_str, dbg_msg.str());
    dbg_msg.str("");

    dbg_msg << "(xi yi zi): ["
            << setw(wdth) << xi << setw(wdth) << yi << setw(wdth) << zi << "]\n"
            << "(xf yf zf): ["
            << setw(wdth) << xf << setw(wdth) << yf << setw(wdth) << zf << "]\n";
    print_wic_dbg(
        dbg_mode, true, rank, "@wellindexcalculator.cpp [ComputeWellBlocks()]: "
            "Bounding box corresponding to xyz segment :\n", dbg_msg.str());
};

/*!
 * \brief Test bounding box after heuristic expansion
 */
inline void dbg_ComputeWellBlocks_bbox_f(
    bool dbg_mode,
    double xi, double yi, double zi,
    double xf, double yf, double zf,
    int rank) {

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
        dbg_mode, true, rank, "@wellindexcalculator.cpp [ComputeWellBlocks()]: "
            "B-box after heuristic increase:\n", dbg_msg.str());

};

// ---------------------------------------------------------------------
// eclgrid.cpp => GetBoundingBoxCellIndices

/*!
 * \brief Debug function vector<int> ECLGrid::GetBoundingBoxCellIndices
 * in eclgrid.cpp
 */
inline void dbg_GetBoundingBoxCellIndices(bool dbg_mode,
                                          vector<int> indices_list,
                                          int rank) {
    ostringstream dbg_msg;
    dbg_msg.precision(3);
    dbg_msg.setf(ios::fixed, ios::floatfield);
    dbg_msg.setf(ios::adjustfield, ios::right);
    dbg_msg << "\tindices_list: [\n";
    for (int ii = 0; ii < indices_list.size(); ++ii) {
        dbg_msg << setw(7) << indices_list[ii] << " ";
        if (remainder(ii+1, 16) == 0) { dbg_msg << "\n"; }
    }
    dbg_msg << "] => Total indices in B-box: "
            << indices_list.size() << "\n";

    print_wic_dbg(
        dbg_mode, true, rank, "@eclgrid.cpp [GetBoundingBoxCellIndices()]: "
            "Indices of current b-box:\n", dbg_msg.str());
};

// ---------------------------------------------------------------------
// wellindexcalculator.cpp

/*!
 * \brief Test overall algorithm for finding intersected cells
 */
inline void dbg_collect_intersected_cells_well_outside_box(
        bool dbg_mode, int error, int rank) {

    stringstream dbg_str;
    if (error == 1) {
        dbg_str << "\tWIC [RANK=" << rank
        << "]: Well or segment is outside of bounding box.";
    } else if(error == 1) {
        dbg_str << "\tWIC [RANK=" << rank
        << "]: Well or segment is outside of logical grid.";
    }
    cout << dbg_str.str() << endl;

    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp "
        "[collect_intersected_cells()]:\nError: \n", dbg_str.str());

};

// ---------------------------------------------------------------------
// wellindexcalculator.cpp
/*!
 * \brief
 */
inline void dbg_FindHeelToeEndPoints(bool dbg_mode,
                                     Reservoir::Grid::Cell &first_cell,
                                     Reservoir::Grid::Cell &last_cell,
                                     string result,
                                     int rank) {

    stringstream dbg_str;
    dbg_str << "\tWIC [RANK=" << rank;

    if (result=="failure") {

        dbg_str << "]: Failed to move well endpoints inside the reservoir.";

    } else if (result=="success") {

        dbg_str << "]: Managed to move well endpoints inside the reservoir, "
                << "or they were already inside."
                << "\tFirst cell global_idx = " << first_cell.global_index() << " -- "
                << "Last cell global_idx = " << last_cell.global_index();
    }

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: "
        "[ if-statement: Find the heel and toe cells ]\n", dbg_str.str());
};


// ---------------------------------------------------------------------
// wellindexcalculator.cpp
/*!
 * \brief
 */
inline void dbg_FindEndPointA(bool dbg_mode,
                             Vector3d start_pt,
                             double step,
                             Reservoir::Grid::Cell &cell,
                             vector<int> bb_cells,
                             int rank) {

    stringstream dbg_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);

    dbg_str << "\tWIC [RANK=" << rank << "]: "
            << "Endpoint (start_pt=" << start_pt.transpose() << ") is in a "
            << "cell that is ACTIVE (cell.is_active()=" << cell.is_active()
            << "). Breaking while loop...\n";

    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: "
        "[ if(grid_->GetCellEnvelopingPoint(cell, start_pt, bb_cells)) ] ",
        dbg_str.str());
};


// ---------------------------------------------------------------------
// wellindexcalculator.cpp
/*!
 * \brief
 */
inline void dbg_FindEndPointB(bool dbg_mode,
                             Vector3d old_start_pt,
                             Vector3d start_pt,
                             double step,
                             string dir,
                             int rank) {

    stringstream dbg_str, step_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);

    step_str.precision(6);
    step_str.setf(ios::fixed, ios::floatfield);
    step_str.setf(ios::adjustfield, ios::right);
    step_str << "(step=" << step << ")";

    dbg_str << "\tWIC [RANK=" << rank;
    if (dir=="inwards") {

        dbg_str << "]: Previous Endpoint (old_start_pt=" << old_start_pt.transpose()
                << ") was found to be in an INACTIVE cell. We now take a step="
                << step_str.str() << " INWARDS along the trajectory. Then test whether "
                << "the new endpoint (start_pt=" << start_pt.transpose() << ") "
                << "is inside an active cell.";

    } else if (dir=="outwards") {

        dbg_str << "]: Endpoint is now in an ACTIVE cell. We now take a "
                << step_str.str() << " OUTWARDS along the trajectory, and "
                << "test whether the new endpoint (start_pt=" << start_pt.transpose()
                << ") is now in an **inactive** cell. This is then the new endpoint "
                << "of the line (well).";
    }

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        // "[ start_pt = org_start_pt * (1 - step) + end_point * step ] -- ",
        dbg_str.str());
};

/*!
 * \brief
 */
inline void dbg_step(bool dbg_mode, Vector3d &start_pt,
                     Vector3d &end_pt, Vector3d &exit_pt,
                     int rank) {

    stringstream dbg_str, step_str, nom_ol, den_ul;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);

    step_str.precision(12);
    step_str.setf(ios::fixed, ios::floatfield);
    step_str.setf(ios::adjustfield, ios::right);

    nom_ol.precision(12);
    nom_ol.setf(ios::fixed, ios::floatfield);
    nom_ol.setf(ios::adjustfield, ios::right);

    den_ul.precision(12);
    den_ul.setf(ios::fixed, ios::floatfield);
    den_ul.setf(ios::adjustfield, ios::right);

    dbg_str << "\tWIC [RANK=" << rank;

    nom_ol << "NOM: (exit_pt - start_pt).norm()="
           << (exit_pt - start_pt).norm() << " -- ISNAN: "
           << isnan((exit_pt - start_pt).norm()) << "";

    if (isnan((exit_pt - start_pt).norm())) {
        nom_ol << "\n(exit_pt=" << exit_pt.transpose() << "); "
               << "\n(start_pt=" << start_pt.transpose() << ");\n";
    }

    den_ul << "DEN: (end_pt - start_pt).norm()="
           << (end_pt - start_pt).norm() << " -- ISNAN: "
           << isnan((end_pt - start_pt).norm()) << "";

    if (isnan((end_pt - start_pt).norm())) {
        nom_ol << "\n(end_pt=" << end_pt.transpose() << "); "
               << "\n(start_pt=" << start_pt.transpose() << ");\n";
    }

    double step_loc = (exit_pt - start_pt).norm() / (end_pt - start_pt).norm();
    step_str << "(step=" << step_loc << ")";

    dbg_str << step_str.str() << nom_ol.str() << den_ul.str();

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        dbg_str.str());

};

/*!
 * \brief
 */
inline void dbg_TraversingCellsA(bool dbg_mode,
                                 Reservoir::Grid::Cell &new_cell,
                                 Reservoir::Grid::Cell &prev_cell,
                                 Vector3d &old_entry_pt, Vector3d &entry_pt,
                                 Vector3d &start_pt, Vector3d &end_pt,
                                 double step,
                                 double epsilon,
                                 int steps,
                                 string activity,
                                 int rank) {


    stringstream dbg_str, step_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);

    step_str.precision(6);
    step_str.setf(ios::fixed, ios::floatfield);
    step_str.setf(ios::adjustfield, ios::right);
    step_str << "(step=" << step << " [eps=" << epsilon << ",  " << steps << "])";

    dbg_str << "\tWIC [RANK=" << rank;
    if (activity=="step-into") {

        dbg_str << "]: We start from (start_pt=" << start_pt.transpose()
                << ") and move to (entry_pt=" << entry_pt.transpose() << ") using "
                << step_str.str();

    } else {

        if (activity=="check-ok") {

            dbg_str << "]: Found new_cell=" << new_cell.global_index()
                    << " enveloping (entry_pt=" << entry_pt.transpose()
                    << ") While-loop now checks idx_new_cell is NOT_EQUAL to idx_prev_cell. "
                    << "If so, make another step into this cell, and try again...";

        } else if (activity=="check-failed") {

            dbg_str << "]: While at (entry_pt=" << entry_pt.transpose() << "), we checked "
                    << "if we can find the cell enveloping this point. But this failed. "
                    << "We discontinue this do-loop." ;


        }

        dbg_str << "new_cell.global_index()=" << new_cell.global_index() << " -- "
        << "new_cell.is_active()=" << new_cell.is_active() << " -- "
        << "prev_cell.global_index()=" << prev_cell.global_index();

    }

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        dbg_str.str());

};

/*!
 * \brief
 */
inline void dbg_TraversingCellsB(bool dbg_mode,
                                 Reservoir::Grid::Cell &new_cell,
                                 Reservoir::Grid::Cell &prev_cell,
                                 Reservoir::Grid::Cell &last_cell,
                                 double step,
                                 vector<Reservoir::WellIndexCalculation::IntersectedCell> &isc_cells,
                                 int isc_cell_idx,
                                 Vector3d &entry_pt, Vector3d &end_pt,
                                 Vector3d &old_exit_pt, Vector3d &exit_pt,
                                 int rank) {

    stringstream dbg_str, step_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);

    step_str.precision(6);
    step_str.setf(ios::fixed, ios::floatfield);
    step_str.setf(ios::adjustfield, ios::right);
    step_str << "(step=" << step << ")";

    dbg_str << "\tWIC [RANK=" << rank;
    dbg_str << "]: At new cell. Finding (exit_pt=" << exit_pt.transpose() << ") -- "
            << "(new_cell.global_index() != prev_cell.global_index()) => "
            << ( new_cell.global_index() != prev_cell.global_index() ) << " -- "
            << "(new_cell.global_index() != last_cell.global_index()) => "
            << ( new_cell.global_index() != last_cell.global_index() )
            << "\n\tNEW CELL";

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        dbg_str.str());

};


/*!
 * \brief
 */
inline void dbg_TraversingCellsC(bool dbg_mode,
                                 Reservoir::Grid::Cell &new_cell,
                                 Reservoir::Grid::Cell &last_cell,
                                 Vector3d &entry_pt, Vector3d &end_pt,
                                 double step,
                                 vector<Reservoir::WellIndexCalculation::IntersectedCell> &isc_cells,
                                 int isc_cell_idx,
                                 string activity,
                                 int rank) {

    stringstream dbg_str, step_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);
    dbg_str << "\tWIC [RANK=" << rank;

    step_str.precision(6);
    step_str.setf(ios::fixed, ios::floatfield);
    step_str.setf(ios::adjustfield, ios::right);
    step_str << "(step=" << step << ")";

    if (activity=="step-into") {

        dbg_str << "]: Stepped beyond orig_well_length " << step_str.str()
                << " have (end_pt=" << end_pt.transpose() << "), or at last_cell: "
                << "(new_cell.global_index() == last_cell.global_index()) => "
                << ( new_cell.global_index() == last_cell.global_index() );

    } else if (activity=="check") {

        dbg_str << "]: ";

        cout << "WIC [RANK=" << rank << "]: Expected last cell does not match "
            "found last cell. (DEBUG: NOT!) Returning empty list." << endl;
        // \TODO Debug:

    }

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        dbg_str.str());

};


/*!
 * \brief
 */
inline void dbg_recover_from_cycle(bool dbg_mode,
                                   string activity,
                                   int rank) {

    stringstream dbg_str;
    dbg_str.precision(3);
    dbg_str.setf(ios::fixed, ios::floatfield);
    dbg_str.setf(ios::adjustfield, ios::right);
    dbg_str << "\tWIC [RANK=" << rank;

    if (activity=="step-into") {


    } else if (activity=="check") {

    }

    dbg_str << endl;
    print_wic_dbg(dbg_mode, true, rank, "@wellindexcalculator.cpp: ",
        dbg_str.str());

};

}
#endif //FIELDOPT_WIC_DEBUG_H
