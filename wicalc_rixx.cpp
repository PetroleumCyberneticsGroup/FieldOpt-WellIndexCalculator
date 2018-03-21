//
// Created by bellout on 3/2/18.
//

// -----------------------------------------------------------------
#include "wicalc_rixx.h"

// -----------------------------------------------------------------
using std::cout;
using std::endl;
using std::list;
using std::pair;
using std::string;
using std::fill;
using std::vector;
using std::stringstream;

#include <memory>

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

// -----------------------------------------------------------------
wicalc_rixx::~wicalc_rixx() {
  //delete RIReaderECL_;
  //delete RICaseData_;
}

// -----------------------------------------------------------------
wicalc_rixx::wicalc_rixx(Settings::Model::Well well_settings,
                         Grid::Grid *grid) {

  well_settings_ = well_settings;
  cl_ = well_settings_.verb_vector_[3]; // current dbg.msg.level
  grid_ = grid;

  RIReaderECL_ = new RIReaderECL();
  RICaseData_ = new RICaseData(grid_->GetFilePath());
  RIReaderECL_->open(grid_->GetFilePathQString(), RICaseData_);
  RICaseData_->computeActiveCellBoundingBoxes();
  RICaseData_->mainGrid()->computeCachedData();

  RIGrid_ = RICaseData_->mainGrid();

  grid_count_ = RICaseData_->mainGrid()->gridCount();
  cell_count_ = RICaseData_->mainGrid()->cellCount();
  gcellarray_sz_ = RICaseData_->mainGrid()->globalCellArray().size();


  intersections_.resize(gcellarray_sz_);
  fill(intersections_.begin(), intersections_.end(), HUGE_VAL);

  // ---------------------------------------------------------------
  // Dbg
  QDateTime tstart = QDateTime::currentDateTime();
  std::stringstream str; str << "Find cell from coords.";
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, true, cl_, 2);

  str.str(""); str << "grid_->gridCount(): " << grid_count_
                   << " -- grid_->cellCount(): " << cell_count_
                   << " -- grid_->globalCellArray().size(): "
                   << gcellarray_sz_;
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, true, cl_, 2);

}

// -----------------------------------------------------------------
enum CompletionType {
  WELL_PATH
};

// -----------------------------------------------------------------
void wicalc_rixx::calculateWellPathIntersections(const WellPath& wellPath,
                                                 const RIGrid *grid,
                                                 vector<double> &isc_values) {

  vector<cvf::HexIntersectionInfo> intersections =
      WellPath::findRawHexCellIntersections(grid,
                                            wellPath.m_wellPathPoints);

  std::stringstream str;
  if (intersections.size() > 0) {
    str << "Set values to WELL_PATH. -- # of intersections = " << intersections.size()
        << " -- intersection[0].m_hexIndex = " << intersections[0].m_hexIndex;
    print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);
  }

  for (auto &intersection : intersections) {

    str.str("");
    str << "intersection.m_hexIndex = "
        << intersection.m_hexIndex
        << " -- values[intersection.m_hexIndex] = "
        << isc_values[intersection.m_hexIndex];
    // print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);
    isc_values[intersection.m_hexIndex] = CompletionType::WELL_PATH;
  }
}

// -----------------------------------------------------------------
void
wicalc_rixx::collectIntersectedCells(vector<IntersectedCell> &isc_cells,
                                     vector<WellPathCellIntersectionInfo> isc_info,
                                     WellDefinition well,
                                     WellPath& wellPath,
                                     int rank) {

  vector<RICompData> completionData;
  Reservoir::Grid::Cell gcell;

  for (auto& cell : isc_info) {

    // -------------------------------------------------------------
    // Get IJK idx
    size_t i, j, k;
    RIGrid_->ijkFromCellIndex(cell.globCellIndex, &i, &j, &k);

    // Check if cell is active, if not, skip
    bool cellIsActive = activeCellInfo_->isActive(cell.globCellIndex);
    if (!cellIsActive) {
      // cout << "Cell is not active" << endl;
      continue;
    }

    // -------------------------------------------------------------
    // Make RI Completion object
    RICompData completion(QString::fromStdString(well.wellname),
                          IJKCellIndex(i, j, k));

    // -------------------------------------------------------------
    // Make FO Cell object + fill values for trans.computation
    gcell = grid_->GetCell(cell.globCellIndex);
    IntersectedCell icell(gcell);

    // -------------------------------------------------------------
    // Calculate direction
    CellDir direction =
        wellPath.calculateDirectionInCell(RICaseData_,
                                          cell,
                                          icell);

    // -------------------------------------------------------------
    // Calculate transmissibility
    double transmissibility =
        wellPath.calculateTransmissibility(RICaseData_,
                                           cell.intersectionLengthsInCellCS,
                                           well.skins[0],
                                           well.radii[0],
                                           cell.globCellIndex,
                                           false, icell);

    // -------------------------------------------------------------
    // Deleted un susbsequent versions
    // Store calculated values in RI completion object (for
    // completeness sake, not necessary for our transfer)
    // completion.setTransAndWPImultBackgroundDataFromPerforation(transmissibility,
    //                                                            well.skins[0],
    //                                                            well.radii[0],
    //                                                            direction);
    // completion.addMetadata("Perforation",
    //                        QString("StartMD: %1 - EndMD: %2")
    //                            .arg(wellPath_->m_measuredDepths[0])
    //                            .arg(wellPath_->m_measuredDepths[1])
    //                            + QString(" : ") + QString::number(transmissibility));
    // completionData.push_back(completion);

    // -------------------------------------------------------------
    // Convert start + exit points + calculated lengths in cell for
    // transfer to FO object
    Vector3d start_pt(cell.startPoint.x(),
                      cell.startPoint.y(),
                      cell.startPoint.z());

    Vector3d exit_pt(cell.endPoint.x(),
                     cell.endPoint.y(),
                     cell.endPoint.z());

    Vector3d isc_lengths(cell.intersectionLengthsInCellCS.x(),
                         cell.intersectionLengthsInCellCS.y(),
                         cell.intersectionLengthsInCellCS.z());

    // -------------------------------------------------------------
    // Transfer segment data, incl. wccf, to FO intersected cell
    icell.add_new_segment(start_pt, exit_pt, well.radii[0], well.skins[0]);
    icell.set_cell_well_index_matrix(transmissibility);

    // Add to vector of intersected cells
    isc_cells.push_back(icell);

  }

}


// -----------------------------------------------------------------
void
wicalc_rixx::ComputeWellBlocks(map<string, vector<IntersectedCell>> &well_indices,
                               vector<WellDefinition> &wells,
                               int rank) {

  stringstream str;


  // Perform well block search for each well
  for (int iWell = 0; iWell < wells.size(); ++iWell) {

    // Intersected cells for well
    vector<IntersectedCell> intersected_cells;

    // Loop through well segments
    for (int iSegment = 0; iSegment < wells[iWell].radii.size(); ++iSegment) {

      cvf::ref<WellPath> wellPath = new WellPath();
      vector<WellPathCellIntersectionInfo> intersectionInfos;

      // Load measuredepths onto wellPath (= current segment)
      wellPath->m_measuredDepths.push_back(wells[iWell].heel_md[iSegment]);
      wellPath->m_measuredDepths.push_back(wells[iWell].toe_md[iSegment]);

      // Load wellpathpoints onto wellPath (= current segment)
      cvf::Vec3d cvf_xyzHeel(wells[iWell].heels[iSegment][0],
                             wells[iWell].heels[iSegment][1],
                             -wells[iWell].heels[iSegment][2]);

      cvf::Vec3d cvf_xyzToe(wells[iWell].toes[iSegment][0],
                            wells[iWell].toes[iSegment][1],
                            -wells[iWell].toes[iSegment][2]);

      wellPath->m_wellPathPoints.push_back(cvf_xyzHeel);
      wellPath->m_wellPathPoints.push_back(cvf_xyzToe);

      // Calculate cells intersected by well path
      calculateWellPathIntersections(*wellPath,
                                     RICaseData_->mainGrid(),
                                     intersections_);

      // Use intersection data to find intersected cell data
      cvf::ref<RIExtractor>
          extractor = new RIECLExtractor(RICaseData_, *wellPath);

      activeCellInfo_ = RICaseData_->activeCellInfo(MATRIX_MODEL);
      vector<WellPathCellIntersectionInfo>
          intersectedCellInfo = extractor->cellIntersectionInfosAlongWellPath();

      collectIntersectedCells(intersected_cells,
                              intersectedCellInfo,
                              wells[iWell],
                              *wellPath,
                              rank);
    }

    // Assign intersected cells to well
    well_indices[wells[iWell].wellname] = intersected_cells;

  } // # of wells

}

// -----------------------------------------------------------------



}
}