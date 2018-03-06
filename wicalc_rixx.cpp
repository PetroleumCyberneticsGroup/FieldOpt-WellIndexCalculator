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
using std::vector;
using std::stringstream;

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

// -----------------------------------------------------------------
wicalc_rixx::~wicalc_rixx() {}

// -----------------------------------------------------------------
wicalc_rixx::wicalc_rixx(Settings::Model::Well well_settings,
                         Grid::Grid *grid) {

  well_settings_ = well_settings;
  grid_ = grid;
}


// -----------------------------------------------------------------
struct WellPathCellIntersectionInfo {

  size_t globCellIndex;
  cvf::Vec3d startPoint;
  cvf::Vec3d endPoint;
  double startMD;
  double endMD;
  cvf::Vec3d intersectionLengthsInCellCS;

//  cvf::StructGridInterface::FaceType intersectedCellFaceIn;
//  cvf::StructGridInterface::FaceType intersectedCellFaceOut;
};



// -----------------------------------------------------------------
enum CompletionType {
  WELL_PATH,
};

// -----------------------------------------------------------------
void wicalc_rixx::calculateWellPathIntersections(const WellPath *wellPath,
                                                 const RIGrid *grid,
                                                 vector<double> &isc_values) {

  vector<cvf::HexIntersectionInfo> intersections =
      WellPath::findRawHexCellIntersections(grid,
                                            wellPath->m_wellPathPoints);

  for (auto &intersection : intersections) {
    isc_values[intersection.m_hexIndex] = WELL_PATH;
  }
}

// -----------------------------------------------------------------
void wicalc_rixx::ComputeWellBlocks(
    map<string, vector<IntersectedCell>> &well_indices,
    vector<WellDefinition> &wells,
    int rank) {

  stringstream str;

  RIGrid* RIGrid_ = new RIGrid("  ");

  // Perform well block search for each well
  for (int iWell = 0; iWell < wells.size(); ++iWell) {

    // Intersected cells for well
    vector<IntersectedCell> intersected_cells;

    // Loop through well segments
    for (int iSegment = 0; iSegment < wells[iWell].radii.size(); ++iSegment) {

      vector<WellPathCellIntersectionInfo> intersectionInfos;

      wellPath_= new WellPath();

      // Add measuredepths
      wellPath_->m_measuredDepths.push_back(wells[iWell].heel_md[iSegment]);
      wellPath_->m_measuredDepths.push_back(wells[iWell].toe_md[iSegment]);

      // Add wellpathpoints
      cvf::Vec3d cvf_xyzHeel(wells[iWell].heels[iSegment][0],
                             wells[iWell].heels[iSegment][1],
                             -wells[iWell].heels[iSegment][2]);

      cvf::Vec3d cvf_xyzToe(wells[iWell].toes[iSegment][0],
                            wells[iWell].toes[iSegment][1],
                            -wells[iWell].toes[iSegment][2]);

      wellPath_->m_wellPathPoints.push_back(cvf_xyzHeel);
      wellPath_->m_wellPathPoints.push_back(cvf_xyzToe);

      calculateWellPathIntersections(wellPath_,
                                     RIGrid_,
                                     intersections_);

    }

  }

}

// -----------------------------------------------------------------



}
}