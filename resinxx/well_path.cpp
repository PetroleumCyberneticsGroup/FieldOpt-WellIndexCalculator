//
// Created by bellout on 3/2/18.
//

// -----------------------------------------------------------------
#include "well_path.h"

namespace cvf {
class HexIntersectionInfo;
};

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

// -----------------------------------------------------------------
WellPath::WellPath()
    : m_hasDatumElevation(false),
      m_datumElevation(0.0) {}

// -----------------------------------------------------------------
vector<cvf::HexIntersectionInfo>
WellPath::findRawHexCellIntersections(const RIGrid* grid,
                                      const vector<cvf::Vec3d>& coords) {

  vector<cvf::HexIntersectionInfo> intersections;

  // ---------------------------------------------------------------
  const QDateTime tstart = QDateTime::currentDateTime();
  std::stringstream str0, str1;
  str0 << "Find raw hexcell intersections. grid->nodes().size() = "
       << grid->nodes().size() << " coords.size() = " << coords.size();
  print_dbg_msg_wic_ri(__func__, str0.str(), 0.0, 1);

  for (size_t i = 0; i < coords.size() - 1; ++i) {

    // Dbg: coord i
    str1.str("");
    str1 << "coord[i=" << i << "]=( "
         << std::setw(10) << std::setprecision(3) << std::fixed
         << coords[i].x() << ", "
         << coords[i].y() << ", "
         << coords[i].z() << " )";
    print_dbg_msg_wic_ri(__func__, str1.str(), 0.0, 0);

    // Dbg: coord i+1
    str1.str("");
    str1 << "coord[i=" << i + 1 << "]=( "
         << std::setw(10) << std::setprecision(3) << std::fixed
         << coords[i+1].x() << ", "
         << coords[i+1].y() << ", "
         << coords[i+1].z() << " )";
    print_dbg_msg_wic_ri(__func__, str1.str(), 0.0, 0);

    // Add coords to bbox
    cvf::BoundingBox bb;
    bb.add(coords[i]);
    bb.add(coords[i + 1]);

    // Find cells close to bbox
    vector<size_t> closeCells = findCloseCells(grid, bb);

    // Loop through cell neighborhood
    array<cvf::Vec3d, 8> hexCorners;
    for (size_t closeCell : closeCells) {

      // Get current cell
      const RICell& cell = grid->globalCellArray()[closeCell];
      if (cell.isInvalid()) {
        print_dbg_msg_wic_ri(__func__, "Cell is invalid", 0.0, 0);
        continue;
      }

      // Get corner vertices of current cell
      grid->cellCornerVertices(closeCell, hexCorners.data());

      // Dbg
      str1.str("");
      for (int i=0; i < hexCorners.size(); i++) {

        if (i < 1 && intersections.size() < 4) {
          str1 << std::setw(10) << std::setprecision(3) << std::fixed
               << "hexCorner[i=" << i << "]=( "
               << hexCorners[i].x() << ", "
               << hexCorners[i].y() << ", "
               << hexCorners[i].z() << " ) -- intersection.sz = "
               << intersections.size();
//          print_dbg_msg_wic_ri(__func__, str1.str(), 0.0, 0);
        }
      }

      //
      cvf::RigHexIntersectionTools::lineHexCellIntersection(
          coords[i], coords[i + 1], hexCorners.data(),
          closeCell, &intersections);
    } // End: for (size_t closeCell : closeCells)
  }

  str1.str("");
  str1 << "# of interections found = " << intersections.size();
  print_dbg_msg_wic_ri(__func__, str1.str(), 0.0, 0);

  print_dbg_msg_wic_ri(__func__, str0.str(), time_since_msecs(tstart), 2);

  return intersections;
}

// -----------------------------------------------------------------
cvf::Vec3d WellPath::calculateLengthInCell(
    const array<cvf::Vec3d, 8>& hexCorners,
    const cvf::Vec3d& startPoint,
    const cvf::Vec3d& endPoint) {

  cvf::Vec3d vec = endPoint - startPoint;
  cvf::Vec3d iAxisDirection;
  cvf::Vec3d jAxisDirection;
  cvf::Vec3d kAxisDirection;

  RigCellGeometryTools::findCellLocalXYZ(hexCorners,
                                         iAxisDirection,
                                         jAxisDirection,
                                         kAxisDirection);

  cvf::Mat3d localCellCoordinateSystem(
      iAxisDirection.x(), jAxisDirection.x(), kAxisDirection.x(),
      iAxisDirection.y(), jAxisDirection.y(), kAxisDirection.y(),
      iAxisDirection.z(), jAxisDirection.z(), kAxisDirection.z());

  return vec.getTransformedVector(localCellCoordinateSystem.getInverted());
}

// -----------------------------------------------------------------
cvf::Vec3d WellPath::calculateLengthInCell(const RIGrid* grid,
                                           size_t cellIndex,
                                           const cvf::Vec3d& startPoint,
                                           const cvf::Vec3d& endPoint) {

  std::array<cvf::Vec3d, 8> hexCorners;
  grid->cellCornerVertices(cellIndex, hexCorners.data());

  return calculateLengthInCell(hexCorners, startPoint, endPoint);
}

// -----------------------------------------------------------------
vector<size_t> WellPath::findCloseCells(const RIGrid* grid,
                                        const cvf::BoundingBox& bb) {

  vector<size_t> closeCells;
  grid->findIntersectingCells(bb, &closeCells);
  return closeCells;
}

// -----------------------------------------------------------------
size_t WellPath::findCellFromCoords(const RIGrid* grid,
                                    const cvf::Vec3d& coords,
                                    bool* foundCell) {

  cvf::BoundingBox bb;
  bb.add(coords);
  vector<size_t> closeCells = findCloseCells(grid, bb);

  array<cvf::Vec3d, 8> hexCorners;
  for (size_t closeCell : closeCells) {

    const RICell& cell = grid->globalCellArray()[closeCell];
    if (cell.isInvalid()) continue;

    grid->cellCornerVertices(closeCell, hexCorners.data());

    if (cvf::RigHexIntersectionTools::isPointInCell(coords, hexCorners.data())) {
      *foundCell = true;
      return closeCell;
    }
  }

  *foundCell = false;
  return 0;
}


}
}
