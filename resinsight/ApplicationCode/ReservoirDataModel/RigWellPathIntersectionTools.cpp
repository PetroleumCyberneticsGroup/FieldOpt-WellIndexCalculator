/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2017 Statoil ASA
//
//  ResInsight is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ResInsight is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.
//
//  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
//  for more details.
//
/////////////////////////////////////////////////////////////////////////////////

// RESINSIGHT: APPLICATIONCODE/RESERVOIRDATAMODEL ------------------
#include "RigWellPathIntersectionTools.h"
#include "RigHexIntersectionTools.h"
#include "RigMainGrid.h"
#include "RigEclipseCaseData.h"
#include "RigCellGeometryTools.h"
#include "cvfGeometryTools.h"
#include "cvfMatrix3.h"
#include "RigEclipseWellLogExtractor.h"

//#include "RiaLogging.h"
//#include "RigWellPath.h"
//#include "RigWellLogExtractionTools.h"
//#include "RigEclipseWellLogExtractor.h"
//#include "RimEclipseCase.h"

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
//std::vector<WellPathCellIntersectionInfo>
//RigWellPathIntersectionTools::findCellIntersectionInfosAlongPath(
//    const RigEclipseCaseData* caseData,
//    const std::vector<cvf::Vec3d>& pathCoords,
//    const std::vector<double>& pathMds)
//{
//  std::vector<WellPathCellIntersectionInfo> intersectionInfos;
//  const RigMainGrid* grid = caseData->mainGrid();
//
//  if (pathCoords.size() < 2) return intersectionInfos;
//
//  cvf::ref<RigWellPath> dummyWellPath = new RigWellPath;
//  dummyWellPath->m_wellPathPoints = pathCoords;
//  dummyWellPath->m_measuredDepths = pathMds;
//
//  cvf::ref<RigEclipseWellLogExtractor> extractor =
//      new RigEclipseWellLogExtractor(
//          caseData,
//          dummyWellPath.p(),
//          caseData->ownerCase()->caseUserDescription().toStdString());
//
//  return extractor->cellIntersectionInfosAlongWellPath();
//}

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
std::vector<HexIntersectionInfo>
RigWellPathIntersectionTools::findRawHexCellIntersections(
    const RigMainGrid* grid,
    const std::vector<cvf::Vec3d>& coords) {

  std::vector<HexIntersectionInfo> intersections;

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
    std::vector<size_t> closeCells = findCloseCells(grid, bb);

    // Loop through cell neighborhood
    std::array<cvf::Vec3d, 8> hexCorners;
    for (size_t closeCell : closeCells) {

      // Get current cell
      const RigCell& cell = grid->globalCellArray()[closeCell];
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
      RigHexIntersectionTools::lineHexCellIntersection(
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
///
// -----------------------------------------------------------------
cvf::Vec3d
RigWellPathIntersectionTools::calculateLengthInCell(
    const std::array<cvf::Vec3d, 8>& hexCorners,
    const cvf::Vec3d& startPoint,
    const cvf::Vec3d& endPoint) {

  cvf::Vec3d vec = endPoint - startPoint;
  cvf::Vec3d iAxisDirection;
  cvf::Vec3d jAxisDirection;
  cvf::Vec3d kAxisDirection;

  RigCellGeometryTools::findCellLocalXYZ(hexCorners,
                                         iAxisDirection, jAxisDirection, kAxisDirection);

  cvf::Mat3d localCellCoordinateSystem(
      iAxisDirection.x(), jAxisDirection.x(), kAxisDirection.x(),
      iAxisDirection.y(), jAxisDirection.y(), kAxisDirection.y(),
      iAxisDirection.z(), jAxisDirection.z(), kAxisDirection.z());

  return vec.getTransformedVector(localCellCoordinateSystem.getInverted());
}

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
cvf::Vec3d
RigWellPathIntersectionTools::calculateLengthInCell(
    const RigMainGrid* grid,
    size_t cellIndex,
    const cvf::Vec3d& startPoint,
    const cvf::Vec3d& endPoint) {

  std::array<cvf::Vec3d, 8> hexCorners;
  grid->cellCornerVertices(cellIndex, hexCorners.data());

  return calculateLengthInCell(hexCorners, startPoint, endPoint);
}

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
std::vector<size_t>
RigWellPathIntersectionTools::findCloseCells(const RigMainGrid* grid,
                                             const cvf::BoundingBox& bb) {

  std::vector<size_t> closeCells;

  // ---------------------------------------------------------------
  const QDateTime tstart = QDateTime::currentDateTime();
  std::string str; str = "Find close cells.";
  print_dbg_msg_wic_ri(__func__, str, 0.0, 1);

  grid->findIntersectingCells(bb, &closeCells);

  print_dbg_msg_wic_ri(__func__, str, time_since_msecs(tstart), 2);
  return closeCells;
}

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
size_t
RigWellPathIntersectionTools::findCellFromCoords
    (const RigMainGrid* grid,
     const cvf::Vec3d& coords,
     bool* foundCell) {

  cvf::BoundingBox bb;
  bb.add(coords);

  // ---------------------------------------------------------------
  const QDateTime tstart = QDateTime::currentDateTime();
  std::string str; str = "Find cell from coords.";
  print_dbg_msg_wic_ri(__func__, str, 0.0, 1);

  std::vector<size_t> closeCells = findCloseCells(grid, bb);

  print_dbg_msg_wic_ri(__func__, str, time_since_msecs(tstart), 2);

  // ---------------------------------------------------------------
  str = "Looping through close cells.";
  print_dbg_msg_wic_ri(__func__, str, 0.0, 1);

  std::array<cvf::Vec3d, 8> hexCorners;
  for (size_t closeCell : closeCells) {

    const RigCell& cell = grid->globalCellArray()[closeCell];
    if (cell.isInvalid()) continue;

    grid->cellCornerVertices(closeCell, hexCorners.data());

    if (RigHexIntersectionTools::isPointInCell(coords, hexCorners.data())) {
      *foundCell = true;
      return closeCell;
    }
  }

  print_dbg_msg_wic_ri(__func__, str, time_since_msecs(tstart), 2);

  // ---------------------------------------------------------------
  *foundCell = false;
  return 0;
}
