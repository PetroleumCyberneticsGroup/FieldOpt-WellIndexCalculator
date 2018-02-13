//
// Created by bellout on 2/13/18.
//

#include "wicalc.h"

wicalc::wicalc(){

};

wicalc::~wicalc(){};

void wicalc::calculateWellPathIntersections(const RigWellPath*   wellPath,
                                            const RigMainGrid*   grid,
                                            std::vector<double>& values) {

  std::vector<HexIntersectionInfo> intersections =
      RigWellPathIntersectionTools::findRawHexCellIntersections(
          grid,
          wellPath->m_wellPathPoints);

  for (auto& intersection : intersections)
  {
    values[intersection.m_hexIndex] = RiaDefines::WELL_PATH;
  }

}

//void wicalc::calculateWellPathIntersections(const RigWellPath*      wellPath,
//                                            const RimMainGrid*      grid,
//                                            std::vector<double>& values) {
//
//  std::vector <HexIntersectionInfo> intersections =
//      RigWellPathIntersectionTools::findRawHexCellIntersections(
//          grid,
//          wellPath->wellPathGeometry()->m_wellPathPoints);
//}