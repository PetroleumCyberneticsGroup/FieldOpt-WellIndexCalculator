//
// Created by bellout on 2/13/18.
//


#include <iostream>
#include "wicalc.h"

wicalc::wicalc(){

  wellPath_ = new RigWellPath();
  grid_ = new RigMainGrid();
  values_ = vector<double>();

  calculateWellPathIntersections(wellPath_, grid_, values_);

  for (unsigned int i = 0; i < values_.size(); i++) {
    cout << values_[i] << " " << endl;
    if (i % 10 == 1) cout << endl;
  }

}

wicalc::~wicalc(){}

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