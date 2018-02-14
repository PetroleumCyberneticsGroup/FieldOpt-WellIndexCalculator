//
// Created by bellout on 2/13/18.
//

#ifndef FIELDOPT_WICALC_H
#define FIELDOPT_WICALC_H

#include "resinsight/ApplicationCode/ReservoirDataModel/RigMainGrid.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigGridBase.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigCell.h"

#include "resinsight/ApplicationCode/ReservoirDataModel/RigWellPath.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigWellPathIntersectionTools.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigHexIntersectionTools.h"
//#include "resinsight/ApplicationCode/ProjectDataModel/RimWellPath.h"
#include "resinsight/ApplicationCode/Application/RiaDefines.h"


#include "array"
#include "vector"
#include "list"

using std::cout;
using std::endl;
using std::vector;

//-----------------------------------------------------------------------------

//class WellPath {
//
//};

//class MainGrid {
//
//};


class wicalc {
  wicalc();
  ~wicalc();

//  void calculateWellPathIntersections(const RimWellPath*      wellPath,
//                                      const MainGrid*      grid,
//                                      std::vector<double>& values);

 private:
  void calculateWellPathIntersections(const RigWellPath*  wellPath,
                                      const RigMainGrid*  grid,
                                      vector<double>&     values);

  RigWellPath*  wellPath_;
  RigMainGrid*  grid_;
  vector<double>     values_;

};

#endif //FIELDOPT_WICALC_H
