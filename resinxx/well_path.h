//
// Created by bellout on 3/2/18.
//

#ifndef FIELDOPT_WELL_PATH_H
#define FIELDOPT_WELL_PATH_H

// -----------------------------------------------------------------
// STD
#include <vector>
#include <iostream>
#include <string>
#include <array>
#include <list>
#include <set>
#include <utility>
#include <cstddef>

// -----------------------------------------------------------------
// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>

// -----------------------------------------------------------------
// Qt
#include "QString"
#include "QDateTime"

// FieldOpt --------------------------------------------------------
#include "Reservoir/grid/cell.h"
#include "Reservoir/grid/grid.h"
#include "Reservoir/grid/eclgrid.h"
#include "resinxx/rixx_grid/rigrid.h"
#include "Settings/model.h"
#include "intersected_cell.h"
#include "FieldOpt-WellIndexCalculator/wellindexcalculator.h"
#include "FieldOpt-WellIndexCalculator/tests/wic_debug.hpp"
#include "Utilities/debug.hpp"
#include "Utilities/time.hpp"

// RESINSIGHT: FWK/VIZFWK/LIBCORE\LIBGEOMETRY ----------------------
#include "resinxx/rixx_core_geom/cvfVector2.h"
#include "resinxx/rixx_core_geom/cvfVector3.h"
#include "resinxx/rixx_core_geom/cvfPlane.h"
#include "resinxx/rixx_core_geom/cvfBoundingBox.h"
#include "resinxx/rixx_core_geom/cvfRay.h"
#include "resinxx/rixx_core_geom/cvfMath.h"
#include "resinxx/rixx_core_geom/cvfAssert.h"
#include "resinxx/rixx_core_geom/cvfBase.h"

// RESINSIGHT: FWK/APPFWK/COMMONCODE\VIZEXT\PROJDATAMOD ------------
#include "resinxx/rixx_app_fwk/cvfStructGrid.h"
#include "resinxx/rixx_app_fwk/cafHexGridIntersectionTools.h"

// RESINSIGHT: APPLICATIONCODE/RESERVOIRDATAMODEL ------------------
#include "resinxx/rixx_res_mod/cvfGeometryTools.h"
#include "resinxx/rixx_res_mod/RigCellGeometryTools.h"

// FieldOpt::RESINXX -----------------------------------------------
#include "resinxx/geometry_tools.h"
#include "resinxx/rixx_grid/rigrid.h"
#include "resinxx/rixx_grid/ricell.h"

namespace cvf {
class HexIntersectionInfo;
};

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

using std::vector;
using std::array;

// =================================================================
class WellPath // : public cvf::Object
{
 public:
  WellPath();

  vector <cvf::Vec3d> m_wellPathPoints;
  vector<double> m_measuredDepths;

//  void                        setDatumElevation(double value);
//  bool                        hasDatumElevation() const;
//  double                      datumElevation() const;
//  cvf::Vec3d                  interpolatedPointAlongWellPath(double measuredDepth) const;
//  double                      wellPathAzimuthAngle(const cvf::Vec3d& position) const;
//  void                        twoClosestPoints(const cvf::Vec3d& position, cvf::Vec3d* p1, cvf::Vec3d* p2) const;

//  std::pair<std::vector<cvf::Vec3d>, std::vector<double> >
//  clippedPointSubset(double startMD, double endMD) const;

//  std::vector<cvf::Vec3d>     wellPathPointsIncludingInterpolatedIntersectionPoint(double intersectionMeasuredDepth) const;

  static vector<cvf::HexIntersectionInfo>
  findRawHexCellIntersections(const RIGrid* grid,
                              const vector<cvf::Vec3d>& coords);

  static cvf::Vec3d calculateLengthInCell(const array<cvf::Vec3d, 8>& hexCorners,
                                          const cvf::Vec3d& startPoint,
                                          const cvf::Vec3d& endPoint);

  static cvf::Vec3d calculateLengthInCell(const RIGrid* grid,
                                          size_t cellIndex,
                                          const cvf::Vec3d& startPoint,
                                          const cvf::Vec3d& endPoint);

  static vector<size_t> findCloseCells(const RIGrid* grid,
                                       const cvf::BoundingBox& bb);

  static size_t findCellFromCoords(const RIGrid* grid,
                                   const cvf::Vec3d& coords,
                                   bool* foundCell);

 private:
  bool m_hasDatumElevation;
  double m_datumElevation;
};

//====================================================================

}
}

#endif //FIELDOPT_WELL_PATH_H
