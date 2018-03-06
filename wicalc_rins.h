//
// Created by bellout on 2/13/18.
//

#ifndef FIELDOPT_WICALC_RI_H
#define FIELDOPT_WICALC_RI_H

// -----------------------------------------------------------------
// STD
#include <vector>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <list>

// -----------------------------------------------------------------
// EIGEN
#include <Eigen/Core>

// -----------------------------------------------------------------
// Qt
#include "QString"

// -----------------------------------------------------------------
// FieldOpt
#include "Reservoir/grid/grid.h"
#include "Settings/model.h"

// RESINSIGHT: FWK/VIZFWK/LIBCORE ----------------------------------
#include "resinsight/Fwk/VizFwk/LibCore/cvfBase.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfObject.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfMath.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfVector3.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfCollection.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfAssert.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfArray.h"
//#include "resinsight/Fwk/VizFwk/LibCore/cvfVector4.h"
//#include "resinsight/Fwk/VizFwk/LibCore/cvfColor4.h"

// RESINSIGHT: FWK/APPFWK/COMMONCODE -------------------------------
#include "resinsight/Fwk/AppFwk/CommonCode/cvfStructGrid.h"

// RESINSIGHT: APPLICATIONCODE/APPLICATION -------------------------
#include "resinsight/ApplicationCode/Application/RiaDefines.h"
#include "resinsight/ApplicationCode/Application/RiaPorosityModel.h"

// RESINSIGHT: APPLICATIONCODE/FILEINTERFACE -----------------------
#include "resinsight/ApplicationCode/FileInterface/RifReaderEclipseOutput.h"
#include "resinsight/ApplicationCode/FileInterface/RifReaderInterface.h"
//#include "resinsight/ApplicationCode/FileInterface/RifEclipseOutputFileTools.h"

// RESINSIGHT: APPLICATIONCODE/RESERVOIRDATAMODEL ------------------
#include "resinsight/ApplicationCode/ReservoirDataModel/RigMainGrid.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigEclipseCaseData.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigCaseCellResultsData.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigWellPathIntersectionTools.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigHexIntersectionTools.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigEclipseWellLogExtractor.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigActiveCellInfo.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigCompletionData.h"
#include "resinsight/ApplicationCode/ReservoirDataModel/RigWellPath.h"

//#include "resinsight/ApplicationCode/ReservoirDataModel/RigGridBase.h"

//#include "resinsight/ApplicationCode/ReservoirDataModel/RigCell.h"

//#include "resinsight/ApplicationCode/ReservoirDataModel/RigWellPath.h"
////#include "resinsight/ApplicationCode/ProjectDataModel/RimWellPath.h"


// -----------------------------------------------------------------
// For RigWellPath.cpp functions
//#include "resinsight/ApplicationCode/ReservoirDataModel/cvfGeometryTools.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

// -----------------------------------------------------------------
//class RigWellPath : public cvf::Object
//{
// public:
//  RigWellPath();
//  vector<cvf::Vec3d>  m_wellPathPoints;
//  vector<double>      m_measuredDepths;
//
////  cvf::Vec3d
////  interpolatedPointAlongWellPath(double measuredDepth) const;
//
//// std::pair<vector<cvf::Vec3d>, vector<double> >
//// clippedPointSubset(double startMD, double endMD) const;
//// void                setDatumElevation(double value);
//// bool                hasDatumElevation() const;
//// double              datumElevation() const;
//// double              wellPathAzimuthAngle(const cvf::Vec3d& position) const;
//// void                twoClosestPoints(const cvf::Vec3d& position, cvf::Vec3d* p1, cvf::Vec3d* p2) const;
//// std::vector<cvf::Vec3d>
//// wellPathPointsIncludingInterpolatedIntersectionPoint(
//// double intersectionMeasuredDepth) const;
//
// private:
//  bool    m_hasDatumElevation;
//  double  m_datumElevation;
//};

// -----------------------------------------------------------------

//class WellPath {
//
//};

//class MainGrid {
//
//};

// -----------------------------------------------------------------

struct WellData
{
  QString                 m_name;
  cvf::ref<RigWellPath>   m_wellPathGeometry;
};

class wicalc_rins {
 public:
  wicalc_rins(Settings::Model::Well well_settings,
            Grid::Grid *grid);
  ~wicalc_rins();

  Settings::Model::Well well_settings_;

//  void calculateWellPathIntersections(const RimWellPath*   wellPath,
//                                      const MainGrid*      grid,
//                                      std::vector<double>& values);

  CellDirection calculateDirectionInCell(RigEclipseCaseData* eclipseCase,
                                         size_t cellIndex,
                                         const cvf::Vec3d& lengthsInCell);

  double calculateTransmissibility(RigEclipseCaseData* eclipseCase,
                                   const RigWellPath* wellPath,
                                   const cvf::Vec3d& internalCellLengths,
                                   double skinFactor,
                                   double wellRadius,
                                   size_t cellIndex,
                                   bool useLateralNTG,
                                   size_t volumeScaleConstant = 1,
                                   CellDirection directionForVolumeScaling =
                                   CellDirection::DIR_I);

  double totalConnectionFactor(double transX,
                               double transY,
                               double transZ);

  double wellBoreTransmissibilityComponent(double cellPerforationVectorComponent,
                                           double permeabilityNormalDirection1,
                                           double permeabilityNormalDirection2,
                                           double cellSizeNormalDirection1,
                                           double cellSizeNormalDirection2,
                                           double wellRadius,
                                           double skinFactor,
                                           double cDarcyForRelevantUnit);

  double peacemanRadius(double permeabilityNormalDirection1,
                        double permeabilityNormalDirection2,
                        double cellSizeNormalDirection1,
                        double cellSizeNormalDirection2);

  void print_dbg_msg(string dbg_str,
                     int vlevel,
                     Eigen::VectorXd eigvec = Eigen::VectorXd(0));



 private:
  void calculateWellPathIntersections(const RigWellPath *wellPath,
                                      const RigMainGrid *grid,
                                      vector<double> &values);

  void generatePerforationsCompdatValues(const RigWellPath* wellPath);

  //
  RigEclipseCaseData *reservoir2_;
  RifReaderEclipseOutput* readerInterfaceEcl_;
  RigWellPath *wellPath_;
  RigMainGrid *grid_;
  vector<double> values_;

};

}
}

#endif //FIELDOPT_WICALC_RI_H
