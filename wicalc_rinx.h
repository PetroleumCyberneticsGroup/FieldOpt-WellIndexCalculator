//
// Created by bellout on 2/24/18.
//

#ifndef FIELDOPT_WICALC_RINX_H
#define FIELDOPT_WICALC_RINX_H

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
#include <Eigen/Dense>

// -----------------------------------------------------------------
// Qt
#include "QString"

// -----------------------------------------------------------------
// FieldOpt
#include "Reservoir/grid/cell.h"
#include "Reservoir/grid/grid.h"
#include "Reservoir/grid/eclgrid.h"
#include "Settings/model.h"
#include "intersected_cell.h"
#include "FieldOpt-WellIndexCalculator/wellindexcalculator.h"

// RESINSIGHT: FWK/VIZFWK/LIBCORE ----------------------------------
#include "resinsight/Fwk/VizFwk/LibCore/cvfBase.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfObject.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfMath.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfVector3.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfCollection.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfAssert.h"
#include "resinsight/Fwk/VizFwk/LibCore/cvfArray.h"

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

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

using std::cout;
using std::endl;
using std::string;
using std::vector;

class wicalc_rinx {
 public:
  wicalc_rinx(Settings::Model::Well well_settings,
              Grid::Grid *grid);
  ~wicalc_rinx();

  Settings::Model::Well well_settings_;
  vector<double> intersections_;
  Grid::Grid* grid_;

  RigEclipseCaseData *eclipseCaseData_;
  RifReaderEclipseOutput* readerEclipseOutput_;
  RigWellPath *wellPath_;
  // RigMainGrid *mainGrid_;

  const RigActiveCellInfo* activeCellInfo_;

  void ComputeWellBlocks(map<string, vector<IntersectedCell>> &well_indices,
                         vector<WellDefinition> &wells, int rank = 0);

  void collectIntersectedCells(vector<IntersectedCell> &isc_cells,
                               vector<WellPathCellIntersectionInfo> isc_info,
                               WellDefinition well,
                               int rank);

  double calculateTransmissibility(RigEclipseCaseData* eclipseCase,
                                   const RigWellPath* wellPath,
                                   const cvf::Vec3d& internalCellLengths,
                                   double skinFactor,
                                   double wellRadius,
                                   size_t cellIndex,
                                   bool useLateralNTG,
                                   IntersectedCell &icell,
                                   size_t volumeScaleConstant = 1,
                                   CellDirection directionForVolumeScaling =
                                   CellDirection::DIR_I);

  double peacemanRadius(double permeabilityNormalDirection1,
                        double permeabilityNormalDirection2,
                        double cellSizeNormalDirection1,
                        double cellSizeNormalDirection2);

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

  void calculateWellPathIntersections(const RigWellPath *wellPath,
                                      const RigMainGrid *grid,
                                      std::vector<double> &values);

  CellDirection calculateDirectionInCell(RigEclipseCaseData* eclipseCase,
                                         size_t cellIndex,
                                         const cvf::Vec3d& lengthsInCell);

 protected:
  int grid_count_;
  int cell_count_;
  int gcellarray_sz_;

};

}
}

#endif //FIELDOPT_WICALC_RINX_H
