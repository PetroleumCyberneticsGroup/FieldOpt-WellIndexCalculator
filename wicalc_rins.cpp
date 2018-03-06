//
// Created by bellout on 2/13/18.
//

// -----------------------------------------------------------------
#include "wicalc_rins.h"
#include <Utilities/time.hpp>

using std::cout;
using std::endl;
using std::string;
using std::vector;

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

wicalc_rins::wicalc_rins(Settings::Model::Well well_settings,
                         Grid::Grid *grid) {

  well_settings_ = well_settings;

//  Eigen::Vector3d xyzHeel, xyzToe;
//  double xyzHeelTVD, xyzToeTVD;
//  double xyzHeelMD, xyzToeMD, wellLength;
//
//  xyzHeel[0] = 5.2482537021119706e+05;
//  xyzHeel[1] = 6.1802875116297025e+06;
//  xyzHeel[2] = 2.0508774871826172e+03;
////  xyzHeelTVD = xyzHeel[2];
//  xyzHeelMD = xyzHeel[2];
//
//  xyzToe[0] = 5.2683537021119706e+05,
//  xyzToe[1] = 6.2002875116297025e+06,
//  xyzToe[2] = 2.0608774871826172e+03;
////  xyzToeTVD = xyzToe[2];
//  wellLength = sqrt((xyzToe - xyzHeel).norm());
//  xyzToeMD = xyzHeelMD + wellLength;

  // Create the data container
//  vector<WellData> fileWellDataArray;
//
//  fileWellDataArray.push_back(WellData());
//  fileWellDataArray.back().m_wellPathGeometry = new RigWellPath();
//
//  cvf::Vec3d wellPointHeel(xyzHeel[0], xyzHeel[1], xyzHeelTVD);
//  fileWellDataArray.back().m_wellPathGeometry->m_wellPathPoints.push_back(wellPointHeel);
//  fileWellDataArray.back().m_wellPathGeometry->m_measuredDepths.push_back(xyzHeelMD);
//
//  cvf::Vec3d wellPointToe(xyzHeel[0], xyzHeel[1], xyzToeTVD);
//  fileWellDataArray.back().m_wellPathGeometry->m_wellPathPoints.push_back(wellPointToe);
//  fileWellDataArray.back().m_wellPathGeometry->m_measuredDepths.push_back(xyzToeMD);

//  RigEclipseCaseData* eclipseData = new RigEclipseCaseData();
//  grid_ = new RigMainGrid();
//  eclipseData->setMainGrid(grid_);
//
//
////  wellPath_ = fileWellDataArray[0].m_wellPathGeometry.;
//
//  cvf::ref<RifReaderEclipseOutput> readerEclipseOutput = new RifReaderEclipseOutput;
//
//  intersections_ = vector<double>();
//
//  calculateWellPathIntersections(wellPath_, grid_, intersections_);
//
//  for (unsigned int i = 0; i < intersections_.size(); i++) {
//    cout << intersections_[i] << " " << endl;
//    if (i % 10 == 1) cout << endl;
//  }

//#include <Utilities/time.hpp>
//  time_t start, end;
//  time(&start);
//  time(&end);


//  cvf::ref<RigEclipseCaseData> reservoir = new RigEclipseCaseData();

  RigEclipseCaseData* reservoir2 = new RigEclipseCaseData();
  reservoir2_ = reservoir2;

//  cvf::ref<RifReaderEclipseOutput>
//      readerInterfaceEcl = new RifReaderEclipseOutput();

  RifReaderEclipseOutput* readerInterfaceEcl = new RifReaderEclipseOutput();
  readerInterfaceEcl_ = readerInterfaceEcl;

  QString filename("/home/bellout/WORK-3/SCRATCH_RUNS-PROJECTS/P11_OLYMPUS-cases+work_scratch-runs/xruns/xruns___EEE/INPUT/OLYMPUS.EGRID");
  readerInterfaceEcl_->open(filename, reservoir2_);

  reservoir2_->computeActiveCellBoundingBoxes();
  reservoir2_->mainGrid()->computeCachedData();
  grid_ = reservoir2_->mainGrid();
//  grid_->computeCachedData();


  // ---------------------------------------------------------------
  QDateTime tstart = QDateTime::currentDateTime();
  std::stringstream str; str << "Find cell from coords.";
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

//  grid_ = reservoir.p()->mainGrid();
//  grid_ = reservoir2->mainGrid();
//  grid_->computeCachedData();


  auto g_count = reservoir2->mainGrid()->gridCount();
  auto c_count = reservoir2->mainGrid()->cellCount();
  auto carr_sz = reservoir2->mainGrid()->globalCellArray().size();

  str.str(""); str << "grid_->gridCount(): " << g_count
                   << " -- grid_->cellCount(): " << c_count
                   << " -- grid_->globalCellArray().size(): "
                   << carr_sz;
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

  cvf::Vector3<double> cvf_xyzHeel, cvf_xyzToe;
  double cvf_xyzHeelTVD, cvf_xyzToeTVD;
  double cvf_xyzHeelMD, cvf_xyzToeMD, cvf_wellLength;


//  polygonExample.push_back(cvf::Vec3d(0.00, 0.00, 0.0));

//  RifReaderSettings readerSettings;
//  readerInterfaceEcl->setReaderSetting(&readerSettings);

  //  RigMainGrid gridd = new RigMainGrid();



//  QStringList staticResults = readerInterfaceEcl->staticResults();
//  qDebug() << "Static results\n" << staticResults;
//
//  QStringList dynamicResults = readerInterfaceEcl->dynamicResults();
//  qDebug() << "Dynamic results\n" << dynamicResults;

//  size_t idx = 0;
//  reservoir2->results(RiaDefines::MATRIX_MODEL)->cellScalarResults(idx);

//  size_t dxResultGridIndex =
//      reservoir2->results(RiaDefines::MATRIX_MODEL)->
//          findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DX");
//
//
//  size_t scalarResultIndexDX =
//      reservoir2->results(RiaDefines::MATRIX_MODEL)->
//          findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DX");
//
//  const std::vector<double>* dxResults =
//      &(reservoir2->results(RiaDefines::MATRIX_MODEL)->
//          cellScalarResults(scalarResultIndexDX, 0));

//  cvf_xyzHeel[0] = 5.2482537021119706e+05 - 000;
//  cvf_xyzHeel[1] = 6.1802875116297025e+06;
//  cvf_xyzHeel[2] = -2.0508774871826172e+03;
//  cvf_xyzHeelMD = abs(cvf_xyzHeel[2]);
//
//  cvf_xyzToe[0] = 5.2482537021119706e+05 + 100;
//  cvf_xyzToe[1] = 6.1802875116297025e+06;
//  cvf_xyzToe[2] = -2.0508774871826172e+03;

  cvf_xyzHeel[0] = 5.2304999558394222e+05;
  cvf_xyzHeel[1] = 6.1801758504716540e+06;
  cvf_xyzHeel[2] = -2.0654362487792969e+03;
  cvf_xyzHeelMD = abs(cvf_xyzHeel[2]);

  cvf_xyzToe[0] = 5.2305999558394222e+05;
  cvf_xyzToe[1] = 6.1802758504716540e+06;
  cvf_xyzToe[2] = -2.0654362487792969e+03;

  auto cvf_diff = cvf_xyzHeel;
  cvf_diff.subtract(cvf_xyzToe);
  cvf_wellLength = sqrt(cvf_diff.dot(cvf_diff));
  cvf_xyzToeMD = cvf_xyzHeelMD + cvf_wellLength;



  std::vector<WellPathCellIntersectionInfo> intersectionInfos;

  RigWellPath *wellPath = new RigWellPath();
  wellPath_ = wellPath;

  wellPath_->m_measuredDepths.push_back(cvf_xyzHeelMD);
  wellPath_->m_measuredDepths.push_back(cvf_xyzToeMD);

  wellPath_->m_wellPathPoints.push_back(cvf_xyzHeel);
  wellPath_->m_wellPathPoints.push_back(cvf_xyzToe);

  values_.resize(carr_sz);
  std::fill(values_.begin(), values_.end(), HUGE_VAL);

  str.str(""); str << "calculateWellPathIntersections.";
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

  calculateWellPathIntersections(wellPath_, grid_, values_);

  str.str(""); str << "intersections_.size(): " << values_.size();
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

//  cvf::ref<RigWellPath> dummyWellPath = new RigWellPath;

//  const RigWellPath* testpath = new RigWellPath;
//  testpath = wellPath_;

//  dummyWellPath->m_wellPathPoints = wellPath_->m_wellPathPoints;
//  dummyWellPath->m_measuredDepths = wellPath_->m_measuredDepths;

  cout << "startMD: " << wellPath_->m_measuredDepths[0] << " -- "
       << "endMD: " << wellPath_->m_measuredDepths[1] << " === ";

  cout << "startPoint: " << wellPath_->m_wellPathPoints[0].x()
       << " " << wellPath_->m_wellPathPoints[0].y()
       << " " << wellPath_->m_wellPathPoints[0].z();

  cout << " -- endPoint: " << wellPath_->m_wellPathPoints[1].x()
       << " " << wellPath_->m_wellPathPoints[1].y()
       << " " << wellPath_->m_wellPathPoints[1].z() << endl;

  const std::string& test = "test";
  cvf::ref<RigEclipseWellLogExtractor>
      extractor = new RigEclipseWellLogExtractor(reservoir2_,
                                                 wellPath_,
                                                 test);

//  RigEclipseWellLogExtractor(
//  const RigEclipseCaseData* aCase,
//  const RigWellPath* wellpath,
//  const std::string& wellCaseErrorMsgName)

  std::vector<RigCompletionData> completionData;

  std::vector<WellPathCellIntersectionInfo>
      intersectedCells = extractor->cellIntersectionInfosAlongWellPath();

  cout << "# of intersectedCells: " << intersectedCells.size() << endl;


  const RigActiveCellInfo* activeCellInfo =
      reservoir2_->activeCellInfo(RiaDefines::MATRIX_MODEL);

  double skinFactor = 0.0;
  double diameter = 2*0.1905;

  for (auto& cell : intersectedCells) {

    size_t i, j, k;
    grid_->ijkFromCellIndex(cell.globCellIndex, &i, &j, &k);

    cout << "i: " << static_cast<int>(i)
         << "  j: " << static_cast<int>(j)
         << "  k: " << static_cast<int>(k) << endl;

    cout << "startMD: " << cell.startMD << " -- "
         << "endMD: " << cell.endMD << " === ";

    cout << "startPoint: "
         << cell.startPoint.x() << " " << cell.startPoint.y() << " " << cell.startPoint.z();
    cout << " -- endPoint: "
         << cell.endPoint.x() << " " << cell.endPoint.y() << " " << cell.endPoint.z() << endl;

    bool cellIsActive = activeCellInfo->isActive(cell.globCellIndex);
    if (!cellIsActive) {
      cout << "Cell is not active" << endl;
      continue;
    }

    RigCompletionData completion("TESTWELL", IJKCellIndex(i, j, k));

    CellDirection direction =
        calculateDirectionInCell(reservoir2_,
                                 cell.globCellIndex,
                                 cell.intersectionLengthsInCellCS);


    double transmissibility = calculateTransmissibility(reservoir2_,
                                                        wellPath_,
                                                        cell.intersectionLengthsInCellCS,
                                                        skinFactor,
                                                        diameter,
                                                        cell.globCellIndex,
                                                        false);

    completion.setTransAndWPImultBackgroundDataFromPerforation(transmissibility,
                                                               skinFactor,
                                                               diameter,
                                                               direction);
    completion.addMetadata("Perforation",
                           QString("StartMD: %1 - EndMD: %2")
                               .arg(wellPath_->m_measuredDepths[0])
                               .arg(wellPath_->m_measuredDepths[1])
                               + QString(" : ") + QString::number(transmissibility));
    completionData.push_back(completion);

  }

  for (auto& comp : completionData) {
    cout << std::fixed << std::setprecision(3)
         << "comp.transmissibility() = " << comp.transmissibility() << endl;
  }


}

wicalc_rins::~wicalc_rins() {}

void wicalc_rins::generatePerforationsCompdatValues(const RigWellPath* wellPath) {


//  std::vector<WellPathCellIntersectionInfo> intersectedCells =
//      RigWellPathIntersectionTools::findCellIntersectionInfosAlongPath(
//          settings.caseToApply->eclipseCaseData(),
//          perforationPointsAndMD.first,
//          perforationPointsAndMD.second);

//  std::vector<WellPathCellIntersectionInfo> intersectedCells =
//      RigWellPathIntersectionTools::findCellIntersectionInfosAlongPath(
//          settings.caseToApply->eclipseCaseData(),
//          perforationPointsAndMD.first,
//          perforationPointsAndMD.second);

//  for (auto& cell : intersectedCells) {
//
//    bool cellIsActive = activeCellInfo->isActive(cell.globCellIndex);
//    if (!cellIsActive) continue;
//
//  }


}

void wicalc_rins::calculateWellPathIntersections(const RigWellPath *wellPath,
                                               const RigMainGrid *grid,
                                               std::vector<double> &values) {

  std::vector<HexIntersectionInfo> intersections =
      RigWellPathIntersectionTools::findRawHexCellIntersections(
          grid, wellPath->m_wellPathPoints);

  std::stringstream str;
  if (intersections.size() > 0) {
    str << "Set values to WELL_PATH. -- # of intersections = " << intersections.size()
      << " -- intersection.m_hexIndex = " << intersections[0].m_hexIndex;
    print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);
  }

  for (auto &intersection : intersections) {

    values[intersection.m_hexIndex] = RiaDefines::WELL_PATH;

    str.str(""); str << "intersection.m_hexIndex = "
                     << intersection.m_hexIndex
                     << " -- values[intersection.m_hexIndex] = "
                     << values[intersection.m_hexIndex];
//    print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);
  }

}

//-----------------------------------------------------------------------------
// RigWellPath.cpp function
//RigWellPath::RigWellPath()
//    : m_hasDatumElevation(false),
//      m_datumElevation(0.0) { }

//-----------------------------------------------------------------------------
// RigWellPath.cpp function
//cvf::Vec3d
//RigWellPath::interpolatedPointAlongWellPath(double measuredDepth) const {
//  cvf::Vec3d wellPathPoint = cvf::Vec3d::ZERO;
//
//  size_t i = 0;
//  while (i < m_measuredDepths.size() && m_measuredDepths.at(i) < measuredDepth ) {
//    i++;
//  }
//
//  if (m_measuredDepths.size() > i) {
//    if (i == 0) {
//      //For measuredDepth same or lower than first point, use this first point
//      wellPathPoint = m_wellPathPoints.at(0);
//    } else {
//      //Do interpolation
//      double stepsize = (measuredDepth - m_measuredDepths.at(i-1)) /
//          (m_measuredDepths.at(i) - m_measuredDepths.at(i - 1));
//      wellPathPoint = m_wellPathPoints.at(i - 1) +
//          stepsize * (m_wellPathPoints.at(i) - m_wellPathPoints.at(i-1));
//    }
//  } else {
//    //Use endpoint if measuredDepth same or higher than last point
//    wellPathPoint = m_wellPathPoints.at(i-1);
//  }
//  return wellPathPoint;
//}

void wicalc_rins::print_dbg_msg(string dbg_str,
                              int vlevel,
                              Eigen::VectorXd eigvec) {

  if (well_settings_.verb_vector_[6] >= vlevel) { // idx:6 -> opt (Optimization)
    if (   dbg_str == "[opt]Init. Abs.Class GSS.---- "
        || dbg_str == "[opt]Generating trial points. "
        || dbg_str == "[opt]Init. CompassSearch.---- "
        ||     dbg_str == "[opt]Launching opt.iteration. ") {
      cout << dbg_str << endl;
    }
  }
}



//void wicalc_rins::calculateWellPathIntersections(const RigWellPath*      wellPath,
//                                            const RimMainGrid*      grid,
//                                            std::vector<double>& values) {
//
//  std::vector <HexIntersectionInfo> intersections =
//      RigWellPathIntersectionTools::findRawHexCellIntersections(
//          grid,
//          wellPath->wellPathGeometry()->m_wellPathPoints);
//}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
CellDirection wicalc_rins::calculateDirectionInCell(RigEclipseCaseData* eclipseCase,
                                       size_t cellIndex,
                                       const cvf::Vec3d& lengthsInCell) {

//  RigEclipseCaseData* eclipseCaseData = eclipseCase->eclipseCaseData();

  std::vector<double> valuesx, valuesy, valuesz;

  readerInterfaceEcl_->staticResult("DX", RiaDefines::MATRIX_MODEL, &valuesx);
  readerInterfaceEcl_->staticResult("DY", RiaDefines::MATRIX_MODEL, &valuesy);
  readerInterfaceEcl_->staticResult("DZ", RiaDefines::MATRIX_MODEL, &valuesz);

//  for (int i = 0; i < 10; i++) {
//    cout << "DX.DY.DZ[i=" << i << "] = "
//         << valuesx[i] << " -- "
//         << valuesy[i] << " -- "
//         << valuesz[i] << endl;
//  }

//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DX");

//  cvf::ref<RigResultAccessor> dxAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCase, 0, RiaDefines::MATRIX_MODEL, 0, "DX");

//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DY");

//  cvf::ref<RigResultAccessor> dyAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCase, 0, RiaDefines::MATRIX_MODEL, 0, "DY");

//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DZ");

//  cvf::ref<RigResultAccessor> dzAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCase, 0, RiaDefines::MATRIX_MODEL, 0, "DZ");

//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->cellScalarResults(cellIndex);

//  double xLengthFraction = fabs(lengthsInCell.x() / dxAccessObject->cellScalarGlobIdx(cellIndex));
//  double yLengthFraction = fabs(lengthsInCell.y() / dyAccessObject->cellScalarGlobIdx(cellIndex));
//  double zLengthFraction = fabs(lengthsInCell.z() / dzAccessObject->cellScalarGlobIdx(cellIndex));

  double xLengthFraction = fabs(lengthsInCell.x() / valuesx[cellIndex]);
  double yLengthFraction = fabs(lengthsInCell.y() / valuesx[cellIndex]);
  double zLengthFraction = fabs(lengthsInCell.z() / valuesx[cellIndex]);

  if (xLengthFraction > yLengthFraction && xLengthFraction > zLengthFraction) {
    return CellDirection::DIR_I;
  }
  else if (yLengthFraction > xLengthFraction && yLengthFraction > zLengthFraction) {
    return CellDirection::DIR_J;
  }
  else {
    return CellDirection::DIR_K;
  }
}


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------

double wicalc_rins::calculateTransmissibility(RigEclipseCaseData* eclipseCase,
                                            const RigWellPath* wellPath,
                                            const cvf::Vec3d& internalCellLengths,
                                            double skinFactor,
                                            double wellRadius,
                                            size_t cellIndex,
                                            bool useLateralNTG,
                                            size_t volumeScaleConstant,
                                            CellDirection directionForVolumeScaling) {

//  RigEclipseCaseData* eclipseCaseData = eclipseCase->eclipseCaseData();

//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DX");
//
//  cvf::ref<RigResultAccessor> dxAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "DX");
//
//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DY");
//
//  cvf::ref<RigResultAccessor> dyAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "DY");
//
//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "DZ");
//
//  cvf::ref<RigResultAccessor> dzAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "DZ");
//
//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "PERMX");
//
//  cvf::ref<RigResultAccessor> permxAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "PERMX");
//
//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "PERMY");
//
//  cvf::ref<RigResultAccessor> permyAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "PERMY");
//
//  eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//      findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "PERMZ");
//
//  cvf::ref<RigResultAccessor> permzAccessObject =
//      RigResultAccessorFactory::createFromUiResultName(
//          eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "PERMZ");

  std::vector<double> valuesx, valuesy, valuesz;

  readerInterfaceEcl_->staticResult("DX", RiaDefines::MATRIX_MODEL, &valuesx);
  readerInterfaceEcl_->staticResult("DY", RiaDefines::MATRIX_MODEL, &valuesy);
  readerInterfaceEcl_->staticResult("DZ", RiaDefines::MATRIX_MODEL, &valuesz);

  std::vector<double> permxv, permyv, permzv;

  readerInterfaceEcl_->staticResult("PERMX", RiaDefines::MATRIX_MODEL, &permxv);
  readerInterfaceEcl_->staticResult("PERMY", RiaDefines::MATRIX_MODEL, &permyv);
  readerInterfaceEcl_->staticResult("PERMZ", RiaDefines::MATRIX_MODEL, &permzv);

  std::vector<double> ntgv;

  readerInterfaceEcl_->staticResult("NTG", RiaDefines::MATRIX_MODEL, &ntgv);

  double ntg = 1.0;

//  size_t ntgResIdx =
//      eclipseCase->results(RiaDefines::MATRIX_MODEL)->
//          findOrLoadScalarResult(RiaDefines::STATIC_NATIVE, "NTG");

//  if (ntgResIdx != cvf::UNDEFINED_SIZE_T) {
//    cvf::ref<RigResultAccessor> ntgAccessObject =
//        RigResultAccessorFactory::createFromUiResultName(
//            eclipseCaseData, 0, RiaDefines::MATRIX_MODEL, 0, "NTG");
//    ntg = ntgAccessObject->cellScalarGlobIdx(cellIndex);
//  }

  ntg = ntgv[cellIndex];
  double latNtg = useLateralNTG ? ntg : 1.0;

//  double dx = dxAccessObject->cellScalarGlobIdx(cellIndex);
//  double dy = dyAccessObject->cellScalarGlobIdx(cellIndex);
//  double dz = dzAccessObject->cellScalarGlobIdx(cellIndex);
//  double permx = permxAccessObject->cellScalarGlobIdx(cellIndex);
//  double permy = permyAccessObject->cellScalarGlobIdx(cellIndex);
//  double permz = permzAccessObject->cellScalarGlobIdx(cellIndex);

  double dx = valuesx[cellIndex];
  double dy = valuesx[cellIndex];
  double dz = valuesx[cellIndex];

  double permx = permxv[cellIndex];
  double permy = permyv[cellIndex];
  double permz = permzv[cellIndex];

  double darcy = 0.008527;

  if (volumeScaleConstant != 1)
  {
    if (directionForVolumeScaling == CellDirection::DIR_I) dx = dx / volumeScaleConstant;
    if (directionForVolumeScaling == CellDirection::DIR_J) dy = dy / volumeScaleConstant;
    if (directionForVolumeScaling == CellDirection::DIR_K) dz = dz / volumeScaleConstant;
  }

  double transx = wellBoreTransmissibilityComponent(internalCellLengths.x() * latNtg, permy, permz, dy, dz, wellRadius, skinFactor, darcy);
  double transy = wellBoreTransmissibilityComponent(internalCellLengths.y() * latNtg, permx, permz, dx, dz, wellRadius, skinFactor, darcy);
  double transz = wellBoreTransmissibilityComponent(internalCellLengths.z() * ntg, permy, permx, dy, dx, wellRadius, skinFactor, darcy);

  return totalConnectionFactor(transx, transy, transz);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
double wicalc_rins::wellBoreTransmissibilityComponent(double cellPerforationVectorComponent,
                                                    double permeabilityNormalDirection1,
                                                    double permeabilityNormalDirection2,
                                                    double cellSizeNormalDirection1,
                                                    double cellSizeNormalDirection2,
                                                    double wellRadius,
                                                    double skinFactor,
                                                    double cDarcyForRelevantUnit)
{
  double K = cvf::Math::sqrt(permeabilityNormalDirection1 * permeabilityNormalDirection2);

  double nominator = cDarcyForRelevantUnit * 2 * cvf::PI_D * K * cellPerforationVectorComponent;

  double peaceManRad = peacemanRadius(permeabilityNormalDirection1,
                                      permeabilityNormalDirection2,
                                      cellSizeNormalDirection1,
                                      cellSizeNormalDirection2);

  double denominator = log(peaceManRad / wellRadius) + skinFactor;

  double trans = nominator / denominator;
  return trans;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
double wicalc_rins::totalConnectionFactor(double transX,
                                        double transY,
                                        double transZ) {
  return  cvf::Math::sqrt(
      pow(transX, 2.0) + pow(transY, 2.0) + pow(transZ, 2.0));
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
double wicalc_rins::peacemanRadius(double permeabilityNormalDirection1,
                                double permeabilityNormalDirection2,
                                 double cellSizeNormalDirection1,
                                 double cellSizeNormalDirection2) {

  double numerator = cvf::Math::sqrt(
      pow(cellSizeNormalDirection2, 2.0) * pow(permeabilityNormalDirection1 / permeabilityNormalDirection2, 0.5)
          + pow(cellSizeNormalDirection1, 2.0) * pow(permeabilityNormalDirection2 / permeabilityNormalDirection1, 0.5) );

  double denominator = pow((permeabilityNormalDirection1 / permeabilityNormalDirection2), 0.25 )
      + pow((permeabilityNormalDirection2 / permeabilityNormalDirection1), 0.25 );

  double r0 = 0.28 * numerator / denominator;

  return r0;
}

}
}

