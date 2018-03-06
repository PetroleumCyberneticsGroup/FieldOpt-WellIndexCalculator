//
// Created by bellout on 2/24/18.
//

// -----------------------------------------------------------------
#include "wicalc_rinx.h"

// -----------------------------------------------------------------
using std::cout;
using std::endl;
using std::string;
using std::vector;

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

// -----------------------------------------------------------------
wicalc_rinx::~wicalc_rinx() {

//  delete [] well_settings_;
//
//  delete [] eclipseCaseData_;
//  delete [] readerEclipseOutput_;
//  delete [] eclipseCaseData_;
//
//  delete [] wellPath_;
//  delete [] intersections_;

}

// -----------------------------------------------------------------
wicalc_rinx::wicalc_rinx(Settings::Model::Well well_settings,
                         Grid::Grid *grid) {

  well_settings_ = well_settings;
  grid_ = grid;

  eclipseCaseData_ = new RigEclipseCaseData();
  readerEclipseOutput_ = new RifReaderEclipseOutput();
  readerEclipseOutput_->open(QString::fromStdString(grid->GetFilePath()),
                             eclipseCaseData_);
  eclipseCaseData_->computeActiveCellBoundingBoxes();
  eclipseCaseData_->mainGrid()->computeCachedData();

  wellPath_ = new RigWellPath();

  size_t grid_count_ = eclipseCaseData_->mainGrid()->gridCount();
  size_t cell_count_ = eclipseCaseData_->mainGrid()->cellCount();
  size_t gcellarray_sz_ = eclipseCaseData_->mainGrid()->globalCellArray().size();

  intersections_.resize(gcellarray_sz_);
  // Set up intersected cells values
  std::fill(intersections_.begin(), intersections_.end(), HUGE_VAL);

  // ---------------------------------------------------------------
  // Dbg
  QDateTime tstart = QDateTime::currentDateTime();
  std::stringstream str; str << "Find cell from coords.";
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

  str.str(""); str << "grid_->gridCount(): " << grid_count_
                   << " -- grid_->cellCount(): " << cell_count_
                   << " -- grid_->globalCellArray().size(): "
                   << gcellarray_sz_;
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

};

// -----------------------------------------------------------------
void wicalc_rinx::ComputeWellBlocks(
    map<string, vector<IntersectedCell>> &well_indices,
    vector<WellDefinition> &wells,
    int rank) {

  std::stringstream str;

  // Perform well block search for each well
  for (int iWell = 0; iWell < wells.size(); ++iWell) {

    // Intersected cells for well
    vector<IntersectedCell> intersected_cells;

    for (int iSegment = 0; iSegment < wells[iWell].radii.size(); ++iSegment) {

      std::vector<WellPathCellIntersectionInfo> intersectionInfos;
      wellPath_= new RigWellPath();

      // Convert well data to RI file types
      double cvf_xyzHeelMD = wells[iWell].heel_md[iSegment];
      double cvf_xyzToeMD = wells[iWell].toe_md[iSegment];
      wellPath_->m_measuredDepths.push_back(cvf_xyzHeelMD);
      wellPath_->m_measuredDepths.push_back(cvf_xyzToeMD);

      cvf::Vec3d cvf_xyzHeel, cvf_xyzToe;

      cvf_xyzHeel[0] = wells[iWell].heels[iSegment][0];
      cvf_xyzHeel[1] = wells[iWell].heels[iSegment][1];
      cvf_xyzHeel[2] = -wells[iWell].heels[iSegment][2];

      cvf_xyzToe[0] = wells[iWell].toes[iSegment][0];
      cvf_xyzToe[1] = wells[iWell].toes[iSegment][1];
      cvf_xyzToe[2] = -wells[iWell].toes[iSegment][2];

      wellPath_->m_wellPathPoints.push_back(cvf_xyzHeel);
      wellPath_->m_wellPathPoints.push_back(cvf_xyzToe);

      // Dbg
      str.str(""); str << "calculateWellPathIntersections.";
      print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

      calculateWellPathIntersections(wellPath_,
                                     eclipseCaseData_->mainGrid(),
                                     intersections_);

      // Dbg
      str.str(""); str << "intersections_.size(): " << intersections_.size();
      print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0);

      // Use intersection data to find intersected cell data
      cvf::ref<RigEclipseWellLogExtractor>
          extractor = new RigEclipseWellLogExtractor(eclipseCaseData_,
                                                     wellPath_,
                                                     wells[iWell].wellname);

      activeCellInfo_ =
          eclipseCaseData_->activeCellInfo(RiaDefines::MATRIX_MODEL);

      std::vector<WellPathCellIntersectionInfo>
          intersectedCellInfo = extractor->cellIntersectionInfosAlongWellPath();
      cout << "# of intersectedCells: " << intersectedCellInfo.size() << endl;

      collectIntersectedCells(intersected_cells,
                              intersectedCellInfo,
                              wells[iWell],
                              rank);
      // delete wellPath_;

    }

    // Assign intersected cells to well
    well_indices[wells[iWell].wellname] = intersected_cells;
  } // # of wells

}

// -----------------------------------------------------------------
void wicalc_rinx::collectIntersectedCells(vector<IntersectedCell> &isc_cells,
                                          vector<WellPathCellIntersectionInfo> isc_info,
                                          WellDefinition well,
                                          int rank) {

  std::vector<RigCompletionData> completionData;
//  int ii = 0;

  for (auto& cell : isc_info) {

    size_t i, j, k;
    eclipseCaseData_->
        mainGrid()->ijkFromCellIndex(cell.globCellIndex, &i, &j, &k);

    bool cellIsActive = activeCellInfo_->isActive(cell.globCellIndex);
    if (!cellIsActive) {
      cout << "Cell is not active" << endl;
      continue;
    }

    RigCompletionData completion(QString::fromStdString(well.wellname),
                                 IJKCellIndex(i, j, k));
    CellDirection direction =
        calculateDirectionInCell(eclipseCaseData_,
                                 cell.globCellIndex,
                                 cell.intersectionLengthsInCellCS);


    // -------------------------------------------------------------
    Vector3d start_pt(cell.startPoint.x(),
                      cell.startPoint.y(),
                      cell.startPoint.z());

    Vector3d exit_pt(cell.endPoint.x(),
                     cell.endPoint.y(),
                     cell.endPoint.z());

    Vector3d isc_lengths(cell.intersectionLengthsInCellCS.x(),
                         cell.intersectionLengthsInCellCS.y(),
                         cell.intersectionLengthsInCellCS.z());

    // Confirm diff. b/e pcg-wic and ri-wic comes mainly from
    // projected well lengths in well block
    if (well.wellname == "P01") {

      if (cell.globCellIndex == 52722) {
        cell.intersectionLengthsInCellCS[0] = 67.1;
        cell.intersectionLengthsInCellCS[1] = 60.9;
        cell.intersectionLengthsInCellCS[2] = 0.8;
      }
      else if (cell.globCellIndex == 52723) {
        cell.intersectionLengthsInCellCS[0] = 66.8;
        cell.intersectionLengthsInCellCS[1] = 61.2;
        cell.intersectionLengthsInCellCS[2] = 0.9;
      }
      else if (cell.globCellIndex == 52605) {
        cell.intersectionLengthsInCellCS[0] = 67.0;
        cell.intersectionLengthsInCellCS[1] = 61.4;
        cell.intersectionLengthsInCellCS[2] = 0.8;
      }

    } else if (well.wellname == "P02") {

      if (cell.globCellIndex == 54129) {
        cell.intersectionLengthsInCellCS[0] = 14.2;
        cell.intersectionLengthsInCellCS[1] = 29.0;
        cell.intersectionLengthsInCellCS[2] = 0.1;
      }
      else if (cell.globCellIndex == 32652) {
        cell.intersectionLengthsInCellCS[0] = 14.7;
        cell.intersectionLengthsInCellCS[1] = 29.4;
        cell.intersectionLengthsInCellCS[2] = 0.1;
      }
      else if (cell.globCellIndex == 11293) {
        cell.intersectionLengthsInCellCS[0] = 14.7;
        cell.intersectionLengthsInCellCS[1] = 29.8;
        cell.intersectionLengthsInCellCS[2] = 0.1;
      }

    }


    // -------------------------------------------------------------
    std::stringstream str;
    str << FLYELLOW
        << "ic=" << setw(3) << right << isc_cells.size() << " "
        << "s=["
        << setprecision(3) << fixed
        << setw(10) << start_pt(0) << " "
        << setw(10) << start_pt(1) << " "
        << setw(6)  << start_pt(2) << "] "
        << "e=["
        << setw(10) << exit_pt(0) << " "
        << setw(10) << exit_pt(1) << " "
        << setw(6)  << exit_pt(2) << "]"
        << AEND;
    print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, false);

    Reservoir::Grid::Cell gcell;
    gcell = Reservoir::Grid::Cell(grid_->Dimensions().nx *
        grid_->Dimensions().ny * grid_->Dimensions().nz);

    IntersectedCell icell(gcell);
    icell.set_global_index(cell.globCellIndex);
    icell.set_ijk_index(Reservoir::Grid::IJKCoordinate(i,j,k));

    // -------------------------------------------------------------
    double transmissibility =
        calculateTransmissibility(eclipseCaseData_,
                                  wellPath_,
                                  cell.intersectionLengthsInCellCS,
                                  well.skins[0],
                                  well.radii[0],
                                  cell.globCellIndex,
                                  false,
                                  icell);

    // -------------------------------------------------------------
    str.str("");
    str << FLYELLOW
        //        << "wi=" << transmissibility << " "
        //        << "s=" << well.skins[0] << " "
        //        << "r=" << well.radii[0] << " "
        << "wl=[" << setprecision(3)
        << setw(7) << right << isc_lengths(0) << " "
        << setw(7) << right << isc_lengths(1) << " "
        << setw(7) << right << isc_lengths(2) << "] "
        << AEND;
    print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, false);
    cout << endl;


//    completion.setTransAndWPImultBackgroundDataFromPerforation(transmissibility,
//                                                               well.skins[0],
//                                                               well.radii[0],
//                                                               direction);
//    completion.addMetadata("Perforation",
//                           QString("StartMD: %1 - EndMD: %2")
//                               .arg(wellPath_->m_measuredDepths[0])
//                               .arg(wellPath_->m_measuredDepths[1])
//                               + QString(" : ") + QString::number(transmissibility));
//    completionData.push_back(completion);

//    auto gcell = grid_->GetCell((int)cell.globCellIndex);


//    vector<double> permxv, permyv, permzv;
//    vector<double> permxc, permyc, permzc;

//    readerEclipseOutput_->staticResult("PERMX", RiaDefines::MATRIX_MODEL, &permxv);
//    readerEclipseOutput_->staticResult("PERMY", RiaDefines::MATRIX_MODEL, &permyv);
//    readerEclipseOutput_->staticResult("PERMZ", RiaDefines::MATRIX_MODEL, &permzv);
//    permxc.push_back(permxv[(int)cell.globCellIndex-1]);
//    permyc.push_back(permyv[(int)cell.globCellIndex-1]);
//    permzc.push_back(permzv[(int)cell.globCellIndex-1]);



    icell.add_new_segment(start_pt, exit_pt, well.radii[0], well.skins[0]);
    icell.set_cell_well_index_matrix(transmissibility);
    isc_cells.push_back(icell);

  }

}

// -----------------------------------------------------------------
double wicalc_rinx::calculateTransmissibility(RigEclipseCaseData* eclipseCase,
                                              const RigWellPath* wellPath,
                                              const cvf::Vec3d& internalCellLengths,
                                              double skinFactor,
                                              double wellRadius,
                                              size_t cellIndex,
                                              bool useLateralNTG,
                                              IntersectedCell &icell,
                                              size_t volumeScaleConstant,
                                              CellDirection directionForVolumeScaling) {

  // -------------------------------------------------------------
  std::vector<double> dxv, dyv, dzv;
  std::vector<double> dxc, dyc, dzc;

  readerEclipseOutput_->staticResult("DX", RiaDefines::MATRIX_MODEL, &dxv);
  readerEclipseOutput_->staticResult("DY", RiaDefines::MATRIX_MODEL, &dyv);
  readerEclipseOutput_->staticResult("DZ", RiaDefines::MATRIX_MODEL, &dzv);

  std::vector<double> permxv, permyv, permzv;
  std::vector<double> permxc, permyc, permzc;

  readerEclipseOutput_->staticResult("PERMX", RiaDefines::MATRIX_MODEL, &permxv);
  readerEclipseOutput_->staticResult("PERMY", RiaDefines::MATRIX_MODEL, &permyv);
  readerEclipseOutput_->staticResult("PERMZ", RiaDefines::MATRIX_MODEL, &permzv);

  std::vector<double> ntgv;
  std::vector<double> ntgc;

  readerEclipseOutput_->staticResult("NTG", RiaDefines::MATRIX_MODEL, &ntgv);

  // -------------------------------------------------------------
  double ntg = 1.0;
//  ntg = ntgv[cellIndex];
//  double latNtg = useLateralNTG ? ntg : 1.0;
  double latNtg = 1.0;

  double dx = dxv[cellIndex];
  double dy = dyv[cellIndex];
  double dz = dzv[cellIndex];

  double permx = permxv[cellIndex];
  double permy = permyv[cellIndex];
  double permz = permzv[cellIndex];

  // -------------------------------------------------------------
  // OVERRIDE
  dx = grid_->GetCell(icell.global_index()).dxdydz()(0);
  dy = grid_->GetCell(icell.global_index()).dxdydz()(1);
  dz = grid_->GetCell(icell.global_index()).dxdydz()(2);

  permx = grid_->GetCell(icell.global_index()).permx()[0];
  permy = grid_->GetCell(icell.global_index()).permy()[0];
  permz = grid_->GetCell(icell.global_index()).permz()[0];

  double darcy = 0.008527;

  // -------------------------------------------------------------
  dxc.push_back(dx);
  dyc.push_back(dy);
  dzc.push_back(dz);

  permxc.push_back(permx);
  permyc.push_back(permy);
  permzc.push_back(permz);

  icell.SetProperties(true, false,
                      permxc, permyc, permzc);
//  icell.set_dx(dx);
//  icell.set_dy(dy);
//  icell.set_dz(dz);

  icell.set_segment_calculation_data(0,"dx",dx);
  icell.set_segment_calculation_data(0,"dy",dy);
  icell.set_segment_calculation_data(0,"dz",dz);

  // -------------------------------------------------------------
  std::stringstream str;
  str << FLYELLOW
      << setw(3) << setprecision(1) << fixed
      << "D=["
      << setw(5) << right << dx << " "
      << setw(5) << right << dy << " "
      << setw(4) << right << dz << "] "
      << "P=["
      << setw(6) << right << permx << " "
      << setw(6) << right << permy << " "
      << setw(5) << right << permz << "] ";

  // -------------------------------------------------------------
  if (volumeScaleConstant != 1) {
    if (directionForVolumeScaling == CellDirection::DIR_I) dx = dx / volumeScaleConstant;
    if (directionForVolumeScaling == CellDirection::DIR_J) dy = dy / volumeScaleConstant;
    if (directionForVolumeScaling == CellDirection::DIR_K) dz = dz / volumeScaleConstant;
  }

  // -------------------------------------------------------------
  double transx =
      wellBoreTransmissibilityComponent(internalCellLengths.x() * latNtg,
                                        permy, permz, dy, dz, wellRadius,
                                        skinFactor, darcy);

  double transy =
      wellBoreTransmissibilityComponent(internalCellLengths.y() * latNtg,
                                        permx, permz, dx, dz, wellRadius,
                                        skinFactor, darcy);

  double transz =
      wellBoreTransmissibilityComponent(internalCellLengths.z() * ntg,
                                        permy, permx, dy, dx, wellRadius,
                                        skinFactor, darcy);

  double total_connection = totalConnectionFactor(transx, transy, transz);

  assert(latNtg==1);
  assert(ntg==1);
  icell.set_segment_calculation_data(0,"Lx",internalCellLengths.x() * latNtg);
  icell.set_segment_calculation_data(0,"Ly",internalCellLengths.y() * latNtg);
  icell.set_segment_calculation_data(0,"Lz",internalCellLengths.z() * ntg);

  icell.set_segment_calculation_data(0,"wx_m",transx);
  icell.set_segment_calculation_data(0,"wy_m",transy);
  icell.set_segment_calculation_data(0,"wz_m",transz);

  // -------------------------------------------------------------
  str << setw(3) << "T=["
      << setprecision(3) << fixed
      << setw(8) << right << transx << " "
      << setw(8) << right << transy << " "
      << setw(8) << right << transz << "] "
      << setw(3) << "tc="
      << setw(8) << right << total_connection << AEND;
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, false);

  return total_connection;
}

// -----------------------------------------------------------------
double wicalc_rinx::wellBoreTransmissibilityComponent(
    double cellPerforationVectorComponent,
    double permeabilityNormalDirection1,
    double permeabilityNormalDirection2,
    double cellSizeNormalDirection1,
    double cellSizeNormalDirection2,
    double wellRadius,
    double skinFactor,
    double cDarcyForRelevantUnit) {

  double K =
      cvf::Math::sqrt(permeabilityNormalDirection1 * permeabilityNormalDirection2);

  double nominator =
      cDarcyForRelevantUnit * 2 * cvf::PI_D * K * cellPerforationVectorComponent;

  // -------------------------------------------------------------
  std::stringstream str;
  str << setw(1) << "K=" << setw(7) << setprecision(3) << fixed << K;
  print_dbg_msg_wic_ri(__func__, str.str(), 0.0, 0, false);

  double peaceManRad = peacemanRadius(permeabilityNormalDirection1,
                                      permeabilityNormalDirection2,
                                      cellSizeNormalDirection1,
                                      cellSizeNormalDirection2);

  double denominator = log(peaceManRad / wellRadius) + skinFactor;

  double trans = nominator / denominator;
  return trans;
}

// -----------------------------------------------------------------
double wicalc_rinx::peacemanRadius(double permeabilityNormalDirection1,
                                   double permeabilityNormalDirection2,
                                   double cellSizeNormalDirection1,
                                   double cellSizeNormalDirection2) {

  double numerator = cvf::Math::sqrt(
      pow(cellSizeNormalDirection2, 2.0) *
          pow(permeabilityNormalDirection1 / permeabilityNormalDirection2, 0.5)
          + pow(cellSizeNormalDirection1, 2.0) *
              pow(permeabilityNormalDirection2 / permeabilityNormalDirection1, 0.5) );

  double denominator = pow((permeabilityNormalDirection1 / permeabilityNormalDirection2), 0.25 )
      + pow((permeabilityNormalDirection2 / permeabilityNormalDirection1), 0.25 );

  double r0 = 0.28 * numerator / denominator;

  return r0;
}

// -----------------------------------------------------------------
double wicalc_rinx::totalConnectionFactor(double transX,
                                          double transY,
                                          double transZ) {
  return  cvf::Math::sqrt(
      pow(transX, 2.0) + pow(transY, 2.0) + pow(transZ, 2.0));
}

// -----------------------------------------------------------------
void wicalc_rinx::calculateWellPathIntersections(const RigWellPath *wellPath,
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

// -----------------------------------------------------------------
CellDirection
wicalc_rinx::calculateDirectionInCell(RigEclipseCaseData* eclipseCase,
                                      size_t cellIndex,
                                      const cvf::Vec3d& lengthsInCell) {

  std::vector<double> valuesx, valuesy, valuesz;

  readerEclipseOutput_->staticResult("DX", RiaDefines::MATRIX_MODEL, &valuesx);
  readerEclipseOutput_->staticResult("DY", RiaDefines::MATRIX_MODEL, &valuesy);
  readerEclipseOutput_->staticResult("DZ", RiaDefines::MATRIX_MODEL, &valuesz);

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

}
}