//
// Created by bellout on 3/2/18.
//

#ifndef FIELDOPT_WICALC_RIXX_H
#define FIELDOPT_WICALC_RIXX_H

// -----------------------------------------------------------------
// FieldOpt::RESINXX
#include "resinxx/well_path.h"

// -----------------------------------------------------------------
namespace Reservoir {
namespace WellIndexCalculation {

using std::cout;
using std::endl;
using std::string;
using std::vector;

//====================================================================
class wicalc_rixx {
 public:
  wicalc_rixx(Settings::Model::Well well_settings,
  Grid::Grid *grid);
  ~wicalc_rixx();

  Settings::Model::Well well_settings_;
  vector<double> intersections_;
  Grid::Grid* grid_;
  RIGrid* RIGrid_;

  WellPath *wellPath_;

  void calculateWellPathIntersections(const WellPath *wellPath,
                                      const RIGrid *grid,
                                      vector<double> &isc_values);

  void ComputeWellBlocks(map<string, vector<IntersectedCell>> &well_indices,
                         vector<WellDefinition> &wells, int rank = 0);

};

}
}

#endif //FIELDOPT_WICALC_RIXX_H
