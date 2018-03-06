//
// Created by bellout on 3/5/18.
//

#ifndef FIELDOPT_GEOMETRY_TOOLS_H
#define FIELDOPT_GEOMETRY_TOOLS_H

// RESINSIGHT: FWK/VIZFWK/LIBCORE\LIBGEOMETRY ----------------------
#include "resinxx/rixx_core_geom/cvfBase.h"

// FieldOpt::RESINXX -----------------------------------------------
//
#include "resinxx/well_path.h"

// -----------------------------------------------------------------
namespace cvf {

using std::vector;
using std::array;
using std::set;
using std::list;
using std::pair;

// =================================================================
// Internal class for intersection point info
struct HexIntersectionInfo {

 public:
  HexIntersectionInfo(cvf::Vec3d intersectionPoint,
                      bool isIntersectionEntering,
                      cvf::StructGridInterface::FaceType face,
                      size_t hexIndex)
      : m_intersectionPoint(intersectionPoint),
        m_isIntersectionEntering(isIntersectionEntering),
        m_face(face),
        m_hexIndex(hexIndex) {}

  cvf::Vec3d m_intersectionPoint;
  bool m_isIntersectionEntering;
  cvf::StructGridInterface::FaceType m_face;
  size_t m_hexIndex;
};

// =================================================================
bool operator<( const HexIntersectionInfo& hi1,
                const HexIntersectionInfo& hi2);

// =================================================================
// SAVE FOR LATER
// Specialized Line - Hex intersection
struct RigHexIntersectionTools {

  static int lineHexCellIntersection(const cvf::Vec3d p1,
                                     const cvf::Vec3d p2,
                                     const cvf::Vec3d hexCorners[8],
                                     const size_t hexIndex,
                                     vector<cvf::HexIntersectionInfo> *intersections);

  static bool isPointInCell(const cvf::Vec3d point,
                            const cvf::Vec3d hexCorners[8]);

  static bool planeHexCellIntersection(
      cvf::Vec3d *hexCorners,
      cvf::Plane fracturePlane,
      list<pair<cvf::Vec3d, cvf::Vec3d > > &intersectionLineSegments);

  static bool planeHexIntersectionPolygons(array<cvf::Vec3d, 8> hexCorners,
                                           cvf::Mat4d transformMatrixForPlane,
                                           vector<vector<cvf::Vec3d > > &polygons);
};

}

//}
//}
#endif //FIELDOPT_GEOMETRY_TOOLS_H
