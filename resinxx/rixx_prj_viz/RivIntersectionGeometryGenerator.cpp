/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) Statoil ASA
//  Copyright (C) Ceetron Solutions AS
//
//  ResInsight is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ResInsight is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.
//
//  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
//  for more details.
//
/////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
#include "RivIntersectionGeometryGenerator.h"

// ---------------------------------------------------------------
// #include "RigMainGrid.h"
#include "../rixx_grid/rigrid.h"
// #include "RigResultAccessor.h"

// ---------------------------------------------------------------
#include "RimIntersection.h"
// #include "Rim3dView.h"

// ---------------------------------------------------------------
#include "RivHexGridIntersectionTools.h"
#include "RivIntersectionPartMgr.h"

// ---------------------------------------------------------------
#include "../rixx_app_fwk/cafHexGridIntersectionTools.h"
//#include "cafDisplayCoordTransform.h"

// ---------------------------------------------------------------
//#include "cvfDrawableGeo.h"
#include "cvfGeometryTools.h"
#include "cvfPlane.h"
//#include "cvfPrimitiveSetDirect.h"
//#include "cvfPrimitiveSetIndexedUInt.h"
//#include "cvfScalarMapper.h"
#include "cvfRay.h"

// ---------------------------------------------------------------
#include "RivSectionFlattner.h"


// ===============================================================
//cvf::ref<caf::DisplayCoordTransform>
//displayCoordTransform(const RimIntersection* intersection) {
//
//  Rim3dView* rimView = nullptr;
//  intersection->firstAncestorOrThisOfType(rimView);
//  CVF_ASSERT(rimView);
//
//  cvf::ref<caf::DisplayCoordTransform> transForm = rimView->displayCoordTransform();
//  return transForm;
//}

// ===============================================================
// isFlattened means to transform each flat section of the
// intersection onto the XZ plane placed adjacent to each
// other as if they were rotated around the common extrusion
// line like a hinge
RivIntersectionGeometryGenerator::
RivIntersectionGeometryGenerator(RimIntersection* crossSection,
                                 std::vector<std::vector<cvf::Vec3d> > &polylines,
                                 const cvf::Vec3d& extrusionDirection,
                                 const RivIntersectionHexGridInterface* grid,
                                 bool isFlattened,
                                 const cvf::Vec3d& flattenedPolylineStartPoint)
    : m_crossSection(crossSection),
      m_polyLines(polylines),
      m_extrusionDirection(extrusionDirection),
      m_hexGrid(grid),
      m_isFlattened(isFlattened),
      m_flattenedPolylineStartPoint(flattenedPolylineStartPoint)
{
  m_triangleVxes = new cvf::Vec3fArray;
  m_cellBorderLineVxes = new cvf::Vec3fArray;
  m_faultCellBorderLineVxes = new cvf::Vec3fArray;

  if (m_isFlattened) m_extrusionDirection = -cvf::Vec3d::Z_AXIS;
}

// ===============================================================
RivIntersectionGeometryGenerator::~RivIntersectionGeometryGenerator() {
}

// ===============================================================
void RivIntersectionGeometryGenerator::
calculateSegementTransformPrLinePoint() {

  if ( m_isFlattened ) {

    // -----------------------------------------------------------
    if ( !(m_polyLines.size() && m_polyLines.back().size()) ) { return; }

    cvf::Vec3d startOffset = m_flattenedPolylineStartPoint;

    for ( size_t pLineIdx = 0; pLineIdx < m_polyLines.size(); ++pLineIdx ) {

      // ---------------------------------------------------------
      const std::vector<cvf::Vec3d>& polyLine = m_polyLines[pLineIdx];
      startOffset.z() = polyLine[0].z();

      // ---------------------------------------------------------
      m_segementTransformPrLinePoint.
          emplace_back(RivSectionFlattner::
                       calculateFlatteningCSsForPolyline(polyLine,
                                                         m_extrusionDirection,
                                                         startOffset,
                                                         &startOffset));
    }

  } else {

    // -----------------------------------------------------------
    m_segementTransformPrLinePoint.clear();

    cvf::Mat4d invSectionCS = cvf::Mat4d::fromTranslation(-m_hexGrid->displayOffset());

    for ( const auto & polyLine : m_polyLines ) {

      // ---------------------------------------------------------
      m_segementTransformPrLinePoint.emplace_back();
      std::vector<cvf::Mat4d>& segmentTransforms =
          m_segementTransformPrLinePoint.back();

      // ---------------------------------------------------------
      for ( size_t lIdx = 0; lIdx < polyLine.size(); ++lIdx ) {
        segmentTransforms.push_back(invSectionCS);
      }
    }
  }
}

// ===============================================================
void RivIntersectionGeometryGenerator::calculateFlattenedOrOffsetedPolyline() {

  CVF_ASSERT(m_segementTransformPrLinePoint.size() == m_polyLines.size());

  // -------------------------------------------------------------
  for ( size_t pLineIdx = 0; pLineIdx < m_polyLines.size(); ++pLineIdx ) {

    // -----------------------------------------------------------
    m_flattenedOrOffsettedPolyLines.emplace_back();
    const std::vector<cvf::Vec3d>& polyLine = m_polyLines[pLineIdx];

    // -----------------------------------------------------------
    CVF_ASSERT(polyLine.size() == m_polyLines[pLineIdx].size());

    // -----------------------------------------------------------
    for ( size_t pIdx = 0; pIdx < polyLine.size(); ++pIdx ) {
      m_flattenedOrOffsettedPolyLines.back().push_back(
          polyLine[pIdx].getTransformedPoint(
              m_segementTransformPrLinePoint[pLineIdx][pIdx]));
    }
  }
  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__);
}

// ===============================================================
void RivIntersectionGeometryGenerator::calculateArrays() {

  // -------------------------------------------------------------
  if (m_triangleVxes->size()) return;

  // -------------------------------------------------------------
  m_extrusionDirection.normalize();
  std::vector<cvf::Vec3f> triangleVertices;
  std::vector<cvf::Vec3f> cellBorderLineVxes;
  std::vector<cvf::Vec3f> faultCellBorderLineVxes;
  cvf::Vec3d displayOffset = m_hexGrid->displayOffset();
  cvf::BoundingBox gridBBox = m_hexGrid->boundingBox();

  // -------------------------------------------------------------
  calculateSegementTransformPrLinePoint();
  calculateFlattenedOrOffsetedPolyline();

  // -------------------------------------------------------------
  for (size_t pLineIdx = 0; pLineIdx < m_polyLines.size(); ++pLineIdx) {

    // -----------------------------------------------------------
    const std::vector<cvf::Vec3d>& polyLine = m_polyLines[pLineIdx];

    if (polyLine.size() < 2) continue;

    // -----------------------------------------------------------
    size_t lineCount = polyLine.size();
    size_t lIdx = 0;
    while ( lIdx < lineCount - 1) {

      // ---------------------------------------------------------
      size_t idxToNextP =
          RivSectionFlattner::indexToNextValidPoint(
              polyLine, m_extrusionDirection, lIdx);

      if (idxToNextP == size_t(-1)) break;

      // ---------------------------------------------------------
      cvf::Vec3d p1 = polyLine[lIdx];
      cvf::Vec3d p2 = polyLine[idxToNextP];

      // ---------------------------------------------------------
      cvf::BoundingBox sectionBBox;
      sectionBBox.add(p1);
      sectionBBox.add(p2);

      // ---------------------------------------------------------
      cvf::Vec3d maxHeightVec;

      // ---------------------------------------------------------
      double maxSectionHeightUp = 0;
      double maxSectionHeightDown = 0;

//      // ---------------------------------------------------------
//      if (m_crossSection->type == RimIntersection::CS_AZIMUTHLINE) {
//
//        // -------------------------------------------------------
//        maxSectionHeightUp = m_crossSection->lengthUp();
//        maxSectionHeightDown = m_crossSection->lengthDown();
//
//        // -------------------------------------------------------
//        if (maxSectionHeightUp + maxSectionHeightDown == 0) {
//          return;
//        }
//
//        // -------------------------------------------------------
//        cvf::Vec3d maxHeightVecDown = m_extrusionDirection*maxSectionHeightUp;
//        cvf::Vec3d maxHeightVecUp = m_extrusionDirection*maxSectionHeightDown;
//
//        // -------------------------------------------------------
//        sectionBBox.add(p1 + maxHeightVecUp);
//        sectionBBox.add(p1 - maxHeightVecDown);
//        sectionBBox.add(p2 + maxHeightVecUp);
//        sectionBBox.add(p2 - maxHeightVecDown);
//
//        // -------------------------------------------------------
//        maxHeightVec = maxHeightVecUp + maxHeightVecDown;
//
//      } else {

        // -------------------------------------------------------
        maxHeightVec = m_extrusionDirection*gridBBox.radius();

        sectionBBox.add(p1 + maxHeightVec);
        sectionBBox.add(p1 - maxHeightVec);
        sectionBBox.add(p2 + maxHeightVec);
        sectionBBox.add(p2 - maxHeightVec);
//      }

      // ---------------------------------------------------------
      std::vector<size_t> columnCellCandidates;
      m_hexGrid->findIntersectingCells(sectionBBox,
                                       &columnCellCandidates);

      // ---------------------------------------------------------
      cvf::Plane plane;
      plane.setFromPoints(p1,
                          p2,
                          p2 + maxHeightVec);

      // ---------------------------------------------------------
      cvf::Plane p1Plane;
      p1Plane.setFromPoints(p1,
                            p1 + maxHeightVec,
                            p1 + plane.normal());

      // ---------------------------------------------------------
      cvf::Plane p2Plane;
      p2Plane.setFromPoints(p2,
                            p2 + maxHeightVec,
                            p2 - plane.normal());

      // ---------------------------------------------------------
      std::vector<caf::HexGridIntersectionTools::ClipVx>
          hexPlaneCutTriangleVxes;
      hexPlaneCutTriangleVxes.reserve(5*3);

      // ---------------------------------------------------------
      std::vector<int> cellFaceForEachTriangleEdge;
      cellFaceForEachTriangleEdge.reserve(5*3);

      // ---------------------------------------------------------
      cvf::Vec3d cellCorners[8];
      size_t cornerIndices[8];

      // ---------------------------------------------------------
      cvf::Mat4d invSectionCS =
          m_segementTransformPrLinePoint[pLineIdx][lIdx];

      // ---------------------------------------------------------
      for (size_t cccIdx = 0;
           cccIdx < columnCellCandidates.size(); ++cccIdx) {

        size_t globalCellIdx = columnCellCandidates[cccIdx];

        // -------------------------------------------------------
        if (!m_hexGrid->useCell(globalCellIdx)) continue;

        // -------------------------------------------------------
        hexPlaneCutTriangleVxes.clear();
        m_hexGrid->cellCornerVertices(globalCellIdx, cellCorners);
        m_hexGrid->cellCornerIndices(globalCellIdx, cornerIndices);

        // -------------------------------------------------------
        caf::HexGridIntersectionTools::
        planeHexIntersectionMC(plane,
                               cellCorners,
                               cornerIndices,
                               &hexPlaneCutTriangleVxes,
                               &cellFaceForEachTriangleEdge);

        // -------------------------------------------------------
//        if (m_crossSection->type == RimIntersection::CS_AZIMUTHLINE) {
//
//          // -----------------------------------------------------
//          bool hasAnyPointsOnSurface = false;
//
//          // -----------------------------------------------------
//          for (caf::HexGridIntersectionTools::ClipVx
//                vertex : hexPlaneCutTriangleVxes) {
//
//            // ---------------------------------------------------
//            cvf::Vec3d temp = vertex.vx - p1;
//            double dot = temp.dot(m_extrusionDirection);
//            double lengthCheck = 0;
//
//            // ---------------------------------------------------
//            if (dot < 0) {
//              lengthCheck = maxSectionHeightUp;
//
//            } else {
//              lengthCheck = maxSectionHeightDown;
//            }
//
//            // ---------------------------------------------------
//            double distance =
//                cvf::Math::sqrt(
//                    cvf::GeometryTools::linePointSquareDist(p1,
//                                                            p2,
//                                                            vertex.vx));
//
//            // ---------------------------------------------------
//            if (distance < lengthCheck) {
//              hasAnyPointsOnSurface = true;
//              break;
//            }
//          }
//
//          // -----------------------------------------------------
//          if (!hasAnyPointsOnSurface) {
//            continue;
//          }
//        }

        // -------------------------------------------------------
        std::vector<caf::HexGridIntersectionTools::ClipVx>
            clippedTriangleVxes;
        std::vector<int> cellFaceForEachClippedTriangleEdge;

        caf::HexGridIntersectionTools::
        clipTrianglesBetweenTwoParallelPlanes(
            hexPlaneCutTriangleVxes,
            cellFaceForEachTriangleEdge,
            p1Plane,
            p2Plane,
            &clippedTriangleVxes,
            &cellFaceForEachClippedTriangleEdge);

        // -------------------------------------------------------
        size_t clippedTriangleCount = clippedTriangleVxes.size()/3;

        // -------------------------------------------------------
        for (uint tIdx = 0; tIdx < clippedTriangleCount; ++tIdx) {

          uint triVxIdx = tIdx*3;

          // -----------------------------------------------------
          // Accumulate triangle vertices
          cvf::Vec3d p0(clippedTriangleVxes[triVxIdx+0].vx);
          cvf::Vec3d p1(clippedTriangleVxes[triVxIdx+1].vx);
          cvf::Vec3d p2(clippedTriangleVxes[triVxIdx+2].vx);

          // -----------------------------------------------------
          p0 = p0.getTransformedPoint(invSectionCS);
          p1 = p1.getTransformedPoint(invSectionCS);
          p2 = p2.getTransformedPoint(invSectionCS);

          // -----------------------------------------------------
          triangleVertices.emplace_back(p0);
          triangleVertices.emplace_back(p1);
          triangleVertices.emplace_back(p2);

          // -----------------------------------------------------
          // Accumulate mesh lines
#define isFace( faceEnum ) (0 <= faceEnum && faceEnum <= 5 )

          using FaceType = cvf::StructGridInterface::FaceType;

          // -----------------------------------------------------
          if ( isFace(cellFaceForEachClippedTriangleEdge[triVxIdx]) ) {

            // ---------------------------------------------------
            if ( m_hexGrid->findFaultFromCellIndexAndCellFace(
                globalCellIdx,
                (FaceType)cellFaceForEachClippedTriangleEdge[triVxIdx]) ) {

              // -------------------------------------------------
              faultCellBorderLineVxes.emplace_back(p0);
              faultCellBorderLineVxes.emplace_back(p1);

            } else {

              // -------------------------------------------------
              cellBorderLineVxes.emplace_back(p0);
              cellBorderLineVxes.emplace_back(p1);
            }
          }

          // -----------------------------------------------------
          if ( isFace(cellFaceForEachClippedTriangleEdge[triVxIdx+1]) ) {

            // ---------------------------------------------------
            if ( m_hexGrid->findFaultFromCellIndexAndCellFace(
                globalCellIdx,
                (FaceType)cellFaceForEachClippedTriangleEdge[triVxIdx+1]) ) {

              // -------------------------------------------------
              faultCellBorderLineVxes.emplace_back(p1);
              faultCellBorderLineVxes.emplace_back(p2);

            } else {

              // -------------------------------------------------
              cellBorderLineVxes.emplace_back(p1);
              cellBorderLineVxes.emplace_back(p2);
            }
          }

          // -----------------------------------------------------
          if ( isFace(cellFaceForEachClippedTriangleEdge[triVxIdx+2]) ) {

            // ---------------------------------------------------
            if ( m_hexGrid->findFaultFromCellIndexAndCellFace(
                globalCellIdx,
                (FaceType)cellFaceForEachClippedTriangleEdge[triVxIdx+2]) ) {

              // -------------------------------------------------
              faultCellBorderLineVxes.emplace_back(p2);
              faultCellBorderLineVxes.emplace_back(p0);

            } else {

              // -------------------------------------------------
              cellBorderLineVxes.emplace_back(p2);
              cellBorderLineVxes.emplace_back(p0);
            }
          }

          // -----------------------------------------------------
          // Mapping to cell index
          m_triangleToCellIdxMap.push_back(globalCellIdx);

          // -----------------------------------------------------
          // Interpolation from nodes
          for (int i = 0; i < 3; ++i) {

            caf::HexGridIntersectionTools::ClipVx
                cvx = clippedTriangleVxes[triVxIdx + i];

            // ---------------------------------------------------
            if (cvx.isVxIdsNative) {
              m_triVxToCellCornerWeights.push_back(
                  RivIntersectionVertexWeights(cvx.clippedEdgeVx1Id,
                                               cvx.clippedEdgeVx2Id,
                                               cvx.normDistFromEdgeVx1));

            } else {

              // -------------------------------------------------
              caf::HexGridIntersectionTools::ClipVx
                  cvx1 = hexPlaneCutTriangleVxes[cvx.clippedEdgeVx1Id];

              caf::HexGridIntersectionTools::ClipVx
                  cvx2 = hexPlaneCutTriangleVxes[cvx.clippedEdgeVx2Id];

              // -------------------------------------------------
              m_triVxToCellCornerWeights.push_back(
                  RivIntersectionVertexWeights(
                      cvx1.clippedEdgeVx1Id, cvx1.clippedEdgeVx2Id, cvx1.normDistFromEdgeVx1,
                      cvx2.clippedEdgeVx1Id, cvx2.clippedEdgeVx2Id, cvx2.normDistFromEdgeVx1,
                      cvx.normDistFromEdgeVx1));

            }
          }
        }
      }

      // ---------------------------------------------------------
      lIdx = idxToNextP;
    }
  }
  m_triangleVxes->assign(triangleVertices);
  m_cellBorderLineVxes->assign(cellBorderLineVxes);
  m_faultCellBorderLineVxes->assign(faultCellBorderLineVxes);

  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__, "");
}


// ===============================================================
// Generate surface drawable geo from the specified region
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::generateSurface() {
//
//  // -------------------------------------------------------------
//  calculateArrays();
//
//  // -------------------------------------------------------------
//  CVF_ASSERT(m_triangleVxes.notNull());
//  if (m_triangleVxes->size() == 0) return nullptr;
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::DrawableGeo> geo = new cvf::DrawableGeo;
//  geo->setFromTriangleVertexArray(m_triangleVxes.p());
//
//  return geo;
//}


// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createMeshDrawable() {
//
//  // -------------------------------------------------------------
//  if (!(m_cellBorderLineVxes.notNull()
//      && m_cellBorderLineVxes->size() != 0)) return nullptr;
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::DrawableGeo> geo = new cvf::DrawableGeo;
//  geo->setVertexArray(m_cellBorderLineVxes.p());
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::PrimitiveSetDirect>
//      prim = new cvf::PrimitiveSetDirect(cvf::PT_LINES);
//
//  prim->setIndexCount(m_cellBorderLineVxes->size());
//
//  // -------------------------------------------------------------
//  geo->addPrimitiveSet(prim.p());
//  return geo;
//}

// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createFaultMeshDrawable() {
//
//  if (!(m_faultCellBorderLineVxes.notNull()
//      && m_faultCellBorderLineVxes->size() != 0)) return nullptr;
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::DrawableGeo> geo = new cvf::DrawableGeo;
//  geo->setVertexArray(m_faultCellBorderLineVxes.p());
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::PrimitiveSetDirect> prim =
//      new cvf::PrimitiveSetDirect(cvf::PT_LINES);
//  prim->setIndexCount(m_faultCellBorderLineVxes->size());
//
//  // -------------------------------------------------------------
//  geo->addPrimitiveSet(prim.p());
//  return geo;
//}


// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createLineAlongPolylineDrawable() {
//  print_ri_hck(__func__, __FILE__,
//                       "RETURN: createLineAlongPolylineDrawable"
//                           "(m_flattenedOrOffsettedPolyLines)");
//  return createLineAlongPolylineDrawable(m_flattenedOrOffsettedPolyLines);
//}


// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createLineAlongExtrusionLineDrawable(
//    const std::vector<cvf::Vec3d>& extrusionLine) {
//
//  // -------------------------------------------------------------
//  cvf::ref<caf::DisplayCoordTransform> transform =
//      displayCoordTransform(crossSection());
//
//  std::vector<cvf::Vec3d> displayCoords;
//
//  // -------------------------------------------------------------
//  for (const auto& pt : extrusionLine) {
//    displayCoords.push_back(transform->translateToDisplayCoord(pt));
//  }
//
//  // -------------------------------------------------------------
//  print_ri_hck(__func__, __FILE__,
//                       "RETURN: createLineAlongPolylineDrawable("
//                           "vector<vector<cvf::Vec3d>>({ displayCoords })");
//  return createLineAlongPolylineDrawable(
//      std::vector<std::vector<cvf::Vec3d>>({ displayCoords }));
//}


// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createLineAlongPolylineDrawable(
//    const std::vector<std::vector<cvf::Vec3d> >& polyLines) {
//
//  // -------------------------------------------------------------
//  std::vector<cvf::uint> lineIndices;
//  std::vector<cvf::Vec3f> vertices;
//
//  for (size_t pLineIdx = 0; pLineIdx < polyLines.size(); ++pLineIdx) {
//
//    // -----------------------------------------------------------
//    const std::vector<cvf::Vec3d>& polyLine = polyLines[pLineIdx];
//    if (polyLine.size() < 2) continue;
//
//    // -----------------------------------------------------------
//    for (size_t i = 0; i < polyLine.size(); ++i) {
//
//      // ---------------------------------------------------------
//      vertices.push_back(cvf::Vec3f(polyLine[i]));
//      if (i < polyLine.size() - 1) {
//        lineIndices.push_back(static_cast<cvf::uint>(i));
//        lineIndices.push_back(static_cast<cvf::uint>(i + 1));
//      }
//    }
//  }
//
//  // -------------------------------------------------------------
//  if (vertices.size() == 0) return nullptr;
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::Vec3fArray> vx = new cvf::Vec3fArray;
//  vx->assign(vertices);
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::UIntArray> idxes = new cvf::UIntArray;
//  idxes->assign(lineIndices);
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::PrimitiveSetIndexedUInt>
//      prim = new cvf::PrimitiveSetIndexedUInt(cvf::PT_LINES);
//  prim->setIndices(idxes.p());
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::DrawableGeo> polylineGeo = new cvf::DrawableGeo;
//  polylineGeo->setVertexArray(vx.p());
//  polylineGeo->addPrimitiveSet(prim.p());
//
//  // -------------------------------------------------------------
//  print_ri_hck(__func__, __FILE__, "RETURN polylineGeo");
//
//  // -------------------------------------------------------------
//  return polylineGeo;
//}

// ---------------------------------------------------------------
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createPointsFromPolylineDrawable() {
//  // -------------------------------------------------------------
//  // RETURN FUNCTION => print_ri_hck(__func__, __FILE__);
//  return createPointsFromPolylineDrawable(
//      m_flattenedOrOffsettedPolyLines);
//}

// ---------------------------------------------------------------
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createPointsFromExtrusionLineDrawable(
//    const std::vector<cvf::Vec3d>& extrusionLine) {
//
//  // -------------------------------------------------------------
//  cvf::ref<caf::DisplayCoordTransform>
//      transform = displayCoordTransform(crossSection());
//
//  // -------------------------------------------------------------
//  std::vector<cvf::Vec3d> displayCoords;
//
//  // -------------------------------------------------------------
//  for (const auto& pt : extrusionLine) {
//    displayCoords.push_back(transform->translateToDisplayCoord(pt));
//  }
//
//  // -------------------------------------------------------------
//  return createPointsFromPolylineDrawable(
//      std::vector<std::vector<cvf::Vec3d>>({displayCoords}));
//}

// ===============================================================
//cvf::ref<cvf::DrawableGeo>
//RivIntersectionGeometryGenerator::createPointsFromPolylineDrawable(
//    const std::vector<std::vector<cvf::Vec3d> >& polyLines) {
//
//  std::vector<cvf::Vec3f> vertices;
//
//  // -------------------------------------------------------------
//  for (size_t pLineIdx = 0; pLineIdx < polyLines.size(); ++pLineIdx) {
//
//    const std::vector<cvf::Vec3d>& polyLine = polyLines[pLineIdx];
//
//    // -----------------------------------------------------------
//    for (size_t i = 0; i < polyLine.size(); ++i) {
//      vertices.push_back(cvf::Vec3f(polyLine[i]));
//    }
//  }
//
//  // -------------------------------------------------------------
//  if (vertices.size() == 0) return nullptr;
//
//  cvf::ref<cvf::PrimitiveSetDirect>
//      primSet = new cvf::PrimitiveSetDirect(cvf::PT_POINTS);
//
//  primSet->setStartIndex(0);
//  primSet->setIndexCount(vertices.size());
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::DrawableGeo> geo = new cvf::DrawableGeo;
//
//  // -------------------------------------------------------------
//  cvf::ref<cvf::Vec3fArray> vx = new cvf::Vec3fArray(vertices);
//
//  geo->setVertexArray(vx.p());
//  geo->addPrimitiveSet(primSet.p());
//
//  // -------------------------------------------------------------
//  print_ri_hck(__func__, __FILE__);
//  return geo;
//}

// ===============================================================
const std::vector<size_t>&
RivIntersectionGeometryGenerator::triangleToCellIndex() const {

  CVF_ASSERT(m_triangleVxes->size());
  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__);
  return m_triangleToCellIdxMap;
}

// ===============================================================
const std::vector<RivIntersectionVertexWeights>&
RivIntersectionGeometryGenerator::
triangleVxToCellCornerInterpolationWeights() const {

  CVF_ASSERT(m_triangleVxes->size());
  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__);
  return m_triVxToCellCornerWeights;
}

// ===============================================================
const cvf::Vec3fArray*
RivIntersectionGeometryGenerator::triangleVxes() const {
  CVF_ASSERT(m_triangleVxes->size());
  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__);
  return m_triangleVxes.p();
}

// ===============================================================
RimIntersection*
RivIntersectionGeometryGenerator::crossSection() const {
  // -----------------------------------------------------------
  print_ri_hck(__func__, __FILE__);
  return m_crossSection;
}

// ===============================================================
cvf::Mat4d
RivIntersectionGeometryGenerator::
unflattenTransformMatrix(const cvf::Vec3d& intersectionPointFlat) {

  // -------------------------------------------------------------
  cvf::Mat4d flattenMx =  cvf::Mat4d::IDENTITY;

  for ( size_t pLineIdx = 0;
        pLineIdx < m_flattenedOrOffsettedPolyLines.size(); pLineIdx++ ) {

    // -----------------------------------------------------------
    const std::vector<cvf::Vec3d>& polyLine =
        m_flattenedOrOffsettedPolyLines[pLineIdx];

    // -----------------------------------------------------------
    for(size_t pIdx = 0; pIdx < polyLine.size(); pIdx++) {

      if (polyLine[pIdx].x() >= intersectionPointFlat.x() ) {

        // -------------------------------------------------------
        size_t csIdx = pIdx > 0 ? pIdx - 1: 0;
        flattenMx = m_segementTransformPrLinePoint[pLineIdx][csIdx];
        break;

      } else if (pIdx == polyLine.size() - 1) {
        flattenMx = m_segementTransformPrLinePoint[pLineIdx][pIdx];
      }
    }
  }

  // -------------------------------------------------------------
  return flattenMx.getInverted();
}

// ===============================================================
bool
RivIntersectionGeometryGenerator::isAnyGeometryPresent() const {

  // -------------------------------------------------------------
  print_ri_hck(__func__, __FILE__);

  // -------------------------------------------------------------
  if (m_triangleVxes->size() == 0) {
    return false;

  } else {
    return true;
  }

}

