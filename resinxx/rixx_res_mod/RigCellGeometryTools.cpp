////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2017     Statoil ASA
//
// ResInsight is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// ResInsight is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// See the GNU General Public License at
// <http://www.gnu.org/licenses/gpl.html> for more details.
//
////////////////////////////////////////////////////////////////////
//
// Modified by M.Bellout on 3/5/18.
//

// RESINSIGHT: FWK/APPFWK/COMMONCODE\VIZEXT\PROJDATAMOD ------------
#include "../rixx_app_fwk/cvfStructGrid.h"
#include "../rixx_app_fwk/cafHexGridIntersectionTools.h"

// RESINSIGHT: APPLICATIONCODE/RESERVOIRDATAMODEL ------------------
#include "RigCellGeometryTools.h"
#include "cvfGeometryTools.h"

//#ifdef USE_PROTOTYPE_FEATURE_FRACTURES
//#include "clipper/clipper.hpp"
//#endif // USE_PROTOTYPE_FEATURE_FRACTURES

// ╦═╗  ╦  ╔═╗  ╔═╗  ╔═╗  ╦    ╦    ╔═╗  ╔═╗  ╔═╗  ╔╦╗  ╔═╗  ╔╦╗  ╦═╗  ╦ ╦
// ╠╦╝  ║  ║ ╦  ║    ║╣   ║    ║    ║ ╦  ║╣   ║ ║  ║║║  ║╣    ║   ╠╦╝  ╚╦╝
// ╩╚═  ╩  ╚═╝  ╚═╝  ╚═╝  ╩═╝  ╩═╝  ╚═╝  ╚═╝  ╚═╝  ╩ ╩  ╚═╝   ╩   ╩╚═   ╩
// =================================================================
void
RigCellGeometryTools::createPolygonFromLineSegments(
    list<pair<cvf::Vec3d, cvf::Vec3d>> &intersectionLineSegments,
    vector<vector<cvf::Vec3d>> &polygons) {

  bool startNewPolygon = true;
  while (!intersectionLineSegments.empty()) {

    if (startNewPolygon) {

      vector<cvf::Vec3d> polygon;

      //Add first line segments to polygon and remove from list
      pair<cvf::Vec3d, cvf::Vec3d >
          linesegment = intersectionLineSegments.front();

      polygon.push_back(linesegment.first);
      polygon.push_back(linesegment.second);
      intersectionLineSegments.remove(linesegment);
      polygons.push_back(polygon);
      startNewPolygon = false;
    }

    vector<cvf::Vec3d>& polygon = polygons.back();

    //Search remaining list for next point...
    bool isFound = false;
    float tolerance = 0.0001f;

    for (list<pair<cvf::Vec3d, cvf::Vec3d > >::iterator
             lIt = intersectionLineSegments.begin();
         lIt != intersectionLineSegments.end(); lIt++) {

      cvf::Vec3d lineSegmentStart = lIt->first;
      cvf::Vec3d lineSegmentEnd = lIt->second;
      cvf::Vec3d polygonEnd = polygon.back();

      double lineSegmentLength =
          (lineSegmentStart - lineSegmentEnd).lengthSquared();

      if (lineSegmentLength < tolerance*tolerance) {
        intersectionLineSegments.erase(lIt);
        isFound = true;
        break;
      }

      double lineSegmentStartDiff =
          (lineSegmentStart - polygonEnd).lengthSquared();

      if (lineSegmentStartDiff < tolerance*tolerance) {
        polygon.push_back(lIt->second);
        intersectionLineSegments.erase(lIt);
        isFound = true;
        break;
      }

      double lineSegmentEndDiff =
          (lineSegmentEnd - polygonEnd).lengthSquared();

      if (lineSegmentEndDiff < tolerance*tolerance) {
        polygon.push_back(lIt->first);
        intersectionLineSegments.erase(lIt);
        isFound = true;
        break;
      }
    }

    if (isFound) {
      continue;
    } else {
      startNewPolygon = true;
    }
  }
}

// -----------------------------------------------------------------
void
RigCellGeometryTools::findCellLocalXYZ(const array<cvf::Vec3d, 8>& hexCorners,
                                       cvf::Vec3d& localXdirection,
                                       cvf::Vec3d& localYdirection,
                                       cvf::Vec3d& localZdirection) {

  cvf::ubyte faceVertexIndices[4];
  cvf::StructGridInterface::FaceEnum face;

  face = cvf::StructGridInterface::NEG_I;
  cvf::StructGridInterface::cellFaceVertexIndices(face, faceVertexIndices);

  cvf::Vec3d faceCenterNegI =
      cvf::GeometryTools::computeFaceCenter(hexCorners[faceVertexIndices[0]],
                                            hexCorners[faceVertexIndices[1]],
                                            hexCorners[faceVertexIndices[2]],
                                            hexCorners[faceVertexIndices[3]]);
  //TODO: Should we use face centroids instead of face centers?

  face = cvf::StructGridInterface::POS_I;
  cvf::StructGridInterface::cellFaceVertexIndices(face, faceVertexIndices);

  cvf::Vec3d faceCenterPosI =
      cvf::GeometryTools::computeFaceCenter(hexCorners[faceVertexIndices[0]],
                                            hexCorners[faceVertexIndices[1]],
                                            hexCorners[faceVertexIndices[2]],
                                            hexCorners[faceVertexIndices[3]]);

  face = cvf::StructGridInterface::NEG_J;
  cvf::StructGridInterface::cellFaceVertexIndices(face, faceVertexIndices);

  cvf::Vec3d faceCenterNegJ =
      cvf::GeometryTools::computeFaceCenter(hexCorners[faceVertexIndices[0]],
                                            hexCorners[faceVertexIndices[1]],
                                            hexCorners[faceVertexIndices[2]],
                                            hexCorners[faceVertexIndices[3]]);

  face = cvf::StructGridInterface::POS_J;
  cvf::StructGridInterface::cellFaceVertexIndices(face, faceVertexIndices);

  cvf::Vec3d faceCenterPosJ =
      cvf::GeometryTools::computeFaceCenter(hexCorners[faceVertexIndices[0]],
                                            hexCorners[faceVertexIndices[1]],
                                            hexCorners[faceVertexIndices[2]],
                                            hexCorners[faceVertexIndices[3]]);

  cvf::Vec3d faceCenterCenterVectorI = faceCenterPosI - faceCenterNegI;
  cvf::Vec3d faceCenterCenterVectorJ = faceCenterPosJ - faceCenterNegJ;

  localZdirection.cross(faceCenterCenterVectorI,
                        faceCenterCenterVectorJ);
  localZdirection.normalize();

  cvf::Vec3d crossPoductJZ;
  crossPoductJZ.cross(faceCenterCenterVectorJ, localZdirection);

  localXdirection = faceCenterCenterVectorI + crossPoductJZ;
  localXdirection.normalize();

  cvf::Vec3d crossPoductIZ;
  crossPoductIZ.cross(faceCenterCenterVectorI, localZdirection);

  localYdirection = faceCenterCenterVectorJ - crossPoductIZ;
  localYdirection.normalize();

  //TODO: Check if we end up with 0-vectors, and handle this case...
}

// -----------------------------------------------------------------
double
RigCellGeometryTools::polygonLengthInLocalXdirWeightedByArea(
    vector<cvf::Vec3d> polygonToCalcLengthOf) {

  //Find bounding box
  cvf::BoundingBox polygonBBox;
  for (cvf::Vec3d nodeCoord : polygonToCalcLengthOf) polygonBBox.add(nodeCoord);
  cvf::Vec3d bboxCorners[8];
  polygonBBox.cornerVertices(bboxCorners);

  //Split bounding box in multiple polygons (2D)
  int resolutionOfLengthCalc = 20;
  double widthOfPolygon = polygonBBox.extent().y() / resolutionOfLengthCalc;

  vector<double> areasOfPolygonContributions;
  vector<double> lengthOfPolygonContributions;

  cvf::Vec3d directionOfLength(1, 0, 0);

  for (int i = 0; i < resolutionOfLengthCalc; i++) {

    cvf::Vec3d pointOnLine1(bboxCorners[0].x(),
                            bboxCorners[0].y() + i*widthOfPolygon, 0);

    cvf::Vec3d pointOnLine2(bboxCorners[0].x(),
                            bboxCorners[0].y() + (i + 1)*widthOfPolygon, 0);

    pair<cvf::Vec3d, cvf::Vec3d> line1 =
        getLineThroughBoundingBox(directionOfLength,
                                  polygonBBox,
                                  pointOnLine1);

    pair<cvf::Vec3d, cvf::Vec3d> line2 =
        getLineThroughBoundingBox(directionOfLength,
                                  polygonBBox,
                                  pointOnLine2);

    vector<cvf::Vec3d> polygon;
    polygon.push_back(line1.first);
    polygon.push_back(line1.second);
    polygon.push_back(line2.second);
    polygon.push_back(line2.first);

    //Use clipper to find overlap between bbpolygon and fracture
    vector<vector<cvf::Vec3d> > clippedPolygons =
        intersectPolygons(polygonToCalcLengthOf, polygon);

    double area = 0;
    double length = 0;
    cvf::Vec3d areaVector = cvf::Vec3d::ZERO;

    //Calculate length (max-min) and area
    for (vector<cvf::Vec3d> clippedPolygon : clippedPolygons) {
      areaVector = cvf::GeometryTools::polygonAreaNormal3D(clippedPolygon);
      area += areaVector.length();
      length += (getLengthOfPolygonAlongLine(line1, clippedPolygon)
          + getLengthOfPolygonAlongLine(line2, clippedPolygon)) / 2;
    }
    areasOfPolygonContributions.push_back(area);
    lengthOfPolygonContributions.push_back(length);
  }

  //Calculate area-weighted length average.
  double totalArea = 0.0;
  double totalAreaXlength = 0.0;

  for (size_t i = 0; i < areasOfPolygonContributions.size(); i++) {
    totalArea += areasOfPolygonContributions[i];
    totalAreaXlength += (areasOfPolygonContributions[i] * lengthOfPolygonContributions[i]);
  }

  double areaWeightedLength = totalAreaXlength / totalArea;
  return areaWeightedLength;
}

#ifdef USE_PROTOTYPE_FEATURE_FRACTURES

double clipperConversionFactor = 10000; //For transform to clipper int

ClipperLib::IntPoint toClipperPoint(const cvf::Vec3d& cvfPoint)
{
    int xInt = cvfPoint.x()*clipperConversionFactor;
    int yInt = cvfPoint.y()*clipperConversionFactor;
    int zInt = cvfPoint.z()*clipperConversionFactor;
    return  ClipperLib::IntPoint(xInt, yInt, zInt);
}

cvf::Vec3d fromClipperPoint(const ClipperLib::IntPoint& clipPoint)
{
    double zDValue;

    if (clipPoint.Z == numeric_limits<int>::max())
    {
        zDValue = HUGE_VAL;
    }
    else
    {
        zDValue = clipPoint.Z;
    }

    return cvf::Vec3d (clipPoint.X, clipPoint.Y, zDValue ) /clipperConversionFactor;
}
#endif // USE_PROTOTYPE_FEATURE_FRACTURES

// -----------------------------------------------------------------
vector<vector<cvf::Vec3d> >
RigCellGeometryTools::intersectPolygons(vector<cvf::Vec3d> polygon1,
                                        vector<cvf::Vec3d> polygon2) {

  vector<vector<cvf::Vec3d> > clippedPolygons;

#ifdef USE_PROTOTYPE_FEATURE_FRACTURES
  // Convert to int for clipper library and store as clipper "path"
    ClipperLib::Path polygon1path;
    for (cvf::Vec3d& v : polygon1)
    {
        polygon1path.push_back(toClipperPoint(v));
    }

    ClipperLib::Path polygon2path;
    for (cvf::Vec3d& v : polygon2)
    {
        polygon2path.push_back(toClipperPoint(v));
    }

    ClipperLib::Clipper clpr;
    clpr.AddPath(polygon1path, ClipperLib::ptSubject, true);
    clpr.AddPath(polygon2path, ClipperLib::ptClip, true);

    ClipperLib::Paths solution;
    clpr.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd);

    // Convert back to vector<vector<cvf::Vec3d> >
    for (ClipperLib::Path pathInSol : solution)
    {
        vector<cvf::Vec3d> clippedPolygon;
        for (ClipperLib::IntPoint IntPosition : pathInSol)
        {
            clippedPolygon.push_back(fromClipperPoint(IntPosition));
        }
        clippedPolygons.push_back(clippedPolygon);
    }
#endif // USE_PROTOTYPE_FEATURE_FRACTURES

  return clippedPolygons;
}


// -----------------------------------------------------------------
#ifdef USE_PROTOTYPE_FEATURE_FRACTURES
void fillInterpolatedSubjectZ(ClipperLib::IntPoint& e1bot,
                              ClipperLib::IntPoint& e1top,
                              ClipperLib::IntPoint& e2bot,
                              ClipperLib::IntPoint& e2top,
                              ClipperLib::IntPoint& pt)
{
    ClipperLib::IntPoint ePLbot;
    ClipperLib::IntPoint ePLtop;

    if (e1top.Z == numeric_limits<int>::max())
    {
        ePLtop = e2top;
        ePLbot = e2bot;
    }
    else
    {
        ePLtop = e1top;
        ePLbot = e1bot;
    }

    double ePLXRange = (ePLtop.X - ePLbot.X);
    double ePLYRange = (ePLtop.Y - ePLbot.Y);

    double ePLLength =  sqrt(ePLXRange*ePLXRange + ePLYRange*ePLYRange);

    if (ePLLength <= 1)
    {
        pt.Z = ePLbot.Z;
        return;
    }

    double ePLBotPtXRange = pt.X - ePLbot.X;
    double ePLBotPtYRange = pt.Y - ePLbot.Y;

    double ePLBotPtLength =  sqrt(ePLBotPtXRange*ePLBotPtXRange + ePLBotPtYRange*ePLBotPtYRange);

    double fraction = ePLBotPtLength/ePLLength;

    pt.Z = std::nearbyint( ePLbot.Z + fraction*(ePLtop.Z - ePLbot.Z) );
}

// -----------------------------------------------------------------
///
// -----------------------------------------------------------------
void fillUndefinedZ(ClipperLib::IntPoint& e1bot,
                              ClipperLib::IntPoint& e1top,
                              ClipperLib::IntPoint& e2bot,
                              ClipperLib::IntPoint& e2top,
                              ClipperLib::IntPoint& pt)
{
   pt.Z = std::numeric_limits<int>::max();
}
#endif // USE_PROTOTYPE_FEATURE_FRACTURES

// -----------------------------------------------------------------
// Assumes x.y plane polygon. Polyline might have a Z, and the
// returned Z is the polyline Z, interpolated if it is clipped.
std::vector<std::vector<cvf::Vec3d> >
RigCellGeometryTools::clipPolylineByPolygon(const std::vector<cvf::Vec3d>& polyLine,
                                            const std::vector<cvf::Vec3d>& polygon,
                                            ZInterpolationType interpolType) {

  std::vector<std::vector<cvf::Vec3d> > clippedPolyline;

#ifdef USE_PROTOTYPE_FEATURE_FRACTURES
  //Adjusting polygon to avoid clipper issue with interpolating z-values when lines crosses though polygon vertecies
    std::vector<cvf::Vec3d> adjustedPolygon = ajustPolygonToAvoidIntersectionsAtVertex(polyLine, polygon);

    //Convert to int for clipper library and store as clipper "path"
    ClipperLib::Path polyLinePath;
    for (const cvf::Vec3d& v : polyLine)
    {
        polyLinePath.push_back(toClipperPoint(v));
    }

    ClipperLib::Path polygonPath;
    for (const cvf::Vec3d& v : adjustedPolygon)
    {
        ClipperLib::IntPoint intp = toClipperPoint(v);
        intp.Z = std::numeric_limits<int>::max();
        polygonPath.push_back(intp);
    }

    ClipperLib::Clipper clpr;
    clpr.AddPath(polyLinePath, ClipperLib::ptSubject, false);
    clpr.AddPath(polygonPath, ClipperLib::ptClip, true);

    if ( interpolType == INTERPOLATE_LINE_Z )
    {
        clpr.ZFillFunction(&fillInterpolatedSubjectZ);
    }
    else if ( interpolType == USE_HUGEVAL )
    {
        clpr.ZFillFunction(&fillUndefinedZ);
    }

    ClipperLib::PolyTree solution;
    clpr.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd);

    // We only expect open paths from this method (unless the polyline is self intersecting, a condition we do not handle)
    ClipperLib::Paths solutionPaths;
    ClipperLib::OpenPathsFromPolyTree(solution, solutionPaths);

    //Convert back to std::vector<std::vector<cvf::Vec3d> >
    for (ClipperLib::Path pathInSol : solutionPaths)
    {
        std::vector<cvf::Vec3d> clippedPolygon;
        for (ClipperLib::IntPoint IntPosition : pathInSol)
        {
            clippedPolygon.push_back(fromClipperPoint(IntPosition));
        }
        clippedPolyline.push_back(clippedPolygon);
    }
#endif // USE_PROTOTYPE_FEATURE_FRACTURES

  return clippedPolyline;
}

// -----------------------------------------------------------------
std::pair<cvf::Vec3d, cvf::Vec3d>
RigCellGeometryTools::getLineThroughBoundingBox(cvf::Vec3d lineDirection,
                                                cvf::BoundingBox polygonBBox,
                                                cvf::Vec3d pointOnLine) {

  cvf::Vec3d bboxCorners[8];
  polygonBBox.cornerVertices(bboxCorners);

  cvf::Vec3d startPoint = pointOnLine;
  cvf::Vec3d endPoint = pointOnLine;


  // To avoid doing many iterations in loops
  // below linedirection should be quite large.
  lineDirection.normalize();
  lineDirection =
      lineDirection * polygonBBox.extent().length() / 5;

  // Extend line in positive direction
  while (polygonBBox.contains(startPoint)) {
    startPoint = startPoint + lineDirection;
  }

  // Extend line in negative direction
  while (polygonBBox.contains(endPoint)) {
    endPoint = endPoint - lineDirection;
  }

  std::pair<cvf::Vec3d, cvf::Vec3d> line;
  line = { startPoint, endPoint };
  return line;
}

// -----------------------------------------------------------------
double
RigCellGeometryTools::getLengthOfPolygonAlongLine(
    const std::pair<cvf::Vec3d, cvf::Vec3d>& line,
    const std::vector<cvf::Vec3d>& polygon) {

  cvf::BoundingBox lineBoundingBox;

  for (cvf::Vec3d polygonPoint : polygon) {
    cvf::Vec3d pointOnLine =
        cvf::GeometryTools::projectPointOnLine(line.first,
                                               line.second,
                                               polygonPoint,
                                               nullptr);
    lineBoundingBox.add(pointOnLine);
  }

  double length = lineBoundingBox.extent().length();

  return length;
}

// -----------------------------------------------------------------
std::vector<cvf::Vec3d>
RigCellGeometryTools::ajustPolygonToAvoidIntersectionsAtVertex(
    const std::vector<cvf::Vec3d>& polyLine,
    const std::vector<cvf::Vec3d>& polygon) {

  std::vector<cvf::Vec3d> adjustedPolygon;

  //5 times polygonScaleFactor for converting to int for clipper
  double treshold =  (1.0 / 10000.0) * 5;

  for (cvf::Vec3d polygonPoint : polygon) {

    for (size_t i = 0; i < polyLine.size() - 1; i++) {

      cvf::Vec3d linePoint1(polyLine[i].x(), polyLine[i].y(), 0.0);
      cvf::Vec3d linePoint2(polyLine[i + 1].x(), polyLine[i + 1].y(), 0.0);

      double pointDistanceFromLine =
          cvf::GeometryTools::linePointSquareDist(linePoint1,
                                                  linePoint2,
                                                  polygonPoint);

      if (pointDistanceFromLine < treshold) {

        // Calculate new polygonPoint
        cvf::Vec3d directionOfLineSegment = linePoint2 - linePoint1;

        // Finding normal to the direction of the line segment
        // in the XY plane (z=0)
        cvf::Vec3d normalToLine(-directionOfLineSegment.y(),
                                directionOfLineSegment.x(),
                                0.0);

        normalToLine.normalize();
        polygonPoint = polygonPoint + normalToLine * 0.005;
      }
    }
    adjustedPolygon.push_back(polygonPoint);
  }

  return adjustedPolygon;
}






