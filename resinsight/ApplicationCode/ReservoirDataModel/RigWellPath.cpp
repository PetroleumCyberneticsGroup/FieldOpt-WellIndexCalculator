/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2011-2012 Statoil ASA, Ceetron AS
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

#include "RigWellPath.h"
#include "cvfGeometryTools.h"

//#include "resinsight/ApplicationCode/ReservoirDataModel/cvfGeometryTools.h"

#include "cvfGeometryTools.h"

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
RigWellPath::RigWellPath()
    : m_hasDatumElevation(false),
    m_datumElevation(0.0)
{
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RigWellPath::setDatumElevation(double value)
{
    m_hasDatumElevation = true;
    m_datumElevation = value;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
bool RigWellPath::hasDatumElevation() const
{
    return m_hasDatumElevation;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
double RigWellPath::datumElevation() const
{
    return m_datumElevation;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
cvf::Vec3d RigWellPath::interpolatedPointAlongWellPath(double measuredDepth) const
{
    cvf::Vec3d wellPathPoint = cvf::Vec3d::ZERO;

    size_t i = 0;
    while (i < m_measuredDepths.size() && m_measuredDepths.at(i) < measuredDepth )
    {
        i++;
    }

    if (m_measuredDepths.size() > i)
    {
        if (i == 0)
        {
            //For measuredDepth same or lower than first point, use this first point
            wellPathPoint = m_wellPathPoints.at(0);
        }
        else
        {
            //Do interpolation
            double stepsize = (measuredDepth - m_measuredDepths.at(i-1)) / 
                                        (m_measuredDepths.at(i) - m_measuredDepths.at(i - 1));
            wellPathPoint = m_wellPathPoints.at(i - 1) + stepsize * (m_wellPathPoints.at(i) - m_wellPathPoints.at(i-1));
        }
    }
    else
    {
        //Use endpoint if measuredDepth same or higher than last point
        wellPathPoint = m_wellPathPoints.at(i-1);
    }


    return wellPathPoint;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
double RigWellPath::wellPathAzimuthAngle(const cvf::Vec3d& position) const
{
    size_t closestIndex = cvf::UNDEFINED_SIZE_T;
    double closestDistance = cvf::UNDEFINED_DOUBLE;

    for (size_t i = 1; i < m_wellPathPoints.size(); i++)
    {
        cvf::Vec3d p1 = m_wellPathPoints[i - 1];
        cvf::Vec3d p2 = m_wellPathPoints[i - 0];

        double candidateDistance = cvf::GeometryTools::linePointSquareDist(p1, p2, position);
        if (candidateDistance < closestDistance)
        {
            closestDistance = candidateDistance;
            closestIndex = i;
        }
    }

    //For vertical well (x-component of direction = 0) returned angle will be 90. 
    double azimuthAngleDegrees = 90.0;

    if (closestIndex != cvf::UNDEFINED_DOUBLE)
    {
        cvf::Vec3d p1;
        cvf::Vec3d p2;

        if (closestIndex > 0)
        {
            p1 = m_wellPathPoints[closestIndex - 1];
            p2 = m_wellPathPoints[closestIndex - 0];
        }
        else
        {
            p1 = m_wellPathPoints[closestIndex + 1];
            p2 = m_wellPathPoints[closestIndex + 0];
        }

        cvf::Vec3d direction = p2 - p1;

        if (fabs(direction.y()) > 1e-5)
        {
            double atanValue = direction.x() / direction.y();
            double azimuthRadians = atan(atanValue);
            azimuthAngleDegrees = cvf::Math::toDegrees(azimuthRadians);
        }
    }

    return azimuthAngleDegrees;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
void RigWellPath::twoClosestPoints(const cvf::Vec3d& position, cvf::Vec3d* p1, cvf::Vec3d* p2) const
{
    CVF_ASSERT(p1 && p2);

    size_t closestIndex = cvf::UNDEFINED_SIZE_T;
    double closestDistance = cvf::UNDEFINED_DOUBLE;

    for (size_t i = 1; i < m_wellPathPoints.size(); i++)
    {
        cvf::Vec3d p1 = m_wellPathPoints[i - 1];
        cvf::Vec3d p2 = m_wellPathPoints[i - 0];

        double candidateDistance = cvf::GeometryTools::linePointSquareDist(p1, p2, position);
        if (candidateDistance < closestDistance)
        {
            closestDistance = candidateDistance;
            closestIndex = i;
        }
    }

    if (closestIndex != cvf::UNDEFINED_DOUBLE)
    {
        if (closestIndex > 0)
        {
            *p1 = m_wellPathPoints[closestIndex - 1];
            *p2 = m_wellPathPoints[closestIndex - 0];
        }
        else
        {
            *p1 = m_wellPathPoints[closestIndex + 1];
            *p2 = m_wellPathPoints[closestIndex + 0];
        }
    }
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
std::pair<std::vector<cvf::Vec3d>, std::vector<double> > RigWellPath::clippedPointSubset(double startMD, double endMD) const
{
    std::pair<std::vector<cvf::Vec3d>, std::vector<double> >  pointsAndMDs;
    if (m_measuredDepths.empty()) return pointsAndMDs;
    if (startMD > endMD) return pointsAndMDs;

    pointsAndMDs.first.push_back(interpolatedPointAlongWellPath(startMD));
    pointsAndMDs.second.push_back(startMD);

    for (size_t i = 0; i < m_measuredDepths.size(); ++i)
    {
        double measuredDepth = m_measuredDepths[i];
        if (measuredDepth > startMD && measuredDepth < endMD)
        {
            pointsAndMDs.first.push_back(m_wellPathPoints[i]);
            pointsAndMDs.second.push_back(measuredDepth);
        }
    }
    pointsAndMDs.first.push_back(interpolatedPointAlongWellPath(endMD));
    pointsAndMDs.second.push_back(endMD);


    return pointsAndMDs;
}

//--------------------------------------------------------------------------------------------------
/// 
//--------------------------------------------------------------------------------------------------
std::vector<cvf::Vec3d> RigWellPath::wellPathPointsIncludingInterpolatedIntersectionPoint(double intersectionMeasuredDepth) const
{
    std::vector<cvf::Vec3d> points;
    if (m_measuredDepths.empty()) return points;

    cvf::Vec3d interpolatedWellPathPoint = interpolatedPointAlongWellPath(intersectionMeasuredDepth);

    for (size_t i = 0; i < m_measuredDepths.size() - 1; i++)
    {
        if (m_measuredDepths[i] == intersectionMeasuredDepth)
        {
            points.push_back(m_wellPathPoints[i]);
        }
        else if (m_measuredDepths[i] < intersectionMeasuredDepth)
        {
            points.push_back(m_wellPathPoints[i]);
            if (m_measuredDepths[i + 1] > intersectionMeasuredDepth)
            {
                points.push_back(interpolatedWellPathPoint);
            }
        }
        else if (m_measuredDepths[i] > intersectionMeasuredDepth)
        {
            if (i == 0)
            {
                points.push_back(interpolatedWellPathPoint);
            }
            else
            {
                points.push_back(m_wellPathPoints[i]);
            }
        }
    }
    points.push_back(m_wellPathPoints.back());

    return points;
}

