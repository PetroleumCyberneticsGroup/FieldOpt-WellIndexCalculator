/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2017     Statoil ASA
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

#pragma once

#include <qglobal.h>
#include <qnamespace.h>

#include <string>

class QDateTime;
class QDate;
class QTime;

//==================================================================================================
// 
//==================================================================================================
class RiaQDateTimeTools
{
public:
    static Qt::TimeSpec currentTimeSpec();

    static QDateTime fromString(const QString& dateString, const QString& format);
    static QDateTime fromYears(double years);
    
    static QDateTime addMSecs(const QDateTime& dt, double msecs);
    static QDateTime addDays(const QDateTime& dt, double days);
    static QDateTime addYears(const QDateTime& dt, double years);

    static QDateTime epoch();

    static QDateTime createUtcDateTime();
    static QDateTime createUtcDateTime(const QDate& date);
    static QDateTime createUtcDateTime(const QDate& date, const QTime& time);

private:
    static quint64  secondsInDay();
    static quint64  secondsInYear();
};
