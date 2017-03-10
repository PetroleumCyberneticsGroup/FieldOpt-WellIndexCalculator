/******************************************************************************
   Copyright (C) 2016-2017 Mathias C. Bellout <mathias.bellout@ntnu.no>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef FIELDOPT_TEST_RESOURCE_WIC_WIDATA_H
#define FIELDOPT_TEST_RESOURCE_WIC_WIDATA_H

#include <QStringList>
#include <QString>
#include <QDir>
#include <QDirIterator>
#include <QtCore/QString>
#include <QList>
#include <Eigen/Dense>

#include <QFile>
#include <QTextStream>
#include <fstream>
#include <QProcess>

using namespace Eigen;

namespace TestResources {
namespace WIBenchmark {

class WIData {
 public:
  // Constructor
  WIData(){};

  // Functions
  void ReadCOMPDAT(QString file_name, QString dir_list);
  void ReadXYZ(QString file_name);
  void PrintIJKData(QString file_name);
  void PrintWCFData(QString file_name);
  void CalculateWCF(QString file_root);
  void PrintCOMPDATPlot(QString file_root);

  // Variables:
  Matrix<int, Dynamic, 4> IJK;
  Matrix<double,Dynamic,1> WCF;

  Matrix<int, Dynamic, 4> IJKN;
  Matrix<double,Dynamic,1> WCFN;

  Matrix<double,1,6> XYZd;
  QStringList XYZc;
  QString XYZh, XYZt;

  QStringList name;
  QString dir_name;
  QString data_tag;

  QString grid_file;
  QString tex_file;
  QString well_name = "TW01";
  QString radius = QString::number(0.1905/2);

  bool debug_ = false;

 private:
};

void WIData::PrintCOMPDATPlot(QString file_root){

    QString csv_file = file_root + ".csv";
    QString pdf_file = file_root + ".pdf";
    // CHANGE AFTER COPYING TOOLS DIR TO BUILD FOLDER AT COMPILE TIME
    QString command = "/home/bellout/git/PCG/FieldOpt/tools"
        "/python_scripts/compdat_plot/create_compdat_plot.py "
        + csv_file + " " + pdf_file + " 60 60";

    // LAUNCH PLOT MAKER
    QProcess plot_process;
    plot_process.start(command);
    plot_process.waitForFinished();
}

void WIData::CalculateWCF(QString file_root){

    bool debug_ = true;

    // CSV FORMAT
    XYZh = XYZc[0] + " " + XYZc[1] + " " + XYZc[2];
    XYZt = XYZc[3] + " " + XYZc[4] + " " + XYZc[5];
    QString command_csv = "./wicalc --grid "
        + grid_file
        + " --heel " + XYZh
        + " --toe " + XYZt
        + " --radius " + radius;

    // LAUNCH WELL INDEX CALCULATOR (CSV FORMAT)
    QProcess wic_process_csv;
    wic_process_csv.start(command_csv);
    wic_process_csv.waitForFinished();

    // READ OUTPUT FROM QProcess COMMAND + CLOSE PROCESSES
    QString wic_process_csv_all_output = QString::fromLatin1(
        wic_process_csv.readAll().data());
    wic_process_csv.close();
    Utilities::FileHandling::WriteStringToFile(
        wic_process_csv_all_output,
        file_root + ".csv");

    // COMPDAT FORMAT
    QString command = "./wicalc --grid "
        + grid_file
        + " --heel " + XYZc[0] + " " + XYZc[1] + " " + XYZc[2]
        + " --toe "  + XYZc[3] + " " + XYZc[4] + " " + XYZc[5]
        + " --radius " + radius
        + " --compdat "
        + " --well " + well_name;

    if (debug_){
        std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl;
        std::cout << "command:" << command.toStdString() << std::endl;
    }

    // LAUNCH WELL INDEX CALCULATOR
    QProcess wic_process;
    wic_process.start(command);
    wic_process.waitForFinished();

    // READ OUTPUT FROM QProcess COMMAND + CLOSE PROCESSES
    QByteArray wic_all_output  = wic_process.readAll();
    QByteArray wic_error_output = wic_process.readAllStandardError();
    QByteArray wic_standard_output = wic_process.readAllStandardOutput();
    wic_process.close();

    // CONVERT OUTPUT TO QSTRING
    QString all_output = QString::fromLatin1(wic_all_output.data());
    QString standard_output = QString::fromLatin1(wic_standard_output.data());
    QString error_output = QString::fromLatin1(wic_error_output.data());
    Utilities::FileHandling::WriteStringToFile(
        all_output, file_root + ".compdat");

    QStringList lines = all_output.split(QRegExp("[\r\n]"));
    QStringList fields;
    Matrix<int, 1, 4> temp_IJK;
    Matrix<int, Dynamic, 4> IJK_stor;
    std::vector<double> wcf;

    foreach(QString line, lines){

        if (line.contains("OPEN")) {

            // Read IJK values from current line
            fields = line.split(QRegExp("\\s+"));
            temp_IJK << fields[2].toInt(), fields[3].toInt(),
                fields[4].toInt(), fields[5].toInt();

            // Store IJK values
            Matrix<int, Dynamic, 4>
                IJK_curr(IJK_stor.rows() + temp_IJK.rows(), 4);
            IJK_curr << IJK_stor, temp_IJK;
            IJK_stor = IJK_curr;

            // Store well connection factor values
            wcf.push_back(fields[8].toDouble());
        }
    }

    IJKN = IJK_stor;
    WCFN = Map<Matrix<double, Dynamic, 1>>(wcf.data(), wcf.size());

    if (debug_){
        std::cout << "all output:" << all_output.toStdString() << std::endl;
        std::cout << "standard output:" << standard_output.toStdString() << std::endl;
        std::cout << "error output:" << error_output.toStdString() << std::endl;

        std::cout << "IJKN:" << IJKN << std::endl;
        std::cout << "WCFN:" << WCFN << std::endl;
        std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
    }
}

void WIData::PrintIJKData(QString file_name) {

    std::ofstream file;
    file.open(file_name.toStdString());

    if (!file.is_open()){
        std::cerr << "Cannot open file '"
                  << file_name.toStdString()
                  << "' for writing." << std::endl;
    }

    file << std::fixed;
    file << IJK;
    file.close();
}

void WIData::PrintWCFData(QString file_name) {

    std::ofstream file;
    file.open(file_name.toStdString());

    if (!file.is_open()){
        std::cerr << "Cannot open file '"
                  << file_name.toStdString()
                  << "' for writing." << std::endl;
    }

    file << std::fixed;
    file << WCF;
    file.close();
}

void WIData::ReadCOMPDAT(QString file_name,
                         QString dir_list){

    // READ IJK AND WCF DATA
    QFile file(file_name);
    file.open(QIODevice::ReadOnly|QIODevice::Text);

    QFileInfo file_info(file_name);
    dir_name = file_info.dir().dirName();
    // std::cout << "dir_name:" << dir_name.toStdString();

    QTextStream in(&file);
    QStringList in_fields;

    Matrix<int, 1, 4> temp_IJK;
    Matrix<int, Dynamic, 4> IJK_stor;
    std::vector<double> wcf;
    QString well_name;

    while(!in.atEnd()) {

        QString line = in.readLine();

        if (line.contains("OPEN")) {

            // TODO: is there a more robust way to read line such that
            // indices are not wrong if a change in spacing, for example?
            // how to remove the space that the line begins with?
            in_fields = line.split(QRegExp("\\s+"));

            // Read & store well name from current line
            well_name.append(in_fields[1]);

            // Read IJK values from current line
            temp_IJK << in_fields[2].toInt(), in_fields[3].toInt(),
                in_fields[4].toInt(), in_fields[5].toInt();

            // Store IJK values
            Matrix<int, Dynamic, 4>
                IJK_curr(IJK_stor.rows() + temp_IJK.rows(), 4);
            IJK_curr << IJK_stor, temp_IJK;
            IJK_stor = IJK_curr;

            // Store well connection factor values
            wcf.push_back(in_fields[8].toDouble());
        };
    }

    file.close();

    IJK = IJK_stor;
    WCF = Map<Matrix<double, Dynamic, 1>>(wcf.data(), wcf.size());
}

void WIData::ReadXYZ(QString file_name){

    // READ XYZ DATA
    QFile xyz_file(file_name + ".xyz");
    xyz_file.open(QIODevice::ReadOnly|QIODevice::Text);

    QTextStream xyz_in(&xyz_file);
    QStringList xyz_in_fields;
    std::vector<double> xyz_d;
    XYZc.clear();

    while(!xyz_in.atEnd()) {

        QString line = xyz_in.readLine();
        xyz_in_fields = line.split(QRegExp("[\r\n\t]"));

        if (!line.contains("TW01")) {
            // Store xyz values
            for(int ii = 0; ii < xyz_in_fields.size(); ++ii){
                xyz_d.push_back(xyz_in_fields[ii].toDouble());
                XYZc << xyz_in_fields[ii];
            }
        }
    }

    xyz_file.close();

    XYZd = Map<Matrix<double,1,6>>(xyz_d.data(), xyz_d.size());

    // Debug: check read process for XYZ data is OK
    if (debug_){
        std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl;
        std::cout << "--" << data_tag.toStdString() << " data--" << std::endl;
        std::cout << std::setfill(' ');
        // DOUBLE DATA
        IOFormat CleanFmt(FullPrecision, 0, " ", "\n", "[", "]");
        std::cout << "XYZd: " << XYZd.format(CleanFmt) << std::endl;
        // CHAR DATA
        std::cout << "XYZc: ";
        for(int ii = 0; ii < XYZc.size(); ++ii){
            std::cout << XYZc[ii].toStdString() << " ";
        }
        std::cout << std::endl;
        std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
    }
}
}
}

#endif //FIELDOPT_TEST_RESOURCE_WIC_WIDATA_H