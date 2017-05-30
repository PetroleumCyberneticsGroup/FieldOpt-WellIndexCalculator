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



#ifndef FIELDOPT_TEST_RESOURCE_DIFF_FUNCTIONS_H
#define FIELDOPT_TEST_RESOURCE_DIFF_FUNCTIONS_H

#include "test_resource_wic_widata.h"
#include <Eigen/Dense>
#include <typeinfo>

#include <iomanip>
#include <QProcess>
#include <QVector>

using namespace Eigen;

namespace TestResources {
namespace WIBenchmark {


/*!
 * \brief
 * \param
 * \return
 */
std::vector<double> ConvertEigenToStd(Matrix<double, Dynamic, 1> eigen_vector){

    // feed eigen vector values to std:vector
    std::vector<double> std_vector;
    std_vector.resize(eigen_vector.size());
    Matrix<double, Dynamic, 1>::Map(&std_vector[0], eigen_vector.size()) = eigen_vector;

    return std_vector;
};

/*!
 * \brief
 * \param
 * \return
 */
Matrix<double, Dynamic, 1> ConvertStdToEigen(std::vector<double> std_vector){

    // feed std:vector values to eigen vector
    Matrix<double, Dynamic, 1> eigen_vector =
        Matrix<double, Dynamic,1>::Map(std_vector.data(), std_vector.size());

    return eigen_vector;
};

/*!
 * \brief
 * \param
 * \return
 */
double GetEpsIJK(){
    double eps = 1e-12; // tolerance for comparison, NOT for removing rows
    return eps;
}

/*!
 * \brief
 * \param
 * \return
 */
double GetEpsWCF(){
    double eps = 0.01; // tolerance for comparison, NOT for removing rows
    return eps;
}

/*!
 * \brief
 * \param
 * \return
 */
bool DiffVectorLength(WIData va, WIData vb){

    bool vector_diff = true;
    if (va.IJK.rows() != vb.IJK.rows()){
        vector_diff = false;
    }
    return vector_diff;
}

/*!
 * \brief
 * \param
 * \return
 */
WIData GetLongestVector(WIData va, WIData vb){

    if (va.IJK.rows() > vb.IJK.rows()){
        return va;
    }else{
        return vb;
    }
}

/*!
 * \brief
 * \param
 * \return
 */
WIData GetShortestVector(WIData va, WIData vb){

    if (va.IJK.rows() < vb.IJK.rows()){
        return va;
    }else{
        return vb;
    }
}

/*!
 * \brief
 * \param
 * \return
 */
double GetColumnOffset(Matrix<double,Dynamic,1> va,
                       Matrix<double,Dynamic,1> vb,
                       Matrix<double,Dynamic,1> vdiff){

    // accuracy_magnitude: norm of difference vector (column offset)
    return vdiff.norm();
}

/*!
 * \brief
 * \param
 * \return
 */
double GetColumnCosine(Matrix<double,Dynamic,1> va,
                       Matrix<double,Dynamic,1> vb,
                       Matrix<double,Dynamic,1> vdiff){

    // cosine similarity (cosine measure)
    return va.dot(vb) / (va.norm() * vb.norm());
}

/*!
 * \brief Remove values from vdiff that are above threshold value
 * \param
 * \return
 */
Matrix<double,Dynamic,1> ApplyThresh(Matrix<double, Dynamic, 1> vvector,
                                     double threshold){

    // tranfer only values below threshold
    std::vector<double> vout_std;
    for (int ii=0; ii < vvector.rows(); ++ii){
        double vtemp = vvector(ii);
        if (vtemp < threshold && vtemp > 1/threshold){
            vout_std.push_back(vtemp);
        }
    }

    return ConvertStdToEigen(vout_std);
}

/*!
 * \brief
 * \param
 * \return
 */
double GetColumnMedian(Matrix<double, Dynamic, 1> va,
                       Matrix<double, Dynamic, 1> vb,
                       Matrix<double, Dynamic, 1> vdiff,
                       int apply_th,
                       double threshold) {

    Matrix<double, Dynamic, 1> vratio = va.cwiseQuotient(vb);
    if (apply_th > 0){
        vratio = ApplyThresh(vratio, threshold);
    }

    auto std_vratio = ConvertEigenToStd(vratio);

    // median of RMS/PCG ratio
    int vsz = (int)std_vratio.size();
    std::sort(std_vratio.data(), std_vratio.data()+vsz);

    // indexing ok since index numbering is 1 less than vector length
    auto median = (vsz % 2 == 1) ?
                  // odd; sz=7->vz/2=3 in [0 1 2 (3) 4 5 6]->i=3
                  std_vratio[vsz / 2] :
                  // even; sz=6->vz/2=3 in [0 1 (2 3) 4 5]->i1=2,i2=3
                  (std_vratio[vsz / 2 - 1] + std_vratio[vsz / 2]) / 2;


    if (apply_th < 2) {
        return median;
    }
    else {
        return (double)va.rows() - (double)vratio.rows();
    }
}

/*!
 * \brief
 * \param
 * \return
 */
double GetColumnMean(Matrix<double, Dynamic, 1> va,
                     Matrix<double, Dynamic, 1> vb,
                     Matrix<double, Dynamic, 1> vdiff,
                     int apply_th,
                     double threshold) {

    Matrix<double,Dynamic,1> vratio = va.cwiseQuotient(vb);
    if (apply_th > 0){
        vratio = ApplyThresh(vratio, threshold);
    }

    // mean
    if (apply_th < 2) {
        return vratio.mean();
    }
    else {
        return (double)va.rows() - (double)vratio.rows();
    }
}

/*!
 * \brief
 * \param
 * \return
 */
void RemoveRowsLowWCF(WIData &data) {

    bool debug = false;
    int jj = 0;
    int rem_count = 0;
    int neg_cnt = 0;
    int tol_cnt = 0;

    // MSG OUT
    QString str_out, wcf_str, temp_str, msg, nl, lstr_out;
    lstr_out = "\n-----------------------------------"
        "---------------------------------------------";
    str_out = lstr_out + "\n>>> If any, start removing "
        "rows with low WCF for " + data.data_tag + " data.";

    // DEBUG
    if(debug){
        std::cout << std::endl << "size of " << data.data_tag.toStdString()
                  << " data entering: " << "nRowsWCF=" << data.WCF.rows()
                  << ", nRowsIJK=" << data.IJK.rows();
    }

    // FIND ROW INDICES OF ORIG MATRIX THAT ARE ABOVE TOL, OUTPUT MSG IF NOT
    std::vector<int> tol_ind;
    double tol = .01; // tolerance for removing row

    for (int ii = 0; ii < data.WCF.rows(); ++ii) {

        auto wcf_num = data.WCF.row(ii).value();
        wcf_str.sprintf("%5.4f", wcf_num);

        if (wcf_num > tol) {
            tol_ind.push_back(ii);
        }else{
            // MSG OUT
            msg = (wcf_num >= 0) ? "is lower than set tolerance" : "is negative!";
            if (wcf_num >= 0) { tol_cnt += 1; } else { neg_cnt += 1; }
            if (rem_count <= 10) {
                temp_str = "Removing row " + QString::number(ii) + " from "
                    + data.data_tag + " data b/c WCF " + msg + " ("
                    + wcf_str + " < " + QString::number(tol) + ")";
                str_out.append("\n" + temp_str);
            }
            rem_count += 1;
        }
    }

    // MAKE TEMP MATRICES TO STORE ORIG ROW ABOVE TOL
    WIData data_temp;
    data_temp.IJK = Matrix<int, Dynamic, 4>::Zero(tol_ind.size(),4);
    data_temp.WCF = Matrix<double,Dynamic,1>::Zero(tol_ind.size());

    // COPY SELECTED ROWS
    for (int ii = 0; ii < tol_ind.size(); ++ii) {
        // std::cout << std::endl << "row:" << ii << ", tol_ind[ii] " << tol_ind[ii];
        data_temp.IJK.row(ii) = data.IJK.row(tol_ind[ii]);
        data_temp.WCF.row(ii) = data.WCF.row(tol_ind[ii]);
    }

    // TRANSFER VALUES
    data.IJK = data_temp.IJK;
    data.WCF = data_temp.WCF;

    // DEBUG
    if(debug){
        std::cout << std::endl << "size of " << data.data_tag.toStdString()
                  << " data leaving: " << "nRowsWCF=" << data.WCF.rows()
                  << ", nRowsIJK=" << data.IJK.rows();
    }

    // MSG OUT
    if (tol_cnt + neg_cnt > 10) {
        str_out = str_out + "\n+" + QString::number(tol_cnt + neg_cnt - 10)
            + " other rows removed b/c WCF " + msg + "\n";
    }
    if (rem_count < 1) { str_out = str_out + " [None removed.]"; }

    std::cout << std::endl << "\033[1;33m" << str_out.toStdString() << "\033[0m";
    Utilities::FileHandling::WriteLineToFile(str_out, data.tex_file);

}

/*!
 * \brief
 * \param
 * \return
 */
double GetColumnAccuracyElements(Matrix<double,Dynamic,1> col_vector){

    // Accuracy_elements: fraction of elements in column which
    // are zero up to given tolerance
    double nrows = col_vector.rows();
    double nrows_nz = 0;

    for( int ii=0; ii < nrows; ++ii ) {
        Matrix<double,1,1> row_element;
        row_element << col_vector[ii];
        if(row_element.isZero(GetEpsWCF())){ nrows_nz += 1; }
    }

    return nrows_nz / nrows;
}

/*!
 * \brief
 * \param
 * \return
 */
QList<Matrix<double,1,4>> CheckColumnwiseDiff(WIData va, WIData vb, WIData vdiff){

    auto va_ = va.IJK.cast<double>();
    auto vb_ = vb.IJK.cast<double>();
    auto vdiff_ = vdiff.IJK.cast<double>();

    Matrix<double,1,4> IJK_accuracy_elements;
    Matrix<double,1,4> IJK_column_offset;
    Matrix<double,1,4> IJK_column_cosine;

    for( int ii=0; ii < vdiff_.cols(); ++ii ) {
        IJK_accuracy_elements(ii) = GetColumnAccuracyElements(vdiff_.col(ii));
        IJK_column_offset(ii) = GetColumnOffset(va_.col(ii), vb_.col(ii), vdiff_.col(ii));
        IJK_column_cosine(ii) = GetColumnCosine(va_.col(ii), vb_.col(ii), vdiff_.col(ii));
    }

    QList<Matrix<double,1,4>> IJK_accuracy_list;
    IJK_accuracy_list.append(IJK_accuracy_elements);
    IJK_accuracy_list.append(IJK_column_offset);
    IJK_accuracy_list.append(IJK_column_cosine);

    QString temp_str, str_out, nl;
    QString lstr_out = "\n--------------------------------------------------------------------------------";

    // Zero element fraction
    str_out = "\nElement accuracy: fraction of zero (<tol) elements in diff. column (1=best)";
    std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    str_out  = "Element accuracy I column: " + QString::number(IJK_accuracy_list[0][0]);
    temp_str = "Element accuracy J column: " + QString::number(IJK_accuracy_list[0][1]);
    str_out.append("\n" + temp_str);
    temp_str = "Element accuracy K column: " + QString::number(IJK_accuracy_list[0][2]);
    str_out.append("\n" + temp_str);
    // temp_str = "Element accuracy K2 column: " + QString::number(IJK_accuracy_list[0][3]);
    // str_out.append("\n" + temp_str);
    std::cout << str_out.toStdString() << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    // Column offset
    str_out = "Column IJK offset: norm of diff. vector [Euclidean dist.]  (0=best)";
    std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    str_out  = "Column offset I column:  " + QString::number(IJK_accuracy_list[1][0]);
    temp_str = "Column offset J column:  " + QString::number(IJK_accuracy_list[1][1]);
    str_out.append("\n" + temp_str);
    temp_str = "Column offset K column:  " + QString::number(IJK_accuracy_list[1][2]);
    str_out.append("\n" + temp_str);
    // temp_str = "Column offset K2 column:  " + QString::number(IJK_accuracy_list[1][3]);
    // str_out.append("\n" + temp_str);
    std::cout << str_out.toStdString() << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    // Cosine measure
    str_out = "Column IJK cosine similarity: angle b/e vectors (1=parallel=best)";
    std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    str_out  = "Column cosine I column: " + QString::number(IJK_accuracy_list[2][0]);
    temp_str = "Column cosine J column: " + QString::number(IJK_accuracy_list[2][1]);
    str_out.append("\n" + temp_str);
    temp_str = "Column cosine K column: " + QString::number(IJK_accuracy_list[2][2]);
    str_out.append("\n" + temp_str);
    // temp_str = "Column cosine K2 column: " + QString::number(IJK_accuracy_list[2][3]);
    // str_out.append("\n" + temp_str);
    std::cout << str_out.toStdString() << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    return IJK_accuracy_list;
}

/*!
 * \brief
 * \param
 * \return
 */
template<typename T, typename V> void CheckRowwiseDiff(
    T& va_, T& vb_, V& vdiff_, QString tag, double tol, WIData data){

    // Check vector has length > 0
    if (!vdiff_.rows()>0)
        throw std::runtime_error("Difference vector (vdiff_) has no elements");

    auto vrel_ = va_.cwiseQuotient(vb_);

    // Output msg
    QString str_out;
    str_out = "Testing: " + tag + " (rowwise) >> values differ at the following rows:";
    std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, data.tex_file);

    // Loop over each row
    for( int ii=0; ii < vdiff_.rows(); ++ii ) {

        // Setup row
        auto vdiff_row = vdiff_.row(ii);

        // print out to double vectors for QString treatment later on
        std::vector<double> va_d, vb_d, vdiff_d, vrel_d;
        QStringList labels;
        labels << "RMS: " << "PCG: " << "DFF: " << "RMS/PCG: ";
        string frmt;
        int ncols;

        for( int jj = 0; jj < vdiff_row.size(); ++jj ) {
            va_d.push_back(va_(ii,jj));
            vb_d.push_back(vb_(ii,jj));
            vdiff_d.push_back(vdiff_(ii,jj));
            vrel_d.push_back(vrel_(ii,jj));
        }

        if (tag.compare("IJK")==0){
            frmt = "%2.0f ";
            ncols = va_d.size() - 1; // SKIP LAST COLUMN
        }else if(tag.compare("WCF")==0){
            frmt = "%7.3f    ";
            ncols = va_d.size();
        }

        // If difference is larger than zero by a given tolerance
        if (!vdiff_row.isZero(tol)){

            QString num_str, txt_str;
            num_str.sprintf("row %3.0f:  ", (double) ii);
            txt_str.append(num_str);

            txt_str.append(labels[0]); // RMS DATA
            for( int jj = 0; jj < ncols; ++jj ) {
                txt_str.append(num_str.sprintf(frmt.c_str(), va_d[jj]));
            }

            txt_str.append(labels[1]); // PCG DATA
            for( int jj = 0; jj < ncols; ++jj ) {
                txt_str.append(num_str.sprintf(frmt.c_str(), vb_d[jj]));
            }

            txt_str.append(labels[2]); // DFF DATA
            for( int jj = 0; jj < ncols; ++jj ) {
                txt_str.append(num_str.sprintf(frmt.c_str(), vdiff_d[jj]));
            }

            txt_str.append(labels[3]); // DFF DATA
            for( int jj = 0; jj < ncols; ++jj ) {
                txt_str.append(num_str.sprintf(frmt.c_str(), vrel_d[jj]));
            }

            std::cout << txt_str.toStdString() << std::endl;
            Utilities::FileHandling::WriteLineToFile(txt_str, data.tex_file);
        }
    }
}

/*!
 * \brief
 * \param
 * \return
 */
void CheckRowwiseDiffIJK(WIData va, WIData vb, WIData vdiff){
    auto vdiff_ = vdiff.IJK;
    auto va_ = va.IJK;
    auto vb_ = vb.IJK;
    QString tag = "IJK";
    CheckRowwiseDiff(va_, vb_, vdiff_, tag, GetEpsIJK(), va);
}

/*!
 * \brief
 * \param
 * \return
 */
void CheckRowwiseDiffWCF(WIData va, WIData vb, WIData vdiff){
    auto vdiff_ = vdiff.WCF;
    auto va_ = va.WCF;
    auto vb_ = vb.WCF;
    QString tag = "WCF";
    CheckRowwiseDiff(va_, vb_, vdiff_, tag, GetEpsWCF(), va);
}

/*!
 * \brief
 * \param
 * \return
 */
WIData CompareIJK(WIData va, WIData vb){

    bool debug_ = false;

    WIData vdiff;
    vdiff.IJK = va.IJK - vb.IJK;
    QList<Matrix<double,1,4>> IJK_accuracy_list;

    // Check if IJK values computed by RMS and PCG WIC are equal
    // If not, output differing rows
    QString str_out;
    QString lstr_out = "\n--------------------------------------------------------------------------------";
    QString tol;
    tol.sprintf("%5.3e",GetEpsIJK());

    if (vdiff.IJK.isZero(GetEpsIJK())){
        str_out = lstr_out + "\nIJK values match exactly for this well.";
        std::cout << "\033[1;32m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

    }else{
        str_out = lstr_out + "\nIJK values are NOT the same for this well.";
        std::cout << "\033[1;35m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Output row differences (i.e., I, J and K columns combined)
        CheckRowwiseDiffIJK(va,vb,vdiff);
        // Output difference for individual columns
        IJK_accuracy_list = CheckColumnwiseDiff(va,vb,vdiff);
    }

    if (debug_){
        int nrows = (vdiff.IJK.rows() > 10) ? 10 : vdiff.IJK.rows();
        std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl
                  << "(WIDataPCG.IJK - WIDataPCG.IJK)= " << std::endl
                  << vdiff.IJK.block(0,0,nrows,4)
                  << std::endl << "..." << std::endl;
        std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
    }

    return vdiff;
}

/*!
 * \brief
 * \param
 * \return
 */
WIData CompareWCF(WIData &va, WIData &vb, int ii) {

    double threshold = 2.0;
    WIData vdiff;
    vdiff.WCF = va.WCF - vb.WCF;

    QString str_out, str_smry;
    QString lstr_out = "\n--------------------------------------------------------------------------------";
    QString tol;
    tol.sprintf("%5.3f", GetEpsWCF());

    if(va.WCF.isApprox(vb.WCF, GetEpsWCF())){
        str_out = lstr_out + "\nWCF values match exactly for this well (WCF tol = " + tol + ").";
        std::cout << "\033[1;32m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        vdiff.WCF_accuracy_list << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0;
//        std::cout << "ii:" << ii << "-> va.test_IJK_removed.size():" << va.test_IJK_removed.size() << std::endl;
//        std::cout << "ii:" << ii << "-> va.test_IJK_removed[ii].size():" << va.test_IJK_removed[ii].size() << std::endl;
//        std::cout << "ii:" << ii << "-> va.test_IJK_removed[ii][1].size():" << va.test_IJK_removed[ii][1].size() << std::endl;
//        std::cout << "ii:" << ii << "-> va.test_IJK_removed[ii][3].size():" << va.test_IJK_removed[ii][3].size() << std::endl;

        if (va.test_IJK_removed[ii][1].size() + va.test_IJK_removed[ii][3].size() > va.max_sup) {
            str_smry = va.dir_name.replace("_","\\_") + " & " + "0.000 & 0.000 \\\\";
        }
        else {
            str_smry = va.dir_name.replace("_","\\_") + " & " + "1.000 & 1.000 \\\\";
        };

    }else{
        str_out = lstr_out + "\nWCF values are NOT the same for this well (WCF tol = " + tol + ").";
        std::cout << "\033[1;35m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Output general difference (i.e., for I, J and K columns)
        CheckRowwiseDiffWCF(va,vb,vdiff);

        // Output difference for column
        vdiff.WCF_accuracy_list.append(GetColumnAccuracyElements(vdiff.WCF)); // 0
        vdiff.WCF_accuracy_list.append(GetColumnOffset(va.WCF, vb.WCF, vdiff.WCF)); // 1
        vdiff.WCF_accuracy_list.append(GetColumnCosine(va.WCF, vb.WCF, vdiff.WCF)); // 2

        // Mean (w/o and w/ threshold, and and removed values)
        vdiff.WCF_accuracy_list.append(GetColumnMean(va.WCF, vb.WCF, vdiff.WCF, 0, threshold)); // 3
        vdiff.WCF_accuracy_list.append(GetColumnMean(va.WCF, vb.WCF, vdiff.WCF, 1, threshold)); // 4
        vdiff.WCF_accuracy_list.append(GetColumnMean(va.WCF, vb.WCF, vdiff.WCF, 2, threshold)); // 5

        // Median (w/o and w/ threshold, and and removed values)
        vdiff.WCF_accuracy_list.append(GetColumnMedian(va.WCF, vb.WCF, vdiff.WCF, 0, threshold)); // 6
        vdiff.WCF_accuracy_list.append(GetColumnMedian(va.WCF, vb.WCF, vdiff.WCF, 1, threshold)); // 7
        vdiff.WCF_accuracy_list.append(GetColumnMedian(va.WCF, vb.WCF, vdiff.WCF, 2, threshold)); // 8

        // Zero element fraction
        str_out = "\nElement accuracy: fraction of zero (<tol) elements in diff. column (1=best)";
        std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        auto str_val0 = QString::number(vdiff.WCF_accuracy_list[0]);
        str_out  = "Element accuracy:  " + str_val0;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Column offset
        str_out = "Column WCF offset: norm of difference vector (0=best)";
        std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        auto str_val1 = QString::number(vdiff.WCF_accuracy_list[1]);
        str_out  = "Column offset:  " + str_val1;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Cosine measure
        str_out = "Column WCF cosine similarity: angle b/e vectors (1=parallel=best)";
        std::cout << "\033[1;33m" << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        auto str_val2 = QString::number(vdiff.WCF_accuracy_list[2]);
        str_out  = "Column cosine measure:  " + str_val2;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Mean (w/o threshold)
        auto str_val3 = QString::number(vdiff.WCF_accuracy_list[3]).leftJustified(5, '0', true);
        str_out  = "\nMean of RMS/PCG wi-ratios (w/o threshold): " + str_val3;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Mean (w/ threshold)
        auto str_val4 = QString::number(vdiff.WCF_accuracy_list[4]).leftJustified(5, '0', true);
        str_out = "Mean of RMS/PCG wi-ratios (w/ threshold):  " + str_val4;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Mean (num removed values)
        auto str_val5 = QString::number(vdiff.WCF_accuracy_list[5]);
        str_out  = "[Number of values removed by threshold (="
        + QString::number(threshold) + "): " + str_val5 + "]";
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Median (w/o threshold)
        auto str_val6 = QString::number(vdiff.WCF_accuracy_list[6]).leftJustified(5, '0', true);
        str_out  = "\nMedian of RMS/PCG wi-ratios (w/o threshold): " + str_val6;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Median (w/ threshold)
        auto str_val7 = QString::number(vdiff.WCF_accuracy_list[7]).leftJustified(5, '0', true);
        str_out  = "Median of RMS/PCG wi-ratios (w/ threshold): " + str_val7;
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        // Median (num removed values)
        auto str_val8 = QString::number(vdiff.WCF_accuracy_list[8]);
        str_out  = "[Number of values removed by threshold (="
        + QString::number(threshold) + "): " + str_val8 + "]";
        std::cout << str_out.toStdString() << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, va.tex_file);

        str_smry = va.dir_name.replace("_","\\_") + " & " + str_val4 + " & " + str_val6 + " \\\\";
    }

    std::cout << "Writing to file:" << va.tex_smry.toStdString() << std::endl << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_smry, va.tex_smry);

    bool debug_ = false;
    if (debug_){
        int nrows = (vdiff.WCF.rows() > 10) ? 10 : vdiff.WCF.rows();
        std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl
                  << "(WIDataPCG.WCF - WIDataPCG.WCF)= " << std::endl
                  << vdiff.WCF.block(0,0,nrows,1)
                  << std::endl << "..." << std::endl;
        std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
    }

    return vdiff;
}

/*!
 * \brief
 * \param
 * \return
 */
QVector<int> DiffTreatmentA(WIData &WIDataRMS,
                            WIData &WIDataPCG,
                            QStringList &diff_files){

    bool debug_ = false;

    // DIFF TREATMENT: IJK COMPARISON: FIND EXTRA ROWS USING diff COMMAND
    QProcess diff_process, grep_process;
    diff_process.setStandardOutputProcess(&grep_process);
    grep_process.setProcessChannelMode(QProcess::MergedChannels);

    if (debug_) {
        std::cout << "Original number of rows in RMS and PCG data:" << std::endl;
        std::cout << "WIDataPCG.IJK.rows:" << WIDataPCG.IJK.rows() << std::endl;
        std::cout << "WIDataRMS.IJK.rows:" << WIDataRMS.IJK.rows() << std::endl;
    }

    // CALL diff COMMAND USING QProcess; PIPE IT TO grep COMMAND
    // SWITCH COLUMNS SUCH THAT LONGEST COLUMN IS ALWAYS THE COLUMN TO THE RIGHT
    if (WIDataRMS.IJK.rows() > WIDataPCG.IJK.rows()) { // longer vector = right column
        diff_process.start("diff -y " + diff_files[1] + " " + diff_files[0]); // PCG[1] < RMS[0]
    } else {
        diff_process.start("diff -y " + diff_files[0] + " " + diff_files[1]); // RMS[0] < PCG[1]
    }
    grep_process.start("grep \">\" -n");
    diff_process.waitForFinished();
    grep_process.waitForFinished();

    // READ OUTPUT FROM QProcess COMMAND + CLOSE PROCESSES
    QByteArray grep_output = grep_process.readAllStandardOutput();
    diff_process.close();
    grep_process.close();

    // OBTAIN INDICES OF SUPERFLUOUS ROWS IN LONGEST COLUMN:
    // SPLIT QString SUCH THAT EACH LINE IS ONE ELEMENT IN A QStringList
    QString DataAsString = QString::fromLatin1(grep_output.data());
    QStringList textString = DataAsString.split(
        QRegExp("[\r\n]"), QString::SkipEmptyParts);
    QVector<int> sup_indices;

    // OBTAIN INDICES OF SUPERFLUOUS ROWS IN LONGEST COLUMN:
    // READ EACH LINE OF diff OUTPUT AND STORE INDICES -- IMPORTANT: INDICES OBTAINED
    // FROM DIFF START FROM 1; WE THEREFORE SUBSTRACT 1 FROM THESE TO MAKE THESE WORK
    // WITH 0-START INDEXING
    if (debug_) std::cout << std::endl << "superfluous indices in rightmost column:" << std::endl;
    for (int ii = 0; ii < textString.size(); ++ii) {
        sup_indices.append(textString[ii].split(":").first().toInt() - 1); // CONVERT TO 0-START INDEXING
        if (debug_) std::cout << ii << ":" << sup_indices[ii] << std::endl;
    }
    if (debug_) std::cout << "sup_indices.size:" << sup_indices.size() << std::endl;

    return sup_indices;
}

struct IJKData{

  Matrix<int, Dynamic, 4> IJK_PCG, IJK_RMS;
  Matrix<int, Dynamic, 1> IJK_PCG_IDX, IJK_RMS_IDX;

};

/*!
 * \brief
 * \param
 * \return
 */
IJKData DiffGrepSedProcess(QString diff_str,
                           QString grep_str,
                           QString sed_str) {

    bool debug_ = false;

    // CALL diff COMMAND USING QProcess; PIPE IT TO grep AND sed COMMAND
    QProcess diff_process, grep_process, sed_process;
    grep_process.setStandardOutputProcess(&sed_process);
    diff_process.setStandardOutputProcess(&grep_process);
    sed_process.setProcessChannelMode(QProcess::MergedChannels);
    grep_process.setProcessChannelMode(QProcess::MergedChannels);

    diff_process.start(diff_str);
    grep_process.start(grep_str);
    sed_process.start(sed_str);

    diff_process.waitForFinished();
    grep_process.waitForFinished();
    sed_process.waitForFinished();

    QByteArray sed_output = sed_process.readAllStandardOutput();

    // READ OUTPUT FROM QProcess COMMAND + CLOSE PROCESSES
    diff_process.close();
    grep_process.close();
    sed_process.close();

    QString DataAsString = QString::fromLatin1(sed_output.data());

    if (debug_) {
        std::cout << "DataAsString:\n"
                  << DataAsString.toStdString() << std::endl;
    }

    QStringList textString = DataAsString.split(
        QRegExp("[\r\n]"), QString::SkipEmptyParts);

    // SET UP PCG/RMS IJK MATRICES
    Matrix<int, Dynamic, 4> IJK_PCG, IJK_RMS;
    IJK_PCG.resize(textString.size(),4);
    IJK_RMS.resize(textString.size(),4);
    IJK_PCG.fill(0);
    IJK_RMS.fill(0);

    // SET UP PCG/RMS IJK INDEX VECTORS
    Matrix<int, Dynamic, 1> IJK_PCG_IDX, IJK_RMS_IDX;
    IJK_PCG_IDX.resize(textString.size(),1);
    IJK_RMS_IDX.resize(textString.size(),1);
    IJK_PCG_IDX.fill(0);
    IJK_RMS_IDX.fill(0);

    // GET IJK DATA FROM STRING
    for (int ii = 0; ii < textString.size(); ++ii) {
        QStringList IJK = textString[ii].split(
            QRegExp("[:\r\t\n ]+"), QString::SkipEmptyParts);


        // INSERT READ IDX VALUES
        IJK_PCG_IDX(ii,0) = IJK[0].toInt() - 1;
        IJK_RMS_IDX(ii,0) = IJK[0].toInt() - 1;

        if (debug_) { std::cout << "IJK.size(): " << IJK.size() << std::endl; }

        int ncols = (IJK.size()==9) ? 5 : 1;
        for (int jj = 0; jj < 4; ++jj) {
            IJK_PCG(ii,jj) = IJK[jj+1].toInt();
            IJK_RMS(ii,jj) = IJK[jj+ncols].toInt();
        }
    }

    if (debug_) {
        std::cout << "IJK_PCG:\n" << IJK_PCG << std::endl;
        std::cout << "IJK_RMS:\n" << IJK_RMS << std::endl;
        std::cout << "IJK_PCG_IDX:\n" << IJK_PCG_IDX << std::endl;
        std::cout << "IJK_RMS_IDX:\n" << IJK_RMS_IDX << std::endl;
    }

    IJKData IJK_IDX_DATA;
    IJK_IDX_DATA.IJK_PCG = IJK_PCG;
    IJK_IDX_DATA.IJK_RMS = IJK_RMS;
    IJK_IDX_DATA.IJK_PCG_IDX = IJK_PCG_IDX;
    IJK_IDX_DATA.IJK_RMS_IDX = IJK_RMS_IDX;

    return IJK_IDX_DATA;
}

/*!
 * \brief
 * \param
 * \return
 */
std::vector< std::vector<int> > DiffTreatmentB(WIData &WIDataRMS,
                                               WIData &WIDataPCG,
                                               QStringList &diff_files){

    bool debug_ = false;
    QString diff_str, grep_str, sed_str;

    // SWITCH COLUMNS SUCH THAT LONGEST COLUMN IS ALWAYS THE COLUMN TO THE RIGHT
    diff_str = "diff -y " + diff_files[1] + " " + diff_files[0]; // PCG[1] < RMS[0]
    grep_str = "egrep \"<|>\" -nv";
    sed_str = "sed \"s/|//g\"";
    // DIFF TREATMENT: IJK COMPARISON: FIND EXTRA ROWS USING diff COMMAND
    IJKData IJK_IDX_UNION = DiffGrepSedProcess(diff_str, grep_str, sed_str);

    grep_str = "egrep \">\" -n";
    sed_str = "sed \"s/>//g\"";
    IJKData IJK_IDX_COMP_PCG = DiffGrepSedProcess(diff_str, grep_str, sed_str);

    grep_str = "egrep \"<\" -n";
    sed_str = "sed \"s/<//g\"";
    IJKData IJK_IDX_COMP_RMS = DiffGrepSedProcess(diff_str, grep_str, sed_str);

    if (debug_) {
        std::cout << "___________________________________________" << std::endl;
        std::cout << "using RemoveSuperfluousRowsB to fix sizes\n" << std::endl;
    }

    auto IJK_PCG_NEW = IJK_IDX_UNION.IJK_PCG;
    auto IJK_RMS_NEW = IJK_IDX_UNION.IJK_RMS;

    // DEFINE INDICES FOR DIFF DATA
    std::vector<int> orig_rows_PCG_dff, added_rows_PCG_dff, orig_rows_RMS_dff, added_rows_RMS_dff;

    Matrix<int, Dynamic, 1> vEigen;
    vEigen = IJK_IDX_UNION.IJK_PCG_IDX;
    orig_rows_PCG_dff.resize((unsigned long)vEigen.size());
    Matrix<int, Dynamic, 1>::Map(&orig_rows_PCG_dff[0], vEigen.size()) = vEigen;

    vEigen = IJK_IDX_UNION.IJK_RMS_IDX;
    orig_rows_RMS_dff.resize((unsigned long)vEigen.size());
    Matrix<int, Dynamic, 1>::Map(&orig_rows_RMS_dff[0], vEigen.size()) = vEigen;

    vEigen = IJK_IDX_COMP_PCG.IJK_PCG_IDX;
    added_rows_PCG_dff.resize((unsigned long)vEigen.size());
    Matrix<int, Dynamic, 1>::Map(&added_rows_PCG_dff[0], vEigen.size()) = vEigen;

    vEigen = IJK_IDX_COMP_RMS.IJK_RMS_IDX;
    added_rows_RMS_dff.resize((unsigned long)vEigen.size());
    Matrix<int, Dynamic, 1>::Map(&added_rows_RMS_dff[0], vEigen.size()) = vEigen;

    // ===============================================================================
    // WORKS BUT DOING EVERTHING THROUGH DIFF COMMANDS

    // DEFINE INDICES FOR COMPUTED DATA
    std::vector<int> orig_rows_PCG, added_rows_PCG, orig_rows_RMS, added_rows_RMS;
    std::vector< std::vector<int> > sup_indices;

    // // FIND WHICH ROWS HAVE BEEN REMOVED FROM PCG DATA
    // PCG: GET INDICES OF ROWS THAT ARE IN ORIGINAL IJK DATA
    for (int jj = 0; jj < IJK_PCG_NEW.rows(); ++jj) {

        for (int ii = 0; ii < WIDataPCG.IJK.rows(); ++ii) {
            // FIND WHICH ORIG ROW FROM PCG DATA THAT ARE NOW IN NEW IJK MATRIX
            if (IJK_PCG_NEW.row(jj).cwiseEqual(WIDataPCG.IJK.row(ii)).count() == 4) {

                if (debug_) { std::cout << "new PCG IJK row (jj=" << jj
                                        << ") == old PCG IJK row (ii=" << ii << ")";
                }
                orig_rows_PCG.push_back(ii);

                if (debug_) { std::cout << "=> orig_rows_PCG is currently: ";
                    for (int kk = 0; kk < orig_rows_PCG.size(); ++kk) {
                        std::cout << orig_rows_PCG[kk] << " "; };
                    std::cout << std::endl;
                }
                break;
            }
        }
    }
    // PCG: GET INDICES OF ROWS THAT HAVE BEEN REMOVED (IF LARGER) /ADDED (IF SMALLER)
    for (int jj = 0; jj < WIDataPCG.IJK.rows(); ++jj) {
        // IF jj IS NOT FOUND IS orig VECTOR
        if (std::find(orig_rows_PCG.begin(), orig_rows_PCG.end(), jj) == orig_rows_PCG.end() ){
            added_rows_PCG.push_back(jj);
            if (debug_) { std::cout << "added_rows_PCG is currently: ";
                for (int kk = 0; kk < added_rows_PCG.size(); ++kk) {
                    std::cout << added_rows_PCG[kk] << " ";}
                std::cout << std::endl;
            }
        }
        else {
//             if (debug_) { std::cout << "no indices added to added_rows_PCG" << std::endl; };
        }
    }

    // FIND WHICH ROWS HAVE BEEN REMOVED FROM RMS DATA
    // RMS: GET INDICES OF ROWS THAT ARE IN ORIGINAL IJK DATA
    for (int jj = 0; jj < IJK_RMS_NEW.rows(); ++jj) {

        for (int ii = 0; ii < WIDataRMS.IJK.rows(); ++ii) {
            // FIND WHICH ORIG ROW FROM RMS DATA THAT ARE NOW IN NEW IJK MATRIX
            if (IJK_RMS_NEW.row(jj).cwiseEqual(WIDataRMS.IJK.row(ii)).count() == 4){

                if (debug_) { std::cout << "new RMS IJK row (jj=" << jj
                                        << ") == old RMS IJK row (ii=" << ii << ")";
                }
                orig_rows_RMS.push_back(ii);

                if (debug_) { std::cout << "=> orig_rows_RMS is currently: ";
                    for (int kk = 0; kk < orig_rows_RMS.size(); ++kk) {
                        std::cout << orig_rows_RMS[kk] << " "; };
                    std::cout << std::endl;
                }
                break;
            }
        }
    }
    // RMS: GET INDICES OF ROWS THAT HAVE BEEN REMOVED (IF LARGER) /ADDED (IF SMALLER)
    for (int jj = 0; jj < WIDataRMS.IJK.rows(); ++jj) {
        // IF jj IS NOT FOUND IS orig VECTOR
        if (std::find(orig_rows_RMS.begin(), orig_rows_RMS.end(), jj) == orig_rows_RMS.end() ){
            added_rows_RMS.push_back(jj);
            if (debug_) { std::cout << "added_rows_RMS is currently: ";
                for (int kk = 0; kk < added_rows_RMS.size(); ++kk) {
                    std::cout << added_rows_RMS[kk] << " ";}
                std::cout << std::endl;
            }
        }
        else {
//             if (debug_) { std::cout << "no indices added to added_rows_PCG" << std::endl; };
        }
    }
    // ===============================================================================

    // FOR CHECK: TRANSFORM STD PCG AND RMS RELATIVE INDICES TO EIGEN FORMAT FOR TESTING
    Matrix<int, Dynamic, 1> orig_rows_PCG_test, orig_rows_RMS_test;
    orig_rows_PCG_test = Matrix<int, Dynamic, 1>::Map(&orig_rows_PCG[0], orig_rows_PCG.size());
    orig_rows_RMS_test = Matrix<int, Dynamic, 1>::Map(&orig_rows_RMS[0], orig_rows_RMS.size());

    // std::cout << IJK_IDX_UNION.IJK_PCG_IDX.cwiseEqual(orig_rows_PCG_test);
    // std::cout << IJK_IDX_UNION.IJK_PCG_IDX.cwiseEqual(orig_rows_RMS_test);
    // NOTE: B/C OF UNION, IJK_IDX_UNION.IJK_PCG_IDX = IJK_IDX_UNION.IJK_RMS_IDX, *BOTH*
    // RELATIVE TO INDEX BASE OF THE LARGEST DATASET

    // THIS IS ONLY TRUE FOR ONE-WAY TRANSFER...
    // THE DIFF INDEX BASE IS DETERMINED BY THE DATASET WITH THE LARGEST NUMBER OF ROWS
    // THIS IS A COMPARISON B/E THE UNION INDEX AGAINST THE INDEX RELATIVE TO PCG AND
    // THE INDEX RELATIVE TO RMS (OF THESE IS THE LARGEST, AND THIS HAS TO MATCH THE
    // DIFF INDEX SET
    // assert(
    //    (IJK_IDX_UNION.IJK_PCG_IDX.cwiseEqual(orig_rows_PCG_test).count()==orig_rows_PCG.size()) ||
    //        (IJK_IDX_UNION.IJK_PCG_IDX.cwiseEqual(orig_rows_RMS_test).count()==orig_rows_RMS.size())
    // );

    // UPDATE PCG/RMS IJK DATA
    if (debug_) {
        std::cout << "old: WIDataPCG.IJK.rows():" << WIDataPCG.IJK.rows() << std::endl;
        std::cout << "old: WIDataRMS.IJK.rows():" << WIDataRMS.IJK.rows() << std::endl;
        std::cout << "updating PCG/RMS IJK data: resizing old IJK vectors... ";
    }

    WIDataPCG.IJK.resize(IJK_PCG_NEW.rows(),4);
    WIDataRMS.IJK.resize(IJK_RMS_NEW.rows(),4);
    WIDataPCG.IJK << IJK_PCG_NEW;
    WIDataRMS.IJK << IJK_RMS_NEW;

    if (debug_) {
        std::cout << "ok." << std::endl;
        std::cout << "new: WIDataPCG.IJK.rows():" << WIDataPCG.IJK.rows() << std::endl;
        std::cout << "new: WIDataPCG.IJK.rows():" << WIDataPCG.IJK.rows() << std::endl << std::endl;

        // SETTING UP TEMPORARY WCF CONTAINER
        std::cout << "old: WIDataPCG.WCF.rows():" << WIDataPCG.WCF.rows() << std::endl;
        std::cout << "old: WIDataRMS.WCF.rows():" << WIDataRMS.WCF.rows() << std::endl;
    }

    Matrix<double, Dynamic, 1> WCF_PCG_TEMP, WCF_RMS_TEMP;
    WCF_PCG_TEMP.resize(orig_rows_PCG.size(),1);
    WCF_RMS_TEMP.resize(orig_rows_RMS.size(),1);
    WCF_PCG_TEMP.fill(0);
    WCF_RMS_TEMP.fill(0);

    if (debug_) {
        // SELECT WCF ROWS FROM ORIGINAL PCG/RMS DATA USING OBTAINED INDICES FROM REMOVED ROWS
        std::cout << "updating PCG/RMS WCF data: filling in temp WCF vectors from old... ";
    }
    for (int jj = 0; jj < orig_rows_PCG.size(); ++jj) {
        WCF_PCG_TEMP.row(jj) = WIDataPCG.WCF.row(orig_rows_PCG[jj]);
    }
    for (int jj = 0; jj < orig_rows_RMS.size(); ++jj) {
        WCF_RMS_TEMP.row(jj) = WIDataRMS.WCF.row(orig_rows_RMS[jj]);
    }

    if (debug_) {
        std::cout << "ok." << std::endl;
        // UPDATE PCG/RMS WCF DATA
        std::cout << "updating PCG/RMS WCF data: resizing and adding to old WCF vectors... ";
    }
    WIDataPCG.WCF.resize(orig_rows_PCG.size(),1);
    WIDataRMS.WCF.resize(orig_rows_RMS.size(),1);
    WIDataPCG.WCF << WCF_PCG_TEMP;
    WIDataRMS.WCF << WCF_RMS_TEMP;

    if (debug_) {
        std::cout << "ok." << std::endl;
        std::cout << "new: WIDataPCG.WCF.rows():" << WIDataPCG.WCF.rows() << std::endl;
        std::cout << "new: WIDataPCG.WCF.rows():" << WIDataPCG.WCF.rows() << std::endl;
    }

    sup_indices.push_back(orig_rows_PCG);
    sup_indices.push_back(added_rows_PCG);
    sup_indices.push_back(orig_rows_RMS);
    sup_indices.push_back(added_rows_RMS);

    return sup_indices;
}

/*!
 * \brief
 * \param
 * \return
 */
void RemoveSuperfluousRowsB(WIData &WIDataRMS,
                            WIData &WIDataPCG,
                            QStringList &diff_files,
                            int jj) {

    bool debug_ = false;
    QString str_out;

    auto sup_indices = DiffTreatmentB(WIDataRMS, WIDataPCG, diff_files);

    if (debug_) {
        std::cout << "\nsup_indices: " << std::endl;
        std::cout << "orig_rows_PCG(sz=" << sup_indices[0].size() << ")" << std::endl;
        std::cout << "added_rows_PCG(sz=" << sup_indices[1].size() << ")" << std::endl;
        std::cout << "orig_rows_RMS(sz=" << sup_indices[2].size() << ")" << std::endl;
        std::cout << "added_rows_RMS(sz=" << sup_indices[3].size() << ")\n" << std::endl;
    }

    // VECTOR LENGTHS HAVE BEEN MADE EQUAL => PRINT RESULTS
    QString rem_str;
    QString ind_str;
    QStringList str_ind;

    str_out.append(">>> Vector lengths have been made equal: ");
    if (sup_indices[1].size() > 0) { // added_rows_PCG
        ind_str = (sup_indices[1].size() > 1) ?
                  QString::number(sup_indices[1].size()) + " rows were" : "1 row was";
        str_out.append(ind_str + " removed from PCG data\nb/c IJK values did not match. ");

        for (int ii = 0; ii < sup_indices[1].size(); ++ii) {
            str_ind.append(QString::number(sup_indices[1][ii]));
        }
        str_out.append("Rows that were removed: [" + str_ind.join(" ") + "].");
    }

    if (sup_indices[3].size() > 0) { // added_rows_RMS
        ind_str = (sup_indices[3].size() > 1) ?
                  QString::number(sup_indices[3].size()) + " rows were" : "1 row was";
        str_out.append(ind_str + " removed from RMS data\nb/c IJK values did not match. ");

        for (int ii = 0; ii < sup_indices[3].size(); ++ii) {
            str_ind.append(QString::number(sup_indices[3][ii]));
        }
        str_out.append("Rows that were removed: [" + str_ind.join(" ") + "].");
    }

    if (sup_indices[1].size() + sup_indices[3].size() > WIDataPCG.max_sup){
        str_out.append("\nWARNING: more than " + QString::number(WIDataPCG.max_sup)
                           + " rows removed, check wells are supposed to be equal. ");
    }

    WIDataRMS.test_IJK_removed[jj] = sup_indices;
    WIDataPCG.test_IJK_removed[jj] = sup_indices;

    str_out.append("\nContinuing comparison.");
    std::cout << "\033[1;36m" << str_out.toStdString() << "\033[0m" << std::endl;
    Utilities::FileHandling::WriteLineToFile(str_out, WIDataPCG.tex_file);
}

/*!
 * \brief
 * \param
 * \return
 */
//void RemoveSuperfluousRowsA(WIData &WIDataRMS,
//                            WIData &WIDataPCG,
//                            QStringList &diff_files,
//                            int jj){
//
//    bool debug_ = false;
//    bool remove_sup = true;
//    QVector<int> sup_indices_total;
//    QString str_out;
//
//    // DIFF TREATMENT: IJK COMPARISON: FIND EXTRA ROWS USING diff COMMAND
//    auto sup_indices = DiffTreatmentA(WIDataRMS, WIDataPCG, diff_files);
//
//    // REMOVE SUPERFLUOUS ROWS
//    WIData WILong = GetLongestVector(WIDataRMS, WIDataPCG);
//    WIData WIShort = GetShortestVector(WIDataRMS, WIDataPCG);
//    WIData WITemp;
//
//    // WIShort.IJK.setZero();
//    WITemp.IJK.resize(WILong.IJK.rows() - sup_indices.size(),4);
//    WITemp.WCF.resize(WILong.IJK.rows() - sup_indices.size(),1);
//    WITemp.IJK.fill(0);
//    WITemp.WCF.fill(0);
//
//    if (debug_){
//        std::cout << "WILong.IJK.rows:" << WILong.IJK.rows() << std::endl;
//        std::cout << "WIShort.IJK.rows:" << WIShort.IJK.rows() << std::endl;
//        std::cout << "WITemp.IJK.rows:" << WITemp.IJK.rows() << std::endl;
//    }
//
//    // LOOP OVER ALL ROWS IN THE LONGEST COLUMN AND INSERT EACH OF THESE INTO A NEW IJK
//    // COLUMN UNLESS THE GIVEN ROW IS A SUPERFLUOUS ONE, IN WHICH CASE WE SKIP IT
//    int kk = 0;
//    for (int ii = 0; ii < WILong.IJK.rows(); ++ii) {
//        if (debug_) std::cout << "ii: " <<  ii
//                              << "[kk: " << kk << "] "
//                              << "is current row superfluous? (1=yes, 0=no): "
//                              << sup_indices.contains(ii) << std::endl;
//        if (! sup_indices.contains(ii)) {
//            WITemp.IJK.row(kk) << WILong.IJK.row(ii);
//            WITemp.WCF.row(kk) << WILong.WCF.row(ii);
//            kk += 1;
//        }else{
//            if (debug_) std::cout << "Row not added to short vector!" << std::endl;
//        }
//    }
//
//    // REMOVE SUPERFLUOUS ROWS FROM ORIGINALLY LONG COLUMN (UPDATE LONG COLUMN)
//    QString rem_str;
//    if (WIDataRMS.IJK.rows() > WIDataPCG.IJK.rows()) {
//        WIDataRMS.IJK = WITemp.IJK;
//        WIDataRMS.WCF = WITemp.WCF;
//        rem_str = "RMS";
//    } else {
//        WIDataPCG.IJK = WITemp.IJK;
//        WIDataPCG.WCF = WITemp.WCF;
//        rem_str = "PCG";
//    }
//
//    if (debug_){
//        std::cout << "WIDataRMS.IJK.rows:" << WIDataRMS.IJK.rows() << std::endl;
//        std::cout << "WIDataPCG.IJK.rows:" << WIDataPCG.IJK.rows() << std::endl;
//        std::cout << "sup_indices.size():" << sup_indices.size() << std::endl;
//    }
//
//    // VECTOR LENGTHS HAVE BEEN MADE EQUAL => COMPARE DIRECTLY
//    QString ind_str = (sup_indices.length() > 1) ?
//                      QString::number(sup_indices.size()) + " rows were" : "1 row was";
//    str_out.append(">>> Vector lengths have been made equal: "
//                       + ind_str + " removed from "
//                       + rem_str + " data\nb/c IJK values did not match. ");
//
//    QStringList str_ind;
//        foreach(int ii, sup_indices){ str_ind.append(QString::number(ii)); }
//    str_out.append("Rows that were removed: [" + str_ind.join(" ") + "].");
//
//    if (sup_indices_total.length()>5){
//        str_out.append("\nWARNING: more than 5 rows removed, "
//                           "check wells are supposed to be equal. ");
//    }
//
//    str_out.append("\nContinuing comparison.");
//    std::cout << "\033[1;36m" << str_out.toStdString() << "\033[0m" << std::endl;
//    Utilities::FileHandling::WriteLineToFile(str_out, WIDataPCG.tex_file);
//}


/*!
 * \brief
 * \param
 * \return
 */
void RemoveSuperfluousRowsWrapper(WIData &WIDataRMS,
                                  WIData &WIDataPCG,
                                  QStringList &dir_list_,
                                  QStringList &dir_names_,
                                  int ii){

//    WIDataRMS.test_IJK_removed.clear();
//    WIDataPCG.test_IJK_removed.clear();

//    for (int ii = 0; ii < 4; ++ii) {
//        WIDataRMS.test_IJK_removed[ii][0] = 0;
//        WIDataPCG.test_IJK_removed[ii][0] = 0;
//    }

    // PRINT IJK, WCF DATA TO INDIVIDUAL FILES TO TREAT WITH diff COMMAND LATER:
    // MAKE DIFF FILE NAMES
    QStringList diff_files = {
        dir_list_[ii] + "/DIFF_" + dir_names_[ii] + "_RMS.IJK",
        dir_list_[ii] + "/DIFF_" + dir_names_[ii] + "_PCG.IJK",
        dir_list_[ii] + "/DIFF_" + dir_names_[ii] + "_RMS.WCF",
        dir_list_[ii] + "/DIFF_" + dir_names_[ii] + "_PCG.WCF"
    };

    // PRINT DIFF FILES
    WIDataRMS.PrintIJKData(diff_files[0]);
    WIDataPCG.PrintIJKData(diff_files[1]);
    WIDataRMS.PrintWCFData(diff_files[2]);
    WIDataPCG.PrintWCFData(diff_files[3]);

    QString str_out;
    QString lstr_out = "\n----------------------------------------"
        "----------------------------------------";
    if (DiffVectorLength(WIDataRMS, WIDataPCG)) {

        // IF VECTOR LENGTHS ARE EQUAL => COMPARE DIRECTLY
        str_out = lstr_out + "\n"
            + ">>> COMPDAT vectors have the same length. Making comparison.";
        std::cout << std::endl << "\033[1;36m"
                  << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, WIDataPCG.tex_file);

        std::vector<int> empty_v;
        std::vector< std::vector<int> > emp_indices;

        emp_indices.push_back(empty_v);
        emp_indices.push_back(empty_v);
        emp_indices.push_back(empty_v);
        emp_indices.push_back(empty_v);

        WIDataRMS.test_IJK_removed[ii] = emp_indices;
        WIDataPCG.test_IJK_removed[ii] = emp_indices;

    } else {

        // IF VECTOR LENGTHS ARE UNEQUAL => MAKE EQUAL, THEN COMPARE DIRECTLY
        str_out = lstr_out + "\n"
            + ">>> COMPDAT vectors have different length. Making them equal.";
        std::cout << std::endl << "\033[1;36m"
                  << str_out.toStdString() << "\033[0m" << std::endl;
        Utilities::FileHandling::WriteLineToFile(str_out, WIDataPCG.tex_file);

        RemoveSuperfluousRowsB(WIDataRMS, WIDataPCG, diff_files, ii);
    }
}

// ///////////////////////////////////////////////////////////////////////////////
// DEBUG MESSAGES
template<typename ZA, typename ZB, typename ZC, typename ZD, typename ZE, typename ZF>
void debug_msg(bool debug_, string msg,
               ZA varA, // dir_names_
               ZB varB, // dir_list_
               ZC varC, // ii
               ZD varD, // WIDataRMS
               ZE varE, // WIDataPCG
               ZF varF  // empty
){

    if (msg.compare("well_dir_list")==0){
        if (debug_) {
            std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl;
            for (int ii = 0; ii < varA.length(); ++ii) {
                std::cout << ii << ":" << varA[ii].toStdString() << std::endl; // list of well dirs (names only)
                std::cout << ii << ":" << varB[ii].toStdString() << std::endl; // list of well dirs (absolute path)
            }
            std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
        }
    }else if(msg.compare("RMS_PCG_IJK_data")==0){
        if (debug_) {
            int nRMS = (varD.IJK.rows() > 5) ? 5 : varD.IJK.rows();
            int nPCG = (varE.IJK.rows() > 5) ? 5 : varE.IJK.rows();

            std::cout << "\033[1;31m<DEBUG:START->\033[0m" << std::endl << std::setfill(' ');

            // RMS-PCG: IJK
            std::cout << "RMS-PCG IJK DATA (well = " << varA[varC].toStdString() << ")" << std::endl;
            // RMS: IJK
            std::cout << "WIDataRMS.IJK (size: " << varD.IJK.size() << "): "
                      << std::endl << varD.IJK.block(0, 0, nRMS, 4)
                      << std::endl << "..." << std::endl;
            // PCG: IJK
            std::cout << "WIDataPCG.IJK (size: " << varE.IJK.size() << "): "
                      << std::endl << varE.IJK.block(0, 0, nPCG, 4)
                      << std::endl << "..." << std::endl;

            // RMS-PCG: WCF
            std::cout << "RMS-PCG WCF DATA (well = " << varA[varC].toStdString() << ")" << std::endl;
            // RMS: WCF
            std::cout << "WIDataRMS.WCF: (size: " << varD.WCF.size() << "): "
                      << std::endl << varD.WCF.block(0, 0, nRMS, 1)
                      << std::endl << "..." << std::endl;
            // PCG: WCF
            std::cout << "WIDataPCG.WCF: (size: " << varE.WCF.size() << "): "
                      << std::endl << varE.WCF.block(0, 0, nPCG, 1)
                      << std::endl << "..." << std::endl;

            std::cout << "\033[1;31m<DEBUG:END--->\033[0m" << std::endl;
        }
    }else if(msg.compare("C")==0){

    }else if(msg.compare("D")==0){

    }else if(msg.compare("E")==0){

    }

}

}
}
#endif //FIELDOPT_TEST_RESOURCE_DIFF_FUNCTIONS_H