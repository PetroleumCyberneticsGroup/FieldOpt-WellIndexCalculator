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

#include <gtest/gtest.h>
#include "Reservoir/grid/grid.h"
#include "Reservoir/grid/eclgrid.h"
#include "Utilities/filehandling.hpp"
#include <FieldOpt-WellIndexCalculator/wellindexcalculator.h>
#include <FieldOpt-WellIndexCalculator/main.hpp>
#include "test_resource_wic_welldir.h"
#include "test_resource_wic_diff_functions.h"

using namespace Reservoir::Grid;
//using namespace WellIndexCalculation;
using namespace TestResources::WIBenchmark;

namespace {

class DeviatedWellIndexTest : public ::testing::Test {

 protected:
  DeviatedWellIndexTest() {
      well_dir_ = new WellDir();
  }

  virtual ~DeviatedWellIndexTest() {
  }

  WellDir *well_dir_;
  QString grid_file_5spot = "../examples/ADGPRS/5spot/ECL_5SPOT.EGRID";
  QString grid_file_norne = "../examples/Flow/norne/OUTPUT/NORNE_ATW2013.EGRID";
  QString tex_smry_5spot, tex_smry_norne;
  QString tab_header_str, tab_tail_str, tab_mean_str;

  bool debug_ = false;

  struct smry_data{
    vector<double> well_mean, well_median;
  };

  smry_data smry_data_;

};

TEST_F(DeviatedWellIndexTest, compareCOMPDAT) {

    // GET LIST OF WELL FOLDERS CONTAINING PCG &
    // RMS COMPDATS (OBTAINED USING WI_BENCHMARK CODE)
    auto file_list_ = well_dir_->GetWellDir();

    // list of well dirs (names only) => dir name only: tw04_04
    auto dir_names_ = file_list_[0];

    // list of well dirs (absolute path) => fullpath: ../tw04_04/
    auto dir_list_ = file_list_[1];

    // fullpath: ../tw04_04/EVENTS_tw04_04_RMS.DATA
    auto rms_files = file_list_[2];

    // fullpath: ../tw04_04/EVENTS_tw04_04_PCG.DATA
    auto pcg_files = file_list_[3];

    // WELL INDEX DATA OBJECTS
    WIData WIDataRMS, WIDataPCG;
    WIDataRMS.data_tag = "RMS";
    WIDataPCG.data_tag = "PCG";

    tex_smry_5spot = dir_list_[0] + "/../summary-table-5spot.tex";
    tex_smry_norne = dir_list_[0] + "/../summary-table-norne.tex";

    tab_header_str = "\\begin{table} \\begin{tabular}{c c c}\n"
        "Well & "
        "{\\bfseries $mean^{th}$} & "
        "{\\bfseries $median$}\\\\"
        "\\toprule";
    tab_tail_str = "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\caption{"
        "$mean^{th}$: the mean of the $\\dfrac{wi_{RMS}}{wi_{PCG}}$ vector for each well,"
        "with outliers removed (above and below threshold values $[\\frac{1/2} 2]$).\n"
        "$median$: the mediam of the $\\\\dfrac{wi_{RMS}}{wi_{PCG}}$ vector for each well,"
        "no values removed.}\n"
        "\\end{table}";

    Utilities::FileHandling::WriteStringToFile(tab_header_str, tex_smry_5spot);
    Utilities::FileHandling::WriteStringToFile(tab_header_str, tex_smry_norne);

    WIDataRMS.test_IJK_removed.resize(rms_files.length());
    WIDataPCG.test_IJK_removed.resize(rms_files.length());

    // DEBUG
    debug_msg(false, "well_dir_list", dir_names_,
              dir_list_, 0, WIDataRMS, WIDataPCG, 0);

    // LOOP THROUGH LIST OF WELL FOLDERS: FOR WELL
    // FOLDER ii: READ PCG & RMS COMPDAT DATA
    int num_files = (debug_) ? 5 : rms_files.length(); //override
    QString str_out;
    QString lstr_out = "================================================================================";

//    WIDataRMS.test_IJK_removed.resize(num_files,1);
//    WIDataPCG.test_IJK_removed.resize(num_files,1);
//    WIDataRMS.test_IJK_removed.fill(0);
//    WIDataPCG.test_IJK_removed.fill(0);

    for (int ii = 0; ii < num_files; ++ii) {

        // USE COMPDAT DATA (RMS) PRODUCED BY BENCHMARK PROGRAM
        WIDataRMS.ReadCOMPDAT(rms_files[ii], dir_list_[ii]);

        // UNCOMMENT IF READING PCG COMPDAT TABLE FROM BENCHMARK PROGRAM
        // WIDataPCG.ReadCOMPDAT(pcg_files[ii],dir_list_[ii]);
        // DEBUG: PRINT OLD IJK/WCF VALUES B/F SHIFTING TO NEW
        // WIDataPCG.PrintIJKData(dir_list_[ii] + "/DBG_" + dir_names_[ii] + "_PCG.IJK");
        // WIDataPCG.PrintWCFData(dir_list_[ii] + "/DBG_" + dir_names_[ii] + "_PCG.WCF");

        // DECIDE ON WHICH GRID FILE TO USE (5SPOT OR NORNE)
        if( QString::compare(dir_names_[ii].at(0), "n", Qt::CaseSensitive) == 0 ) {
            WIDataPCG.grid_file = grid_file_norne;
            WIDataPCG.well_name = "NW01";
            WIDataPCG.tex_smry = tex_smry_norne;
        }
        else if( QString::compare(dir_names_[ii].at(0), "t", Qt::CaseSensitive) == 0 ) {
            WIDataPCG.grid_file = grid_file_5spot;
            WIDataPCG.well_name = "TW01";
            WIDataPCG.tex_smry = tex_smry_5spot;
        }
        WIDataPCG.dir_name = dir_names_[ii];
        WIDataRMS.tex_smry = WIDataPCG.tex_smry;

        // DEFINE TEX FILES + START WRITING TO FILE
        WIDataPCG.tex_file = dir_list_[ii] + "/" + dir_names_[ii] + ".tex";
        WIDataRMS.tex_file = WIDataPCG.tex_file;
        Utilities::FileHandling::WriteStringToFile("\\begin{alltt}", WIDataPCG.tex_file);

        // MAKE NEW COMPDAT DATA USING PRECOMPILED WellIndexCalculator
        WIDataRMS.ReadXYZ(dir_list_[ii] + "/" + dir_names_[ii]);
        WIDataPCG.ReadXYZ(dir_list_[ii] + "/" + dir_names_[ii]);
        WIDataPCG.CalculateWCF(dir_list_[ii] + "/" + dir_names_[ii]);

        if( QString::compare(dir_names_[ii].at(0), "t", Qt::CaseSensitive) == 0 ) {
            WIDataPCG.PrintCOMPDATPlot(dir_list_[ii] + "/" + dir_names_[ii]);
        }

        // USE DATA COMPUTED USING WellIndexCalculator INSTEAD OF OLD DATA
        WIDataPCG.IJK = WIDataPCG.IJKN;
        WIDataPCG.WCF = WIDataPCG.WCFN;

        // WRITE TO TEX FILE
        str_out = "\n" + lstr_out + "\nChecking IJK and WCF data for well: "
            + dir_names_[ii] + " (row numbering uses 0-indexing)\nheel.xyz = "
            + WIDataPCG.XYZh + ", toe.xyz = " + WIDataPCG.XYZt
            + "\n" + lstr_out;
        Utilities::FileHandling::WriteLineToFile(str_out, WIDataPCG.tex_file);
        std::cout << "\033[1;36m" << str_out.toStdString() << "\033[0m";

        // DEBUG
        debug_msg(false, "RMS_PCG_IJK_data",
                  dir_names_, dir_list_, ii,
                  WIDataRMS, WIDataPCG, 0);

        // REMOVE ROW IF LOW WCF
        RemoveRowsLowWCF(WIDataRMS);
        RemoveRowsLowWCF(WIDataPCG);

        // REMOVE EXTRA ROWS IF DATA HAS UNEQUAL LENGTH
        RemoveSuperfluousRowsWrapper(WIDataRMS, WIDataPCG,
                                     dir_list_, dir_names_,
                                     ii);

        // COMPARE IJK AND PCG VALUES (EQUAL LENGTH DATA)
        CompareIJK(WIDataRMS, WIDataPCG);
        auto WIDiff = CompareWCF(WIDataRMS, WIDataPCG, ii);

        // WRITE TO TEX FILE
        Utilities::FileHandling::WriteLineToFile("\\end{alltt}", WIDataPCG.tex_file);

        // COLLECT DATA FOR OVERALL MEAN
        smry_data_.well_mean.push_back(WIDiff.WCF_accuracy_list[4]);
        smry_data_.well_median.push_back(WIDiff.WCF_accuracy_list[6]);
    }

    tab_mean_str = "\\midrule\\\\\nMean & "
        + QString::number(ConvertStdToEigen(
            smry_data_.well_mean).mean()).leftJustified(5, '0', true)
        + " & "
        + QString::number(ConvertStdToEigen(
            smry_data_.well_median).mean()).leftJustified(5, '0', true);

    Utilities::FileHandling::WriteLineToFile(tab_tail_str, tex_smry_5spot);
    Utilities::FileHandling::WriteLineToFile(tab_tail_str, tex_smry_norne);

//    std::cout << "\n\n*********\nTESTING\n*********\n" << std::endl;
//    for (int ii = 0; ii < num_files; ++ii) {
//        std::cout << "Well name: " << dir_names_[ii].toStdString()
//                  << " -- WIDataRMS.test_IJK_removed: "
//                  << WIDataRMS.test_IJK_removed.rows()
//                  << " (n_rows); ii=" << ii << "; "
//                  << WIDataRMS.test_IJK_removed(ii)
//                  << " (n_removed)" << std::endl;
//
//        EXPECT_LT(WIDataRMS.test_IJK_removed(ii), 10);
//        EXPECT_LT(WIDataPCG.test_IJK_removed(ii), 10);
//    }
}
}