cmake_minimum_required(VERSION 2.8)
project(wellindexcalculator)

# Dbg ==================================================================
message("=============================================================")
message(STATUS "PROJECT_NAME => ${PROJECT_NAME}")
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(WIC_LIB_TARGET wellindexcalculator)
set(WIC_LIB_TARGET_ADG wellindexcalculator-adgprs-dbg)

# ResInsight ===========================================================

set(CVF_ENABLE_ASSERTS 1)
set(CVF_ENABLE_TIGHT_ASSERTS 1)

#if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#  add_definitions(-DCVF_LINUX)
#  message(STATUS "CMAKE_SYSTEM_NAME => ${CMAKE_SYSTEM_NAME}")
#endif()
#
#if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DCVF_LINUX -pipe -Wextra ")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#  set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -g3 -O0 -DDEBUG -D_DEBUG")
#  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNO_DEBUG")
#endif()

#include(ResInsight.cmake)

# Include: Eigen =======================================================
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")

# Include: Qt ==========================================================
include_directories(${Qt5Core_INCLUDE_DIRS})
message(STATUS "Qt5Core_INCLUDE_DIRS: ${Qt5Core_INCLUDE_DIRS}")

# Include: ERT =========================================================

#find_library(ert_ecl
#    NAMES ecl
#    HINTS ${CMAKE_BINARY_DIR}/libraries/lib
#    REQUIRED)

find_library(ert_ecl
        NAMES ecl
        HINTS ${PROJECT_SOURCE_DIR}/resinsight/Ert/cmake-build-debug/lib
        REQUIRED)

#find_library(ert_eclxx
#    NAMES ecl
#    HINTS ${CMAKE_BINARY_DIR}/libraries/lib
#    REQUIRED)

find_library(ert_util
        NAMES ert_util
        HINTS ${CMAKE_BINARY_DIR}/libraries/lib
        REQUIRED)

#find_library(ert_utilxx
#    NAMES ert_util
#    HINTS ${CMAKE_BINARY_DIR}/libraries/lib
#    REQUIRED)

find_library(ert_well
        NAMES ecl_well
        HINTS ${CMAKE_BINARY_DIR}/libraries/lib
        REQUIRED)

find_library(ert_geometry
        NAMES ert_geometry
        HINTS ${CMAKE_BINARY_DIR}/libraries/lib
        REQUIRED)

set(ERT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/resinsight/Ert/lib/include)
message(STATUS "ERT_INCLUDE_DIR: ${ERT_INCLUDE_DIR}")
#include_directories(${ERT_INCLUDE_DIR})

# wellindexcalculator lib ==============================================
if (_WIN32) # WINDOWS
    find_library(SHLWAPI_LIBRARY NAMES Shlwapi)

    include_directories(
            ${EIGEN3_INCLUDE_DIR}
            ${ERT_INCLUDE_DIR})

    add_library(${WIC_LIB_TARGET} SHARED dllmain.cpp)

    target_link_libraries(${WIC_LIB_TARGET}
            fieldopt::${WIC_LIB_TARGET}
            ${Boost_LIBRARIES}
            ${SHLWAPI_LIBRARY})

    if(BUILD_WIC_ADGPRS)
        find_library(SHLWAPI_LIBRARY NAMES Shlwapi)

        target_link_libraries(wicalc-adgprs
                fieldopt::wellindexcalculator
                ${Boost_LIBRARIES}
                ${SHLWAPI_LIBRARY})
    endif(BUILD_WIC_ADGPRS)

else() # LINUX --------------------------------------------------------
    if (NOT BUILD_WIC_ADGPRS)

        #    include_directories(
        #        ${RESINSIGHT_FWK_VIZFWK_LIBCORE}
        #        ${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}
        #        ${RESINSIGHT_FWK_VIZFWK_LIBRENDER}
        #        ${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL}
        #        ${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE}
        #        ${RESINSIGHT_FWK_APPFWK_COMMONCODE}
        #        ${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}
        #        ${RESINSIGHT_APPLICATIONCODE_APPLICATION}
        #        ${RESINSIGHT_APPLICATIONCODE_PROJECTDATAMODEL}
        #        ${RESINSIGHT_APPLICATIONCODE_APPLICATION_TOOLS}
        #        ${RESINSIGHT_APPLICATIONCODE_FILEINTERFACE}
        #        )

        include_directories(
                resinxx/rixx_app_fwk
                resinxx/rixx_core_geom
                resinxx/rixx_res_mod
                resinxx/rixx_grid
        )

        include_directories(
                ${EIGEN3_INCLUDE_DIR}
                ${ERT_INCLUDE_DIR})

        add_library(${WIC_LIB_TARGET}
                intersected_cell.cpp
                wellindexcalculator.cpp
                #        wicalc_rins.cpp
                #        ${RESINSIGHT_SOURCE}
                #        wicalc_rinx.cpp
                //
                wicalc_rixx.cpp
                resinxx/rixx_grid/riextractor.cpp
                resinxx/rixx_grid/rifaultncc.cpp
                resinxx/rixx_grid/ricasedata.cpp
                resinxx/rixx_grid/rigrid.cpp
                resinxx/rixx_grid/ricell.cpp
                //
                resinxx/well_path.cpp
                resinxx/geometry_tools.cpp
                //
                resinxx/rixx_app_fwk/cvfStructGrid.cpp
                resinxx/rixx_app_fwk/cvfCellRange.cpp
                resinxx/rixx_app_fwk/cafHexGridIntersectionTools.cpp
                //
                resinxx/rixx_core_geom/cvfAssert.cpp
                resinxx/rixx_core_geom/cvfAtomicCounter.cpp
                resinxx/rixx_core_geom/cvfBoundingBox.cpp
                resinxx/rixx_core_geom/cvfBoundingBoxTree.cpp
                resinxx/rixx_core_geom/cvfCharArray.cpp
                resinxx/rixx_core_geom/cvfMath.cpp
                resinxx/rixx_core_geom/cvfPlane.cpp
                resinxx/rixx_core_geom/cvfObject.cpp
                resinxx/rixx_core_geom/cvfRay.cpp
                resinxx/rixx_core_geom/cvfString.cpp
                resinxx/rixx_core_geom/cvfSystem.cpp
                resinxx/rixx_core_geom/cvfVector2.cpp
                resinxx/rixx_core_geom/cvfVector3.cpp
                resinxx/rixx_core_geom/cvfVector4.cpp
                //
                resinxx/rixx_res_mod/cvfGeometryTools.cpp
                resinxx/rixx_res_mod/RigCellGeometryTools.cpp
                //
                tests/wic_debug.hpp
                )

        add_library(fieldopt::${WIC_LIB_TARGET}
                ALIAS ${PROJECT_NAME})

        message(STATUS "RESINSIGHT_SOURCE: ${RESINSIGHT_SOURCE}")

    elseif(BUILD_WIC_ADGPRS)
        add_library(${WIC_LIB_TARGET_ADG} SHARED
                dllmain.cpp
                intersected_cell.cpp
                wellindexcalculator.cpp
                wicalc_rixx.cpp
                tests/wic_debug.hpp
                )
    endif()

    if (NOT BUILD_WIC_ADGPRS)

        include_directories(${CMAKE_BINARY_DIR}/libraries/include/)

        add_compile_options(-std=c++11)
        set_target_properties(${WIC_LIB_TARGET}
                PROPERTIES COMPILE_FLAGS "-fPIC")

        target_link_libraries(${WIC_LIB_TARGET}
                PUBLIC fieldopt::${WIC_LIB_TARGET}
                PUBLIC fieldopt::reservoir
                PUBLIC ${Boost_LIBRARIES}
                PUBLIC Qt5::Core
                #        PUBLIC ${ert_ecl}
                #        PUBLIC ${ert_eclxx}
                PUBLIC ${ert_util}
                #        PUBLIC ${ert_utilxx}
                PUBLIC ${ert_well}
                PUBLIC ${ert_geometry}
                PUBLIC ${PROJECT_SOURCE_DIR}/resinsight/Ert/cmake-build-debug/lib/libecl.so
                )

        install(TARGETS wellindexcalculator
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib/static
                )

        #    add_executable(wicalc
        #        main.cpp
        #        intersected_cell.cpp
        #        wellindexcalculator.cpp
        #        wicalc_rins.cpp
        #        ${RESINSIGHT_SOURCE}
        #        )
        #
        #    target_link_libraries(wicalc
        #        PUBLIC fieldopt::${WIC_LIB_TARGET}
        #        PUBLIC ${Boost_LIBRARIES}
        #        PUBLIC Qt5::Core
        #        PUBLIC ${ert_ecl}
        #        PUBLIC ${ert_util}
        #        PUBLIC ${RESINSIGHT_HEADERS}
        #        )

        #    install(TARGETS wicalc wellindexcalculator
        #        RUNTIME DESTINATION bin
        #        LIBRARY DESTINATION lib
        #        ARCHIVE DESTINATION lib/static
        #        )

    elseif(BUILD_WIC_ADGPRS)
        add_library(fieldopt::${WIC_LIB_TARGET_ADG}
                ALIAS ${PROJECT_NAME}-adgprs)

        target_link_libraries (${WIC_LIB_TARGET_ADG}
                PUBLIC fieldopt::reservoir
                PUBLIC Qt5::Core
                ${Boost_LIBRARIES})

        add_compile_options(-std=c++11)

        add_executable(wicalc-adgprs dllmain.cpp)

        target_link_libraries(wicalc-adgprs
                fieldopt::${WIC_LIB_TARGET_ADG}
                PUBLIC Qt5::Core
                ${Boost_LIBRARIES}
                )
    endif()
endif()


# Testing ==============================================================
if (BUILD_TESTING AND NOT BUILD_WIC_ADGPRS)

    # Unit tests
    find_package(GTest REQUIRED)

    include_directories(${GTEST_INCLUDE_DIRS}
            ${EIGEN3_INCLUDE_DIR} tests)

    add_executable(test_wellindexcalculator
            tests/test_intersected_cells.cpp
            tests/test_single_cell_wellindex.cpp
            tests/test_standalone_wicalc.cpp)

    target_link_libraries(test_wellindexcalculator
            fieldopt::wellindexcalculator
            ${GTEST_BOTH_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

    add_test(NAME test_wellindexcalculator
            COMMAND $<TARGET_FILE:test_wellindexcalculator>)
endif()


# WI benckmark =========================================================
if (BUILD_BENCHMARK)
    add_executable(benchmark_wellindexcalculator
            tests/well_index_benchmark/test_resource_wic_welldir.h
            tests/well_index_benchmark/test_resource_wic_widata.h
            tests/well_index_benchmark/test_well_index_deviated_well.cpp
            tests/well_index_benchmark/test_resource_wic_diff_functions.h
            resinxx/well_path.cpp
            resinxx/well_path.h
            resinxx/geometry_tools.cpp
            resinxx/geometry_tools.h)

    target_link_libraries(benchmark_wellindexcalculator
            fieldopt::wellindexcalculator
            Qt5::Core
            ${GTEST_BOTH_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT}
            ${GTEST_BOTH_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT}
            ${Boost_LIBRARIES})
endif()

# Debug: include dirs ==================================================
message(".............................................................")
message(STATUS "List current include dirs:")
get_property(DIRS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(DIR ${DIRS})
    message(STATUS "DIR='${DIR}'")
endforeach()