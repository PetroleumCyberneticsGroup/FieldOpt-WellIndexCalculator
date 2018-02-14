cmake_minimum_required(VERSION 2.8)
project(wellindexcalculator)

# ResInsight ===========================================================
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions(-DCVF_LINUX)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS 		"-DCVF_LINUX -pipe -Wextra -Woverloaded-virtual -Wformat -std=c++0x")
    set(CMAKE_CXX_FLAGS_DEBUG 	"-ggdb -g3 -O0 -DDEBUG -D_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNO_DEBUG")
endif()

include(ResInsight.cmake)

# Include: Eigen =======================================================
include_directories(${EIGEN3_INCLUDE_DIR})

# wellindexcalculator lib ==============================================
if (_WIN32) # WINDOWS
    find_library(SHLWAPI_LIBRARY NAMES Shlwapi)

    include_directories(${EIGEN3_INCLUDE_DIR}
            ${ERT_INCLUDE_DIR})

    add_library(wellindexcalculator SHARED dllmain.cpp)

    target_link_libraries(wellindexcalculator
            fieldopt::wellindexcalculator
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

        include_directories(
            ${RESINSIGHT_FWK_VIZFWK_LIBCORE}
            ${RESINSIGHT_FWK_VIZFWK_LIBGEOMETRY}
            ${RESINSIGHT_FWK_VIZFWK_LIBRENDER}
            ${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL}
            ${RESINSIGHT_FWK_APPFWK_CAFPROJECTDATAMODEL_CAFPDMCORE}
            ${RESINSIGHT_FWK_APPFWK_COMMONCODE}
            ${RESINSIGHT_APPLICATIONCODE_RESERVOIRDATAMODEL}
            ${RESINSIGHT_APPLICATIONCODE_APPLICATION}
            ${RESINSIGHT_APPLICATIONCODE_PROJECTDATAMODEL}
            ${RESINSIGHT_APPLICATIONCODE_APPLICATION_TOOLS}
            )

        add_library(wellindexcalculator
                intersected_cell.cpp
                wellindexcalculator.cpp
                #
                wicalc.cpp
                ${LIBCORE_SOURCE}
                ${LIBGEOMETRY_SOURCE}
                ${CAFPDMCORE_SOURCE}
                ${COMMONCODE_SOURCE}
                ${RESERVOIRDATAMODEL_SOURCE}
                ${APPLICATION_SOURCE}
                ${PROJECTDATAMODEL_SOURCE}
                ${TOOLS_SOURCE}
                )

    elseif(BUILD_WIC_ADGPRS)
        add_library(wellindexcalculator-adgprs SHARED
                wicalc.cpp
                dllmain.cpp
                intersected_cell.cpp
                wellindexcalculator.cpp
                )
    endif()


    if (NOT BUILD_WIC_ADGPRS)
        add_library(fieldopt::wellindexcalculator
                ALIAS ${PROJECT_NAME})

        target_link_libraries (wellindexcalculator
                PUBLIC fieldopt::reservoir
                ${Boost_LIBRARIES})

        include_directories(${CMAKE_BINARY_DIR}/libraries/include/)

        add_compile_options(-std=c++11 -fPIC)
               set_target_properties(wellindexcalculator
                       PROPERTIES COMPILE_FLAGS "-fPIC")

        add_executable(wicalc main.cpp)
        target_link_libraries(wicalc
                fieldopt::wellindexcalculator
                ${Boost_LIBRARIES})

        install(TARGETS wicalc wellindexcalculator
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib/static
                )

    elseif(BUILD_WIC_ADGPRS)
        add_library(fieldopt::wellindexcalculator-adgprs
                ALIAS ${PROJECT_NAME}-adgprs)

        target_link_libraries (wellindexcalculator-adgprs
                PUBLIC fieldopt::reservoir
                ${Boost_LIBRARIES})

        add_compile_options(-std=c++11)

        add_executable(wicalc-adgprs dllmain.cpp)

        target_link_libraries(wicalc-adgprs
                fieldopt::wellindexcalculator-adgprs
                ${Boost_LIBRARIES})
    endif()
endif()


# Testing ==============================================================
if (BUILD_TESTING)

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
            )

    target_link_libraries(benchmark_wellindexcalculator
            fieldopt::wellindexcalculator
            Qt5::Core
            ${GTEST_BOTH_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT}
            ${GTEST_BOTH_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT}
            ${Boost_LIBRARIES})
endif()