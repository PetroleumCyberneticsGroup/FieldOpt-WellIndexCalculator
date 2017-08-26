cmake_minimum_required(VERSION 2.8)
project(wellindexcalculator)

# Include: Eigen =======================================================
include_directories(${EIGEN3_INCLUDE_DIR})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(WIC_LIB_TARGET wellindexcalculator-adgprs-dbg)

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
    endif()

else() # LINUX --------------------------------------------------------

    if (NOT BUILD_WIC_ADGPRS)
        add_library(wellindexcalculator
                intersected_cell.cpp
                wellindexcalculator.cpp
                tests/wic_debug.hpp)

    elseif(BUILD_WIC_ADGPRS)
        add_library(${WIC_LIB_TARGET}
                SHARED
                dllmain.cpp
                intersected_cell.cpp
                wellindexcalculator.cpp
                tests/wic_debug.hpp)

    endif()

    if (NOT BUILD_WIC_ADGPRS)
        add_library(fieldopt::wellindexcalculator
                ALIAS ${PROJECT_NAME})

        target_link_libraries (wellindexcalculator
                PUBLIC fieldopt::reservoir
                ${Boost_LIBRARIES})

	    include_directories(${CMAKE_BINARY_DIR}/libraries/include/)

        add_compile_options(-std=c++11)

        add_executable(wicalc
                main.cpp
                tests/wic_debug.hpp)

        target_link_libraries(wicalc
                fieldopt::wellindexcalculator
                ${Boost_LIBRARIES})

        install(TARGETS wicalc wellindexcalculator
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib/static
                )

    elseif(BUILD_WIC_ADGPRS)
        add_library(fieldopt::${WIC_LIB_TARGET}
                ALIAS ${WIC_LIB_TARGET})
            # ALIAS ${PROJECT_NAME}-adgprs

        target_link_libraries(${WIC_LIB_TARGET}
                PUBLIC fieldopt::reservoir
                ${Boost_LIBRARIES})

        add_compile_options(-std=c++11)

        add_executable(wicalc-adgprs
                main.cpp
                dllmain.cpp
                tests/wic_debug.hpp)

        target_link_libraries(wicalc-adgprs
                PUBLIC fieldopt::${WIC_LIB_TARGET}
                ${Boost_LIBRARIES})
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