cmake_minimum_required(VERSION 3.2)
project(wellindexcalculator LANGUAGES CXX)

option(BUILD_WIC_ADGPRS "Build WIC as ADGPRS plugin (shared lib, WIN/LINUX)" ON)

# Include: Eigen =======================================================
include_directories(${EIGEN3_INCLUDE_DIR})

# wellindexcalculator lib ==============================================
if (_WIN32) # WINDOWS
    find_library(SHLWAPI_LIBRARY NAMES Shlwapi)

    add_library(wellindexcalculator SHARED dllmain.cpp)

    target_link_libraries(wellindexcalculator
            fieldopt::wellindexcalculator
            ${Boost_LIBRARIES}
            ${SHLWAPI_LIBRARY})

else() # LINUX --------------------------------------------------------
    if (NOT BUILD_WIC_ADGPRS)
        add_library(wellindexcalculator
                intersected_cell.cpp
                wellindexcalculator.cpp)

    elseif(BUILD_WIC_ADGPRS)
        add_library(wellindexcalculator SHARED
                dllmain.cpp
                intersected_cell.cpp
                wellindexcalculator.cpp)
    endif()

    add_library(fieldopt::wellindexcalculator ALIAS ${PROJECT_NAME})

    target_link_libraries (wellindexcalculator
            PUBLIC fieldopt::reservoir
            ${Boost_LIBRARIES})

    target_compile_features(wellindexcalculator
            PUBLIC cxx_auto_type
            PUBLIC cxx_range_for)
endif()

# Standalone WellIndexCalculator executable (WIN/Linux) ================
add_executable(wicalc main.cpp)

if (_WIN32)
    find_library(SHLWAPI_LIBRARY NAMES Shlwapi)

    target_link_libraries(wicalc
            fieldopt::wellindexcalculator
            ${Boost_LIBRARIES}
            ${SHLWAPI_LIBRARY})
else()
    target_link_libraries(wicalc
            fieldopt::wellindexcalculator
            ${Boost_LIBRARIES})
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
