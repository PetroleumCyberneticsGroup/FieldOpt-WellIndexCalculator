// dllmain.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#if _WIN32
// Including SDKDDKVer.h defines the highest available Windows platform.
// If you wish to build your application for a previous Windows platform, include WinSDKVer.h and
// set the _WIN32_WINNT macro to the platform you wish to support before including SDKDDKVer.h.

#include <SDKDDKVer.h>

#define WIN32_LEAN_AND_MEAN   // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>
#endif

//#ifndef DWORD
//#define WINAPI
//typedef unsigned long DWORD;
//typedef short WCHAR;
//typedef void * HANDLE;
//#define MAX_PATH    PATH_MAX
//typedef unsigned char BYTE;
//typedef unsigned short WORD;
//typedef unsigned int BOOL;
//#endif

// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the WELLINDEXCALCULATOR_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// WELLINDEXCALCULATOR_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef WELLINDEXCALCULATOR_EXPORTS
#if _WIN32
        #ifdef __cplusplus
        #define WELLINDEXCALCULATOR_API  extern "C" __declspec(dllexport)
        #else
        #define WELLINDEXCALCULATOR_API  __declspec(dllexport)
        #endif
    #else
        #ifdef __cplusplus
        #define WELLINDEXCALCULATOR_API  extern "C"
        #else
        #define WELLINDEXCALCULATOR_API
        #endif
    #endif
#else
#if _WIN32
#ifdef __cplusplus
    #define WELLINDEXCALCULATOR_API extern "C" __declspec(dllimport)
    #else
    #define WELLINDEXCALCULATOR_API __declspec(dllimport)
    #endif
#else
#ifdef __cplusplus
#define WELLINDEXCALCULATOR_API extern "C"
#else
#define WELLINDEXCALCULATOR_API
#endif
#endif
#endif

// This is the constructor of a class that has been exported.
// This class is exported from the WellIndexCalculator.dll
//class WELLINDEXCALCULATOR_API CWellIndexCalculator {
//public:
//  CWellIndexCalculator(void);
// TODO: add your methods here.
//};

// This is an example of an exported variable
void* grid;

// This is an example of an exported function.
WELLINDEXCALCULATOR_API int computeWellIndices(
    const char* basepth,
    const double* heel, const double* toe, const double* wellbore_radius,
    int* n, int* i, int* j, int* k, double* wi);

WELLINDEXCALCULATOR_API int getBlockCenters(
    const char* basepth,
    const int* heel,  const int* toe, double* heelxyz, double* toexyz);

WELLINDEXCALCULATOR_API int getBoundaryVertices(
    const char* filepth,
    int* npnts, double* xes, double* yes, double* zes);