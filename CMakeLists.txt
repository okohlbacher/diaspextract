cmake_minimum_required(VERSION 3.16)
project(DIAspeXtract VERSION 1.4.1 LANGUAGES CXX)

# DIAspeXtract is a STANDALONE application that links against an existing OpenMS
# installation. It is deliberately not part of the OpenMS source tree: OpenMS is a
# prerequisite here, not a host.
#
# Point CMake at an OpenMS build or install with:
#   cmake -DOpenMS_DIR=/path/to/OpenMS/build ..

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# src/DnoiseMs1.h reproduces dnoise's MS1 window gate bit for bit only without floating-point
# contraction: an FMA flips 3 of dataset D's 24 boxes, and Apple clang on arm64 contracts by default.
# The header guards its own arithmetic; the flag keeps the translation units around it honest too, and
# matches the cluster's in-tree build, where OpenMS compiles with it but does not export it.
add_compile_options($<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>)

# The tool itself cannot be built without a PATCHED OpenMS: it includes
# <OpenMS/FORMAT/TdfMzCalibration.h>, which scripts/apply_openms_patches.sh installs into the
# OpenMS tree, and calls BrukerTimsFile::lastMzCalibration, which the same patch adds. No released
# OpenMS package contains either. The calibration model, however, is a self-contained header, so
# it can be compiled and tested anywhere -- which is what CI does on every platform. So is the
# dnoise MS1 port (src/DnoiseMs1.h).
option(DIASPEXTRACT_TESTS_ONLY "Build only the OpenMS-independent tests (portability CI)" OFF)
# The option was SPEXTRACTOR_TESTS_ONLY (SpeXtractor) and then DIASPEXTRACTOR_TESTS_ONLY (DIAspeXtractor). CMake would
# ignore an old name and take the full build, which then fails far from the cause in find_package(OpenMS), so either old
# name stops the configure here.
foreach(_old SPEXTRACTOR_TESTS_ONLY DIASPEXTRACTOR_TESTS_ONLY)
  if(DEFINED ${_old})
    message(FATAL_ERROR "${_old} is now DIASPEXTRACT_TESTS_ONLY (the tool was renamed DIAspeXtract): "
      "configure with -DDIASPEXTRACT_TESTS_ONLY=ON, and drop the old entry from an existing build directory with "
      "-U ${_old}.")
  endif()
endforeach()
# The dnoise port ships inside the tool, so both builds compile and register its unit test (it needs no OpenMS).
enable_testing()
add_executable(test_dnoise_ms1 tests/test_dnoise_ms1.cpp)
target_include_directories(test_dnoise_ms1 PRIVATE src)
add_test(NAME dnoise_ms1 COMMAND test_dnoise_ms1)
# src/TdfLoad.h reads analysis.tdf through sqlite3, so its test builds wherever SQLite is found (not on the Windows
# portability leg). It includes the calibration model as the patched OpenMS installs it; the shim maps that include onto
# src/TdfMzCalibration.h, the same file.
find_package(SQLite3 QUIET)
if(SQLite3_FOUND)
  configure_file(src/TdfMzCalibration.h ${CMAKE_CURRENT_BINARY_DIR}/shim/OpenMS/FORMAT/TdfMzCalibration.h COPYONLY)
  add_executable(test_tdf_load tests/test_tdf_load.cpp)
  target_include_directories(test_tdf_load PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/shim)
  target_link_libraries(test_tdf_load PRIVATE SQLite::SQLite3)
  add_test(NAME tdf_load COMMAND test_tdf_load)
endif()
if(DIASPEXTRACT_TESTS_ONLY)
  add_executable(test_calibration_cpp tests/test_calibration_cpp.cpp)
  target_include_directories(test_calibration_cpp PRIVATE src)
  add_test(NAME calibration_golden
           COMMAND test_calibration_cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/calibration_golden.json)
  return()
endif()

find_package(OpenMS REQUIRED)
find_package(OpenMP REQUIRED)   # the window loop and the picker call omp_* unconditionally
find_package(SQLite3 REQUIRED)   # analysis.tdf is read directly; a missing SQLite is a build error, not a fallback

add_executable(diaspextract src/diaspextract.cpp)
target_link_libraries(diaspextract PRIVATE OpenMS SQLite::SQLite3 OpenMP::OpenMP_CXX)
target_compile_definitions(diaspextract PRIVATE DIASPEXTRACT_VERSION="${PROJECT_VERSION}")
# The source carries the same number as a fallback for builds that do not go through this file
# (the in-tree OpenMS build); refuse to configure if the two have drifted.
file(STRINGS src/diaspextract.cpp _spx_fallback REGEX "^#define DIASPEXTRACT_VERSION \"")
if(NOT _spx_fallback MATCHES "\"${PROJECT_VERSION}\"")
  message(FATAL_ERROR "src/diaspextract.cpp's DIASPEXTRACT_VERSION fallback (${_spx_fallback}) != project VERSION ${PROJECT_VERSION}")
endif()

install(TARGETS diaspextract RUNTIME DESTINATION bin)

# Optional streaming .mzpeak input through the mzPeak C++ library (github.com/OpenMS/mzpeak).
# -DMZPEAK_ROOT=<checkout with include/ and a built libmzpeak in build*/>. Without it the tool
# still builds and reads .mzpeak only through OpenMS' non-streaming MzPeakFile.
if(MZPEAK_ROOT)
  find_library(MZPEAK_LIB NAMES mzpeak PATHS ${MZPEAK_ROOT} PATH_SUFFIXES build build-release buildfix REQUIRED)
  find_package(Arrow REQUIRED)
  find_package(Parquet REQUIRED)
  target_compile_definitions(diaspextract PRIVATE DIASPEXTRACT_WITH_MZPEAK)
  target_include_directories(diaspextract PRIVATE ${MZPEAK_ROOT}/include)
  target_link_libraries(diaspextract PRIVATE ${MZPEAK_LIB} Parquet::parquet_static Arrow::arrow_static zip)
endif()
