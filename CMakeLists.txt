cmake_minimum_required(VERSION 3.16)
project(DIAspeXtract VERSION 1.9.2 LANGUAGES CXX)

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
# contraction: an FMA flips 3 of the 24 gate boxes of a real run, and Apple clang on arm64 contracts by default.
# The header guards its own arithmetic; the flag keeps the translation units around it honest too, and
# matches the cluster's in-tree build, where OpenMS compiles with it but does not export it.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
  add_compile_options(/clang:-ffp-contract=off)   # clang-cl ignores the plain flag
else()
  add_compile_options($<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>)
endif()

# The tool itself cannot be built without a PATCHED OpenMS: it includes
# <OpenMS/FORMAT/TdfMzCalibration.h>, which scripts/apply_openms_patches.sh installs into the
# OpenMS tree, and calls BrukerTimsFile::lastMzCalibration, which the same patch adds. No released
# OpenMS package contains either. The calibration model, however, is a self-contained header, so
# it can be compiled and tested anywhere -- which is what CI does on every platform. So is the
# dnoise MS1 port (src/DnoiseMs1.h).
option(DIASPEXTRACT_TESTS_ONLY "Build only the OpenMS-independent tests (portability CI)" OFF)
# The dnoise port ships inside the tool, so both builds compile and register its unit test (it needs no OpenMS).
enable_testing()
add_executable(test_dnoise_ms1 tests/test_dnoise_ms1.cpp)
target_include_directories(test_dnoise_ms1 PRIVATE src)
add_test(NAME dnoise_ms1 COMMAND test_dnoise_ms1)
# The cgroup memory headroom (src/CgroupMem.h) is plain file parsing, tested on a fake /proc and /sys tree.
add_executable(test_cgroup_mem tests/test_cgroup_mem.cpp)
target_include_directories(test_cgroup_mem PRIVATE src)
add_test(NAME cgroup_mem COMMAND test_cgroup_mem)
# src/PortableMath.h: the exp/log the output path uses must give the same bits everywhere; the test pins one hash.
add_executable(test_portable_math tests/test_portable_math.cpp)
add_test(NAME portable_math COMMAND test_portable_math)
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
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # cl's OpenMP lacks taskloop and taskgroup even under /openmp:llvm, and 19.44 mis-shares a lambda's locals in tasks and
  # ends in internal compiler errors on this file. On Windows build the tool with clang-cl (MSVC ABI, full OpenMP).
  message(FATAL_ERROR "diaspextract needs clang-cl on Windows, not cl: configure with -T ClangCL or -DCMAKE_CXX_COMPILER=clang-cl")
endif()
find_package(SQLite3 REQUIRED)   # analysis.tdf is read directly; a missing SQLite is a build error, not a fallback

add_executable(diaspextract src/diaspextract.cpp)
# Windows: allocate through mimalloc instead of the Windows heap, whose locks cost 1.6x at 12 threads (Ryzen 9 7900X,
# HeLa diaPASEF: 8:07 -> 4:58 min, the spectra unchanged). mimalloc.dll is linked first and referenced
# (/INCLUDE:mi_version), so it is the exe's first import, and mimalloc-redirect.dll (which mimalloc.dll imports) sends
# every malloc of the process -- libOpenMS and the CRT included -- to it. Both DLLs go beside the exe.
option(DIASPEXTRACT_MIMALLOC "Windows: use mimalloc for every allocation of the process" ON)
if(WIN32 AND DIASPEXTRACT_MIMALLOC)
  include(FetchContent)
  set(MI_BUILD_SHARED ON CACHE BOOL "" FORCE)
  set(MI_BUILD_STATIC OFF CACHE BOOL "" FORCE)
  set(MI_BUILD_OBJECT OFF CACHE BOOL "" FORCE)
  set(MI_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(mimalloc
    URL https://github.com/microsoft/mimalloc/archive/refs/tags/v3.6.1.tar.gz
    URL_HASH SHA256=2672e6a99d4579fb4bb32d9de0d63aa424ce600d31fbfaf06eb74d940d2e07a6)
  FetchContent_MakeAvailable(mimalloc)
  target_link_libraries(diaspextract PRIVATE mimalloc)   # before every other library
  target_link_options(diaspextract PRIVATE /INCLUDE:mi_version)
  add_custom_command(TARGET diaspextract POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:mimalloc> ${mimalloc_SOURCE_DIR}/bin/mimalloc-redirect.dll
            $<TARGET_FILE_DIR:diaspextract>)
  install(FILES $<TARGET_FILE:mimalloc> ${mimalloc_SOURCE_DIR}/bin/mimalloc-redirect.dll DESTINATION bin)
  install(FILES ${mimalloc_SOURCE_DIR}/LICENSE DESTINATION LICENSES/third-party RENAME mimalloc-LICENSE)
endif()
target_link_libraries(diaspextract PRIVATE OpenMS OpenMS_CLI SQLite::SQLite3 OpenMP::OpenMP_CXX)   # OpenMS_CLI: TOPPBase
target_compile_definitions(diaspextract PRIVATE DIASPEXTRACT_VERSION="${PROJECT_VERSION}")
# Windows: a version resource with the product name and version. SignPath signs the release .exe only when both are
# present and match its artifact configuration (.signpath/artifact-configurations/initial.xml).
if(WIN32)
  enable_language(RC)
  configure_file(src/diaspextract.rc.in ${CMAKE_CURRENT_BINARY_DIR}/diaspextract.rc @ONLY)
  target_sources(diaspextract PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/diaspextract.rc)
endif()
# The source carries the same number as a fallback for builds that do not go through this file
# (the in-tree OpenMS build); refuse to configure if the two have drifted.
file(STRINGS src/diaspextract.cpp _spx_fallback REGEX "^#define DIASPEXTRACT_VERSION \"")
if(NOT _spx_fallback MATCHES "\"${PROJECT_VERSION}\"")
  message(FATAL_ERROR "src/diaspextract.cpp's DIASPEXTRACT_VERSION fallback (${_spx_fallback}) != project VERSION ${PROJECT_VERSION}")
endif()

install(TARGETS diaspextract RUNTIME DESTINATION bin)

# Optional streaming .mzpeak input through the mzPeak C++ library (github.com/OpenMS/mzpeak).
# -DMZPEAK_ROOT=<checkout with include/ and a built libmzpeak in build/>. Without it the tool
# still builds and reads .mzpeak only through OpenMS' non-streaming MzPeakFile.
if(MZPEAK_ROOT)
  find_library(MZPEAK_LIB NAMES mzpeak PATHS ${MZPEAK_ROOT} PATH_SUFFIXES build REQUIRED)
  find_package(Arrow REQUIRED)
  find_package(Parquet REQUIRED)
  target_compile_definitions(diaspextract PRIVATE DIASPEXTRACT_WITH_MZPEAK)
  target_include_directories(diaspextract PRIVATE ${MZPEAK_ROOT}/include)
  target_link_libraries(diaspextract PRIVATE ${MZPEAK_LIB} Parquet::parquet_static Arrow::arrow_static zip)
endif()
