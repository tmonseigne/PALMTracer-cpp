cmake_minimum_required(VERSION 3.15)
project(PALM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output des DLL/.so/.dylib dans le dossier Output/
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Output)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Output)

# Compilation position-independent pour les shared libraries (surtout utile sous Unix)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_SOURCE_DIR}/lib)

# === CPU ===
add_library(CPU SHARED
    CPU/main.cpp
    lib/Constants.hpp
    lib/ExportAPI.hpp
    lib/alignment/Alignment.hpp
    lib/alignment/Alignment.cpp
    lib/fit/Fitting.hpp
    lib/fit/Fitting.cpp
    lib/fit/GaussianFit.hpp
    lib/fit/GaussianFit.cpp
    lib/fit/SplineFit.hpp
    lib/fit/SplineFit.cpp
    lib/fit/Utils.hpp
    lib/fit/Utils.cpp
    lib/palm/ATrous.hpp
    lib/palm/ATrous.cpp
    lib/palm/SegmentationInfo.hpp
    lib/palm/Palm.hpp
    lib/palm/Palm.cpp
    lib/palm/Segmentation2D.hpp
    lib/palm/Segmentation2D.cpp
    lib/tracking/Node.hpp
    lib/tracking/Node.cpp
    lib/tracking/Parameters.hpp
    lib/tracking/TimeNodes.hpp
    lib/tracking/TimeNodes.cpp
    lib/tracking/Tracking.hpp
    lib/tracking/Tracking.cpp
)

# === Gestion manuelle OpenMP cross-plateforme ===
if(APPLE)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
        target_link_libraries(CPU PUBLIC OpenMP::OpenMP_CXX)
    else()
        message(STATUS "Configuring OpenMP manually for macOS (via Homebrew)")
        target_include_directories(CPU PUBLIC /opt/homebrew/include)
        target_compile_options(CPU PUBLIC -Xpreprocessor -fopenmp -I/opt/homebrew/include)
        target_link_libraries(CPU PUBLIC /opt/homebrew/lib/libomp.dylib)
    endif()
elseif(UNIX)
    find_package(OpenMP REQUIRED)
    target_link_libraries(CPU PUBLIC OpenMP::OpenMP_CXX)
elseif(WIN32)
    set_target_properties(CPU PROPERTIES
        LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/CPU/PALMTracer.def"
    )
endif()
