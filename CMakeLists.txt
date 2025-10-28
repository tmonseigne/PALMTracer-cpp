cmake_minimum_required(VERSION 3.15)
project(PALMTracer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output des DLL/.so/.dylib dans le dossier Output/
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Output)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Output)

# Compilation position-independent pour les shared libraries (surtout utile sous Unix)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_SOURCE_DIR}/lib)

# === Fonction de configuration d'OpenMP pour Unix & Mac ===
function(configure_openmp_for_target TARGET_NAME)
    if(APPLE)
        find_package(OpenMP)
        if(OpenMP_CXX_FOUND)
            message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
            target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
        else()
            message(STATUS "Configuring OpenMP manually for macOS (via Homebrew)")
            target_include_directories(${TARGET_NAME} PUBLIC /opt/homebrew/include)
            target_compile_options(${TARGET_NAME} PUBLIC -Xpreprocessor -fopenmp -I/opt/homebrew/include)
            target_link_libraries(${TARGET_NAME} PUBLIC /opt/homebrew/lib/libomp.dylib)
        endif()
    elseif(UNIX)
        find_package(OpenMP REQUIRED)
        target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endfunction()

# === Files ===
file(GLOB_RECURSE CPU_HDRS CONFIGURE_DEPENDS "lib/*.hpp")
list(FILTER CPU_HDRS EXCLUDE REGEX "[/\\]cuda[/\\]")
file(GLOB_RECURSE CPU_SRC CONFIGURE_DEPENDS "lib/*.cpp")
list(FILTER CPU_SRC EXCLUDE REGEX "[/\\]cuda[/\\]")
file(GLOB_RECURSE TESTS_FILES CONFIGURE_DEPENDS "Tests_CPU/*.cpp" "Tests_CPU/*.hpp")

# === CPU ===
add_library(CPU SHARED CPU/main.cpp ${CPU_SRC} ${CPU_HDRS})

# === Gestion manuelle OpenMP cross-plateforme ===
configure_openmp_for_target(CPU)

if(WIN32)
    set_target_properties(CPU PROPERTIES LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/CPU/PALMTracer.def")
endif()

find_package(FFTW3 REQUIRED)
target_link_libraries(CPU PUBLIC FFTW3::fftw3)

# === Google Test (test_cpu) ===
enable_testing()

add_executable(test_cpu ${TESTS_FILES} ${CPU_SRC} ${CPU_HDRS})

# === Gestion manuelle OpenMP cross-plateforme ===
configure_openmp_for_target(test_cpu)


find_package(GTest REQUIRED)
find_package(TIFF REQUIRED)
target_link_libraries(test_cpu
    PRIVATE
        FFTW3::fftw3
        GTest::gtest
        GTest::gtest_main
        TIFF::TIFF
)

include(GoogleTest)

# Args GTest (facultatif)
set(GTEST_ARGS
  --gtest_color=yes
  --gtest_print_time=1
)

# IMPORTANT : NE PAS faire add_test(NAME CPU_Tests COMMAND CPU_Tests ...)
# On laisse CMake découvrir chaque test gtest → tests CTest individuels.
gtest_discover_tests(test_cpu
  EXTRA_ARGS        ${GTEST_ARGS}
  WORKING_DIRECTORY $<TARGET_FILE_DIR:test_cpu>
  DISCOVERY_MODE    PRE_TEST           # robuste en CI
  PROPERTIES        TIMEOUT 120
)

# Copy des éléments dans le dossier de build (en cas de pépin avec le CI, on en créé plus)
file(COPY Tests_CPU/Input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Tests_CPU)
# Créer des dossiers Output au moment de la configuration (en cas de pépin avec le CI, on en créé plus)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Output)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Tests_CPU/Output)
