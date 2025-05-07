cmake_minimum_required(VERSION 3.13)

project(kklib VERSION 1.0.0)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(kklib STATIC)

set(public_headers
    include/kklib/constants.hpp
    include/kklib/graph_data_structure.hpp
    include/kklib/graph.hpp
    include/kklib/metapath.hpp
    include/kklib/metascheme.hpp
    include/kklib/mpi_helper.hpp
    include/kklib/node2vec.hpp
    include/kklib/path.hpp
    include/kklib/ppr.hpp
    include/kklib/static_comp.hpp
    include/kklib/storage.hpp
    include/kklib/type.hpp
    include/kklib/util.hpp
    include/kklib/walk.hpp
)

target_sources(kklib
  PRIVATE
  src/metascheme.cpp
  src/mpi_helper.cpp
  src/util.cpp
)

set_property(TARGET kklib PROPERTY POSITION_INDEPENDENT_CODE ON)

# MPI    
find_package(MPI REQUIRED)
target_link_libraries(kklib PRIVATE MPI::MPI_CXX) 

# OpenMP
find_package(OpenMP)
target_link_libraries(kklib PRIVATE OpenMP::OpenMP_CXX)

# TBB
find_package(TBB REQUIRED)
target_link_libraries(kklib PRIVATE TBB::tbb)

target_include_directories(kklib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS kklib
    EXPORT "kklibTargets"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/kklib)

install(EXPORT "kklibTargets"
  FILE "kklibTargets.cmake"
  NAMESPACE kklib::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kklib
)

add_library(kklib::kklib STATIC IMPORTED)

set_target_properties(kklib::kklib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
)

set_target_properties(kklib PROPERTIES PUBLIC_HEADER "${public_headers}")
set_target_properties(kklib PROPERTIES DEBUG_POSTFIX "d")

# ===================================
# Make a Configuration Package
# ===================================
include(CMakePackageConfigHelpers)

export(EXPORT kklibTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/kklibTargets.cmake"
       NAMESPACE kklib::
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/kklibConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kklib
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/kklibConfigVersion.cmake"
  VERSION "${version}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/kklibConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/kklibConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kklib
)

set_property(TARGET kklib PROPERTY VERSION ${version})

set_target_properties(kklib PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})

# ===================================
# Test Targets
# ===================================
option(WITH_TESTS "Build integration test programs." OFF)
if(WITH_TESTS)
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/lib/gtest")
      message(FATAL_ERROR
      "Missing GTest library in lib/gtest "
      "Please run `git submodule update --init` to fetch the submodule.")
    endif()
    
    add_subdirectory(${PROJECT_SOURCE_DIR}/lib/gtest)
    
    enable_testing()

    find_package(OpenMP REQUIRED)
    find_package(MPI REQUIRED)
    find_package(gdsb REQUIRED)
    
    # TODO: not all tests actually require MPI. Those tests could be decoupled
    # from the MPI requirement and linkage.
    set(TEST_TARGETS test_bound test_deepwalk test_graph test_metapath test_node2vec test_outlier test_path test_ppr test_storage test_util test_walker)
    foreach(test_target ${TEST_TARGETS})
      add_executable(${test_target}
        tests/${test_target}.cpp
      )

      
      # TODO: these compile definitions should be removed. Production code
      # should not contain any pre-processor directives; instead they should
      # only active when testing.
      target_compile_definitions(${test_target} PRIVATE UNIT_TEST)
      target_compile_definitions(${test_target} PRIVATE COLLECT_WALK_SEQUENCE)
      target_compile_definitions(${test_target} PRIVATE COLLECT_WALKER_INIT_STATE)

      target_link_libraries(${test_target} PRIVATE gtest)

      # OpenMP
      target_link_libraries(${test_target} PRIVATE OpenMP::OpenMP_CXX)

      # MPI
      target_link_libraries(${test_target} PRIVATE MPI::MPI_CXX)
    
      target_link_libraries(${test_target} PRIVATE kklib)

      target_link_libraries(${test_target} PRIVATE gdsb::gdsb)
        
      add_test(
        NAME ${test_target}
        COMMAND ${test_target}
      )
    endforeach()
endif()

# ===================================
# Build Applications & Tools
# ===================================
option(KKLIB_APPLICATIONS "Build example applications." OFF)
if(KKLIB_APPLICATIONS)

  find_package(OpenMP)
  find_package(MPI REQUIRED)
  find_package(gdsb REQUIRED)

  # Setup and build all MPI applications
  set(MPI_APPS biased_walk deepwalk dynamic_walk metapath node2vec ppr simple_walk)
  FOREACH(mpi_app ${MPI_APPS})
    add_executable(${mpi_app}
      applications/${mpi_app}.cpp
    )
    target_link_libraries(${mpi_app} PRIVATE kklib)
    target_link_libraries(${mpi_app} PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(${mpi_app} PRIVATE MPI::MPI_CXX)
    target_link_libraries(${mpi_app} PRIVATE gdsb::gdsb)
  
    install(
      TARGETS ${mpi_app}
      RUNTIME DESTINATION bin
    )
  ENDFOREACH()
endif()
