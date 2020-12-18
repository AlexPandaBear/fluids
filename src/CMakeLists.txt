cmake_minimum_required ( VERSION 3.0 )
cmake_policy ( VERSION 3.0 )

project ( fluids CXX )
enable_testing ()


# RPATH settings
set ( CMAKE_SKIP_BUILD_RPATH  FALSE )
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
list ( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir )
if ( "${isSystemDir}" STREQUAL "-1" )
  set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif ( "${isSystemDir}" STREQUAL "-1" )


add_library ( vf SHARED SimManager.cxx IncompressibleKernel.cxx PoissonSolver.cxx ThermalKernel.cxx DataKeeper.cxx DataProcessor.cxx )
install (
  TARGETS vf
  LIBRARY DESTINATION lib
)


add_executable ( DistributionFunction_Test
  test_DistributionFunction.cxx DistributionFunction.cxx
)

add_executable ( Field_Test
  test_Field.cxx Field.cxx
)

add_executable ( LBMKernel_Test
  test_LBMKernel.cxx LBMKernel.cxx CollisionOperator.cxx StreamingOperator.cxx BoundaryCondition.cxx DistributionFunction.cxx Field.cxx
)


find_program ( DOXYGEN_EXECUTABLE
  NAMES doxygen
  DOC /usr/ 
  HINTS $ENV{DOXYGEN_ROOT}
  PATH_SUFFIXES bin
)

FIND_PROGRAM (DOXYGEN_DOT_EXECUTABLE
  NAMES dot
  DOC /usr/bin
  PATHS $ENV{DOT_HOME}
)

add_custom_target ( reference_doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/Doxyfile
)


# Détection des prérequis
find_package ( PythonInterp 3.0 REQUIRED )
find_package ( pybind11 REQUIRED )


# Construction du module Python (import _vortex) avec Pybind11
# NOTA: décommenter au besoin si installation de Pybind11 hors système
#include_directories ( /chemin/vers/pybind11/pybind11.h/si/installation/hors/path/std )
pybind11_add_module ( _vf MODULE Interface.cxx )
target_link_libraries ( _vf PUBLIC vf )
install (
  TARGETS _vf
  LIBRARY DESTINATION lib/python${PYTHON_VERSION}
)


# Local Variables:
# mode: cmake
# coding: utf-8
# End:
