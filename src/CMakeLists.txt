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




if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset  "${Esc}[m")
  set(ColorBold   "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()



message ( "\n${ColorBold}${White}***********************************************************************" )
message ( "*        ${Blue}_______   ___       ___   ___   ___   _______     _______    ${White}*" )
message ( "*       ${Blue}////////  ////      ////  ////  ////  ////////\\   ////////    ${White}*" )
message ( "*      ${Blue}////__    ////      ////  ////  ////  ////  ////  ////__       ${White}*" )
message ( "*     ${Blue}///////   ////      ////  ////  ////  ////  ////   \\\\\\\\\\\\\\      ${White}*" )
message ( "*    ${Blue}////      ////___   ////__////  ////  ////__////  _____////      ${White}*" )
message ( "*   ${Blue}////      ////////   \\////////  ////  /////////   ////////        ${White}*" )
message ( "*                                                                     *" )
message ( "***********************************************************************${ColorReset}\n" )

        

        

### TESTS ###
add_executable ( DistributionFunction_Test
  test_DistributionFunction.cxx DistributionFunction.cxx
)

add_executable ( Field_Test
  test_Field.cxx Field.cxx
)

add_executable ( LBMKernel_Test
  test_LBMKernel.cxx LBMKernel.cxx CollisionOperator.cxx StreamingOperator.cxx BoundaryCondition.cxx DistributionFunction.cxx Field.cxx
)


### LIBRARY ###
add_library ( fluids SHARED SimManager.cxx IncompressibleKernel.cxx PoissonSolver.cxx ThermalKernel.cxx DataKeeper.cxx DataProcessor.cxx
                            LBMKernel.cxx CollisionOperator.cxx StreamingOperator.cxx BoundaryCondition.cxx DistributionFunction.cxx Field.cxx )
install (
  TARGETS fluids
  LIBRARY DESTINATION lib
)



### DOC ###
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



### PYTHON MODULE ###

# Détection des prérequis
find_package ( PythonInterp 3.0 REQUIRED )
find_package ( pybind11 REQUIRED )


# Construction du module Python (import _vortex) avec Pybind11
# NOTA: décommenter au besoin si installation de Pybind11 hors système
#include_directories ( /chemin/vers/pybind11/pybind11.h/si/installation/hors/path/std )
pybind11_add_module ( _fluids MODULE Interface.cxx InterfaceLBM.cxx InterfaceCommon.cxx )
target_link_libraries ( _fluids PUBLIC fluids )
install (
  TARGETS _fluids
  LIBRARY DESTINATION lib/python${PYTHON_VERSION}
)


# Local Variables:
# mode: cmake
# coding: utf-8
# End:
