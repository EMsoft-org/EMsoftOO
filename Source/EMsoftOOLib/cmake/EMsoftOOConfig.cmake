include("${CMAKE_CURRENT_LIST_DIR}/EMsoftOOLibTargets.cmake")


if(@EMsoftOO_ENABLE_HDF5_SUPPORT@)
  include("${CMAKE_CURRENT_LIST_DIR}/EMsoftOOHDFLibTargets.cmake")
endif()

set(EMsoftOO_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")
set(EMsoftOO_LIB_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../lib;${CMAKE_CURRENT_LIST_DIR}/../../../bin")

set(EMsoftOO_Fortran_COMPILER_NAME @Fortran_COMPILER_NAME@)
set(EMsoftOO_BUILD_SHARED_LIBS "@BUILD_SHARED_LIBS@")

if (EMsoftOO_Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(EMsoftOO_Fortran_RT_Libs gfortran @EMsoftOO_FORTRAN_SUPPORT_LIBS@)
  set(EMsoftOO_Fortran_Lib_Dir @GFortran_LIB_DIR@)
endif()
