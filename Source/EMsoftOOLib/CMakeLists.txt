project(EMsoftOOLib)

include (${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake
${EMsoftOOLib_BINARY_DIR}/StringConstants.cmake)

# configure_file(${EMsoftLib_SOURCE_DIR}/StringConstants.cmake
#                ${EMsoftLib_BINARY_DIR}/StringConstants.cmake @ONLY)

cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOOLib_BINARY_DIR}/StringConstants.cmake
CONFIGURED_TEMPLATE_PATH   ${EMsoftOOLib_SOURCE_DIR}/StringConstants.cmake )


add_custom_target(EMsoftOOLib_Generate_StringConstants ALL
                  "${CMAKE_COMMAND}" -P "${EMsoftOOLib_BINARY_DIR}/StringConstants.cmake"
                  COMMENT "Generating StringConstants files...."
                  )


#------------------------------------------------------------------------------
# Setup all the source files
set(EMsoftOOLib_SRCS

  ${EMsoftOOLib_SOURCE_DIR}/mod_kinds.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_global.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_memory.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_io.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_timing.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_quaternions.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_rng.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_Lambert.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_symmetry.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_symmetry2D.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_crystallography.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_HDFsupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_HDFnames.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_postscript.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_misc.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_image.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_image_io.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_bmp.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_tif.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_lzw.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_pgm.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_math.f90
  ${EMsoftOOLib_SOURCE_DIR}/spsort.f
  ${EMsoftOOLib_SOURCE_DIR}/mod_rotations.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_povray.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_so3.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_dirstats.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_notifications.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_others.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_TB.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_diffraction.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_HOLZ.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_Kikuchi.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_gvectors.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_kvectors.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_MRC.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_initializers.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_filters.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_fftw3.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_fft_wrap.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_imageOPs.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_DSHT.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_patterns.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_vendors.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_DIsupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_OMPsupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_bobyqa.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_bobyqa_refinement.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_FitOrientations.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_JSONsupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_defect.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_MuellerCalculus.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_LaueSupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_NLPAR.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_SEMwrappers.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_wrappers.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_xrd.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_DREAM3D.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_PGA3D.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_PGA3Dsupport.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_polyhedra.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_MCA.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_STL.f90
  ${EMsoftOOLib_SOURCE_DIR}/mod_demag.f90
  # ${EMsoftOOLib_SOURCE_DIR}/mod_shapes.f90

  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_sampleRFZ.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ConvertOrientations.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_Disorientations.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_MCfiles.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_MPfiles.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_EBSDmaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_EBSD.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ECPmaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ECP.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_TKDmaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_MPoverlap.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_reflectors.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_Kosselmaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_kinematical.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_OrientationViz.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_DIfiles.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ADP.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ANG.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_CTF.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_OSM.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_dpmerge.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_DIsetting.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_getEulers.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_DIpreview.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_getMP.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_HDFFileInfo.f90
#  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_FitOrientation.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ECCI.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_CPLMmaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_Lauemaster.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_Laue.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_LaueDCT.f90
  ${EMsoftOOLib_SOURCE_DIR}/program_mods/mod_ShapeAmplitude.f90

  ${EMsoftOOLib_Additional_SRCS}
)

# set(EMsoftOOLib_C_SRCS
#   ${EMsoftOOLib_SOURCE_DIR}/mbir.c
#   ${EMsoftOOLib_SOURCE_DIR}/mbirHeader.h
#   ${EMsoftOOLib_SOURCE_DIR}/denoise.c
#   ${EMsoftOOLib_SOURCE_DIR}/bclsf90wrapper.c
#   ${EMsoftOOLib_SOURCE_DIR}/cs_malloc.c
#   ${EMsoftOOLib_SOURCE_DIR}/cs_util.c
#   ${EMsoftOOLib_SOURCE_DIR}/cs.h
#   ${EMsoftOOLib_Additional_C_SRCS}
# )

# set(EMsoftOOLib_Cpp_SRCS
#   ${EMsoftOOLib_SOURCE_DIR}/sht_file.cpp
#   ${SHTfile_BINARY_DIR}/sht_file.hpp
# )

# if(NOT WIN32 AND BUILD_SHARED_LIBS)
#   set_source_files_properties(${EMsoftOOLib_C_SRCS} PROPERTIES COMPILE_FLAGS -fPIC)
#   set_source_files_properties(${EMsoftOOLib_Cpp_SRCS} PROPERTIES COMPILE_FLAGS -fPIC)
# endif()

if(EMsoftOO_ENABLE_TESTING)
  set(EMsoftOO_TESTING_DIR "${EMsoftOO_BINARY_DIR}/Testing")
  file(MAKE_DIRECTORY ${EMsoftOO_TESTING_DIR})
endif()

message(STATUS "HDF5_VERSION_MAJOR: ${HDF5_VERSION_MAJOR}")

if(${HDF5_VERSION_MAJOR} EQUAL 1.10)
set(EMsoft_HDF5_110_DEFINES "
! For HDF5 1.10, these are needed
! Define the \"kind\" parameters for Fortran_REAL_4
  integer,parameter                     :: Fortran_REAL_4 = 4

! Define the \"kind\" parameters for Fortran_REAL_4
  integer,parameter                     :: Fortran_REAL_8 = 8
")
endif()

cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOOLib_BINARY_DIR}/mod_EMsoft.f90
                              CONFIGURED_TEMPLATE_PATH   ${EMsoftOOLib_SOURCE_DIR}/mod_EMsoft.in.f90 )
# -------------------------------------------------------------------------------
# Add files GENERATED DURING CMAKE TIME
set(EMsoftOOLib_SRCS
  ${EMsoftOOLib_SRCS}
  ${EMsoftOOLib_BINARY_DIR}/mod_EMsoft.f90
  ${EMsoftOOLib_BINARY_DIR}/stringconstants.f90
)

#-------------------------------------------------------------------------------
# Put just the *.c files into their own library so Visual Studio will compile
# everything correctly. VS can NOT have both fortran and C files in the same target
# add_library(EMsoftOOLib_C STATIC ${EMsoftOOLib_C_SRCS})
# set_target_properties (EMsoftOOLib_C PROPERTIES
#   LINKER_LANGUAGE Fortran
#    # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
#   FOLDER EMsoftOOPublic
# )
# add_library(EMsoftOOLib_Cpp STATIC ${EMsoftOOLib_Cpp_SRCS})
# set_target_properties (EMsoftOOLib_Cpp PROPERTIES
#   LINKER_LANGUAGE Fortran
#    # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
#   FOLDER EMsoftOOPublic
# )
# target_include_directories(EMsoftOOLib_Cpp PRIVATE ${SHTfile_BINARY_DIR})

GetHDF5LinkLibraries(EMSOFTOO)

add_library(EMsoftOOLib ${LIB_TYPE} ${EMsoftOOLib_SRCS})

target_link_libraries(EMsoftOOLib
                # EMsoftOOLib_C
                # EMsoftOOLib_Cpp
                ${EMSOFTOO_hdf5LinkLibs}
                jsonfortran
                bcls::bcls
)

if (Fortran_COMPILER_NAME MATCHES "ifort.*")
  target_link_libraries(EMsoftOOLib
                  ${MKL_LIBRARIES}
  )
endif()

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  target_link_libraries(EMsoftOOLib
                  ${FFTW3_LIBRARY}
                  ${EMsoftOO_BLAS_LAPACK_LIBS}
                  ${EMsoftOO_FORTRAN_SUPPORT_LIBS}
                  gomp
  )
endif()

target_include_directories(EMsoftOOLib
                          PUBLIC
                            ${CLFortran_INCLUDE_DIRS}
                            ${FFTW3_INCLUDE_DIR}
                            ${HDF5_INCLUDE_DIR}
                            ${HDF5_INCLUDE_DIR_FORTRAN}
                          PRIVATE
                            "${EMsoftOO_SOURCE_DIR}/Source"
                          )

set_target_properties (EMsoftOOLib PROPERTIES
                  LINKER_LANGUAGE Fortran
                   # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
                  FOLDER EMsoftOOPublic
                  )

set_source_files_properties(${EMsoftOOLib_BINARY_DIR}/stringconstants.f90 PROPERTIES GENERATED TRUE)
add_dependencies(EMsoftOOLib EMsoftOOLib_Generate_StringConstants)


# On Windows we need a few extra compile and linker flags. By making these PUBLIC they should
# flow down to every dependent target
if(WIN32)
  target_compile_options(EMsoftOOLib PUBLIC
                          $<$<COMPILE_LANGUAGE:Fortran>:/assume:byterecl>
                          $<$<COMPILE_LANGUAGE:Fortran>:/assume:nounderscore>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qopenmp>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:11082>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:7025>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qip>
                          $<$<COMPILE_LANGUAGE:Fortran>:/traceback>
                          #$<$<COMPILE_LANGUAGE:Fortran>:/warn:all>
                        #  $<$<COMPILE_LANGUAGE:Fortran>:/heap-arrays>
                          )
  target_link_options(EMsoftOOLib PUBLIC /STACK:2000000000 /machine:x64)
else()
  if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # First, Flags that should be used under any build condition
    target_compile_options(EMsoftOOLib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp>
        $<$<COMPILE_LANGUAGE:Fortran>:-fno-sign-zero>
        $<$<COMPILE_LANGUAGE:Fortran>:-fno-f2c>
        $<$<COMPILE_LANGUAGE:Fortran>:-std=gnu>
        $<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>
    )
    # Next, DEBUG compiler options.
    target_compile_options(EMsoftOOLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=bounds>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
    )
    # Next, RELEASE compiler options.
    target_compile_options(EMsoftOOLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-funroll-all-loops>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    )

  elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # First, Flags that should be used under any build condition
    target_compile_options(EMsoftOOLib PUBLIC
                                            #$<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-qopenmp>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-assume byterecl>
                                            $<$<COMPILE_LANGUAGE:Fortran>:/assume:nounderscore>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-stand f08 -diag-disable=7025>
                                            $<$<COMPILE_LANGUAGE:Fortran>:/traceback>
    )
    # Next, DEBUG compiler options.
    target_compile_options(EMsoftOOLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-debug all>
                                            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
                                            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
    )
    # Next, RELEASE compiler options.
    target_compile_options(EMsoftOOLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>)
  else()
    message(WARNING "Fortran compiler not recognized. ${Fortran_COMPILER_NAME}")
  endif()
endif()


#-------------------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
EMsoftOO_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
  INSTALL(TARGETS EMsoftOOLib # EMsoftOOLib_C EMsoftOOLib_Cpp
    COMPONENT Applications
    EXPORT EMsoftOOLibTargets
    RUNTIME DESTINATION ${install_dir}
    LIBRARY DESTINATION ${lib_install_dir}
    ARCHIVE DESTINATION ${lib_install_dir}
    BUNDLE DESTINATION "."
  )
#endif()

install(FILES ${EMsoftOOLib_BINARY_DIR}/EMsoftOOStringConstants.h
      DESTINATION ${include_install_dir}/EMsoftOOLib
      COMPONENT Applications
)

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libEMsoftOOLib.dylib")
  get_property(EMsoftOOBinaryDir GLOBAL PROPERTY EMsoftOOBinaryDir)
  set(OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT "${EMsoftOO_BINARY_DIR}/EMsoftOOLib/EMsoftOOLib_CompleteTool.cmake")
  set(OPTIMIZE_BUNDLE_SHELL_SCRIPT            "${EMsoftOO_BINARY_DIR}/EMsoftOOLib/EMsoftOOLib_OptimizeTool.sh")
  set(PROJECT_INSTALL_DIR "lib")
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.cmake.in"
                "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" @ONLY IMMEDIATE)
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.sh.in"
                "${OPTIMIZE_BUNDLE_SHELL_SCRIPT}" @ONLY IMMEDIATE)
  #install(SCRIPT "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" COMPONENT ${QAB_COMPONENT})
endif()


# --------------------------------------------------------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOLibTargetsConfigVersion.cmake"
  VERSION ${EMsoftOO_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT EMsoftOOLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOLibTargets.cmake"
    NAMESPACE EMsoftOO::
  )
#endif()

configure_file(cmake/EMsoftOOConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOConfig.cmake"
  @ONLY
)

set(ConfigPackageLocation ${lib_install_dir}/cmake/EMsoftOO)

#if(BUILD_SHARED_LIBS)
  install(EXPORT EMsoftOOLibTargets
    FILE
      EMsoftOOLibTargets.cmake
    NAMESPACE
      EMsoftOO::
    DESTINATION
      ${ConfigPackageLocation}
  )

#endif()

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOLibTargetsConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
