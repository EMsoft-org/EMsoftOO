include(${PROJECT_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake)

project(EMsoftOOHDFLib)



set(install_dir "bin")
set(lib_install_dir "lib")

if(APPLE)
  get_property(EMsoftOO_PACKAGE_DEST_PREFIX GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX)
  set(install_dir "${EMsoftOO_PACKAGE_DEST_PREFIX}/bin")
  set(lib_install_dir "${EMsoftOO_PACKAGE_DEST_PREFIX}/lib")
elseif(WIN32)
  set(install_dir ".")
  set(lib_install_dir ".")
endif()



set(EMsoftOOHDFLib_SRCS
  # ${EMsoftOOHDFLib_SOURCE_DIR}/EBSDmod.f90
	# ${EMsoftOOHDFLib_SOURCE_DIR}/initializersHDF.f90
 # 	${EMsoftOOHDFLib_SOURCE_DIR}/EBSDDImod.f90
 	${EMsoftOOHDFLib_SOURCE_DIR}/mod_HDFsupport.f90
 # 	${EMsoftOOHDFLib_SOURCE_DIR}/NameListHDFwriters.f90
 #  ${EMsoftOOHDFLib_Additional_SRCS}
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/commonmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/EBSDiomod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/ECPiomod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/ECPmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/ECPmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/EMh5ebsd.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/FitOrientations.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/InitializersQCHDF.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/Lauemod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/patternmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/PFInversionHDF.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/QCmodHDF.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/TKDDImod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/TKDmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/utilities.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/hhHDFmod.f90
 #  ${EMsoftOOHDFLib_SOURCE_DIR}/MDsubroutines.f90

# 	${EMsoftOOHDFLib_SOURCE_DIR}/EMdymodHDF.f90
)

GetHDF5LinkLibraries(EMSOFTOO)

add_library(EMsoftOOHDFLib ${LIB_TYPE} ${EMsoftOOHDFLib_SRCS})

target_include_directories(EMsoftOOHDFLib 
                          PUBLIC 
                            ${HDF5_INCLUDE_DIR}
                            ${HDF5_INCLUDE_DIR_FORTRAN}
                            $<BUILD_INTERFACE:${EMsoftOOLib_BINARY_DIR}>
                          PRIVATE
                            "${EMsoftOO_SOURCE_DIR}/Source"                            
                    )

target_link_libraries(EMsoftOOHDFLib 
                          ${EMSOFT_hdf5LinkLibs} 
                          EMsoftOOLib 
)

set_target_properties (EMsoftOOHDFLib PROPERTIES
  LINKER_LANGUAGE Fortran
  FOLDER EMsoftOOPublic
  BUILD_RPATH "${EMsoftOO_OpenMP_LIB_DIR}"
)

EMsoftOO_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
    INSTALL(TARGETS EMsoftOOHDFLib
      EXPORT EMsoftOOHDFLibTargets
      COMPONENT Applications
      RUNTIME DESTINATION ${install_dir}
      LIBRARY DESTINATION ${lib_install_dir}
      ARCHIVE DESTINATION ${lib_install_dir}
      BUNDLE DESTINATION "."
  )
#endif()

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libEMsoftOOHDFLib.dylib")
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
#
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOLibTargetsConfigVersion.cmake"
  VERSION ${EMsoftOO_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT EMsoftOOHDFLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOHDFLibTargets.cmake"
    NAMESPACE EMsoftOO::
  )
#endif()

set(ConfigPackageLocation ${lib_install_dir}/cmake/EMsoftOO)

#if(BUILD_SHARED_LIBS)
  install(EXPORT EMsoftOOHDFLibTargets
    FILE
      EMsoftOOHDFLibTargets.cmake
    NAMESPACE
      EMsoftOO::
    DESTINATION
      ${ConfigPackageLocation}
  )
#endif()

