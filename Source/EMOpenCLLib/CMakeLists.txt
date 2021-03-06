include(${PROJECT_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake)

project(EMOpenCLLib)

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


set(EMOpenCLLib_SRCS
	${EMOpenCLLib_SOURCE_DIR}/mod_CLsupport.f90
  # ${EMOpenCLLib_SOURCE_DIR}/Indexingmod.f90

  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_MCOpenCL.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_SEMCLwrappers.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_DI.f90

# 	${EMOpenCLLib_SOURCE_DIR}/EMdymodHDF.f90
)

GetHDF5LinkLibraries(EMSOFTOO)

add_library(EMOpenCLLib ${LIB_TYPE} ${EMOpenCLLib_SRCS})

target_include_directories(EMOpenCLLib
                          PUBLIC
                            ${CLFortran_INSTALL}/include
                            $<BUILD_INTERFACE:${EMsoftOOLib_BINARY_DIR}>
                          PRIVATE
                            "${EMsoftOO_SOURCE_DIR}/Source"
                    )

target_link_libraries(EMOpenCLLib
                          ${OpenCL_LIBRARY}
                          clfortran
                          EMsoftOOLib
)


set_target_properties (EMOpenCLLib PROPERTIES
  LINKER_LANGUAGE Fortran
  FOLDER EMsoftOOPublic
  BUILD_RPATH "${EMsoftOO_OpenMP_LIB_DIR}"
)

EMsoftOO_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
    INSTALL(TARGETS EMOpenCLLib
      EXPORT EMOpenCLLibTargets
      COMPONENT Applications
      RUNTIME DESTINATION ${install_dir}
      LIBRARY DESTINATION ${lib_install_dir}
      ARCHIVE DESTINATION ${lib_install_dir}
      BUNDLE DESTINATION "."
  )
#endif()

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libEMOpenCLLib.dylib")
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
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMOpenCLLibTargetsConfigVersion.cmake"
  VERSION ${EMsoftOO_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT EMOpenCLLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMOpenCLLibTargets.cmake"
    NAMESPACE EMsoftOO::
  )
#endif()

set(ConfigPackageLocation ${lib_install_dir}/cmake/EMsoftOO)

#if(BUILD_SHARED_LIBS)
  install(EXPORT EMOpenCLLibTargets
    FILE
      EMOpenCLLibTargets.cmake
    NAMESPACE
      EMsoftOO::
    DESTINATION
      ${ConfigPackageLocation}
  )
#endif()
